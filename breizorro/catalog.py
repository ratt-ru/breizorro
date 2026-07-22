import logging
from multiprocessing import Pool
from operator import itemgetter

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import wcs_to_celestial_frame
from photutils import centroids
from regions import PolygonSkyRegion
from tqdm import tqdm

from breizorro.utils import calculate_beam_area, get_source_size

logger = logging.getLogger(__name__)


def centroid_gaussian(data):
    """
    Gaussian-weighted centroid fitting.
    Can hang on poorly conditioned data.
    """
    try:
        return centroids.centroid_2dg(data)
    except Exception as e:
        # Fallback to center of mass if Gaussian fitting fails
        logger.warning(f"Gaussian centroid fitting failed: {e}.\nFalling back to center of mass.")
        return centroids.centroid_com(data)


def centroid_com(data):
    """
    Center of mass (fastest, most robust).
    Recommended as default.
    """
    return centroids.centroid_com(data)


def centroid_moments(data):
    """
    Image moments-based centroid (fast and robust).
    Uses weighted center calculation.
    """
    return centroids.centroid_1dg(data)


def centroid_windowed(data, window_size=None):
    """
    Windowed centroid around peak (balances speed and accuracy).
    Finds peak, then calculates centroid in local window.
    """
    # Find peak location
    peak_idx = np.unravel_index(np.argmax(data), data.shape)

    # Set window size based on data size if not provided
    if window_size is None:
        window_size = min(data.shape) // 4

    # Extract window around peak
    y_start = max(0, peak_idx[0] - window_size)
    y_end = min(data.shape[0], peak_idx[0] + window_size)
    x_start = max(0, peak_idx[1] - window_size)
    x_end = min(data.shape[1], peak_idx[1] + window_size)

    windowed_data = data[y_start:y_end, x_start:x_end]

    # Calculate COM on window
    com_y, com_x = centroids.centroid_com(windowed_data)

    # Adjust back to original coordinates
    return (com_y + y_start, com_x + x_start)


def get_centroid_method(method_name="centroid"):
    """
    Return the centroid fitting method function.

    Available methods:
    - gaussian: Gaussian-weighted (can hang, has COM fallback)
    - centroid: Center of mass (fast, recommended)
    - moments: Image moments (fast, robust)
    - windowed: Windowed COM around peak (fast, moderate accuracy)
    """
    methods = {
        "gaussian": centroid_gaussian,
        "centroid": centroid_com,
        "com": centroid_com,
        "moments": centroid_moments,
        "windowed": centroid_windowed,
    }

    method = methods.get(method_name.lower(), centroid_com)
    return method


def estimate_position_uncertainty(peak_flux, noise_out, mean_beam, dec_deg, source_size=None):
    """Estimate RA/DEC position errors from beam size and peak SNR.

    The uncertainty is estimated as FWHM / (2 * SNR), with the RA error
    corrected for declination. For resolved sources we use the larger of the
    fitted source size and the restoring beam as a conservative width.
    """
    safe_noise = max(abs(noise_out), np.finfo(float).tiny)
    snr = max(abs(peak_flux) / safe_noise, 1.0)

    if source_size and len(source_size) >= 2 and (source_size[0] > 0.0 or source_size[1] > 0.0):
        major_fwhm = max(float(source_size[0]), float(mean_beam))
        minor_fwhm = max(float(source_size[1]), float(mean_beam))
    else:
        major_fwhm = float(mean_beam)
        minor_fwhm = float(mean_beam)

    cos_dec = max(np.cos(np.deg2rad(dec_deg)), np.finfo(float).eps)
    ra_err_deg = (major_fwhm / (2.0 * snr)) / 3600.0 / cos_dec
    dec_err_deg = (minor_fwhm / (2.0 * snr)) / 3600.0
    return ra_err_deg, dec_err_deg


def format_scientific(value):
    return f"{value:.4e}"


def process_contour(contour, image_data, fitsinfo, noise_out, source_fitting="centroid"):
    use_max = 0
    ra = -np.inf
    catalog_out = ""
    pix_size = fitsinfo["ddec"] * 3600.0
    bmaj, bmin, _ = np.array(fitsinfo["b_size"]) * 3600.0
    mean_beam = 0.5 * (bmaj + bmin)
    pix_beam = calculate_beam_area(bmaj, bmin, pix_size)
    wcs = fitsinfo["wcs"]
    while len(wcs.array_shape) > 2:
        wcs = wcs.dropaxis(len(wcs.array_shape) - 1)

    contour_sky = wcs.pixel_to_world(contour[:, 1], contour[:, 0])
    polygon_region = PolygonSkyRegion(vertices=contour_sky)
    pix_region = polygon_region.to_pixel(wcs)
    mask = pix_region.to_mask().to_image(image_data.shape[-2:])
    # Calculate the number of pixels in the source region (where mask > 0)
    source_area_pix = np.sum(mask > 0)  # Count of pixels in the masked source region
    # Now calculate the number of beams covering the source
    source_beams = source_area_pix / pix_beam  # Number of beams covering the source
    try:
        data = mask * image_data
        nndata = data  # np.flip(data, axis=0)
        # nndata = nndata[~np.isnan(nndata)]
        total_flux = np.sum(nndata[nndata != -0.0]) / pix_beam
        peak_flux = nndata.max()
    except (ValueError, ZeroDivisionError):
        total_flux = 0.0
        peak_flux = 0.0
    if not np.isfinite(total_flux) or total_flux == 0:
        return (ra, catalog_out, use_max)

    if total_flux:
        total_peak_ratio = np.abs((total_flux - peak_flux) / total_flux)
        # Flux density error estimation
        ten_pc_error = 0.1 * total_flux  # a 10% error term as an additional conservative estimate
        beam_error = np.sqrt(source_beams) * noise_out
        flux_density_error = np.sqrt(ten_pc_error**2 + beam_error**2)  # combined error
        peak_error = np.sqrt((0.1 * peak_flux) ** 2 + beam_error**2)

        # Calculate weighted centroid using selected method
        centroid_method = get_centroid_method(source_fitting)
        try:
            _centroids = centroid_method(data)
        except Exception as e:
            # Fallback to center of mass if any fitting method fails
            logger.warning(f"Centroid fitting with {source_fitting} failed: {e}. Falling back to center of mass.")
            _centroids = centroid_com(data)

        centroid_x, centroid_y = _centroids
        world_x, world_y = wcs.all_pix2world(centroid_x, centroid_y, 0)
        # The WCS's native frame isn't necessarily equatorial
        # Converting through SkyCoord ensures ra/dec are genuine ICRS coordinates
        native_frame = wcs_to_celestial_frame(wcs)
        world_coord = SkyCoord(world_x, world_y, unit="deg", frame=native_frame).icrs
        ra, dec = world_coord.ra.deg, world_coord.dec.deg
        # Ensure RA is positive
        if ra < 0:
            ra += 360
        source_size = get_source_size(contour, pix_size, mean_beam, image_data, total_peak_ratio, _centroids)
        # For unresolved (point) sources, total_flux should equal peak_flux
        if source_size[0] == 0.0 and source_size[1] == 0.0:
            total_flux = peak_flux
            flux_density_error = peak_error
        ra_error, dec_error = estimate_position_uncertainty(peak_flux, noise_out, mean_beam, dec, source_size)
        source = (
            f"{ra:.4f}",
            f"{dec:.4f}",
            format_scientific(ra_error),
            format_scientific(dec_error),
            format_scientific(total_flux),
            format_scientific(flux_density_error),
            format_scientific(peak_flux),
            format_scientific(peak_error),
        ) + source_size
        catalog_out = " ".join(str(src_prop) for src_prop in source)
    return (ra, catalog_out, use_max)


def _worker_task(args):
    """Module-level worker function for multiprocessing (must be at module level to be picklable)."""
    # legacy signature kept for compatibility but prefer new single-arg use
    return process_contour(args, _shared_image_data, _shared_fitsinfo, _shared_noise_out, _shared_source_fitting)


def init_worker(image_data, fitsinfo, noise_out, source_fitting):
    """Initializer for worker processes to set large shared objects once per process."""
    global _shared_image_data, _shared_fitsinfo, _shared_noise_out, _shared_source_fitting
    _shared_image_data = image_data
    _shared_fitsinfo = fitsinfo
    _shared_noise_out = noise_out
    _shared_source_fitting = source_fitting


def multiprocess_contours(contours, image_data, fitsinfo, noise_out, ncpu=None, source_fitting="centroid"):
    """Process contours in parallel with progress bar."""
    # Determine number of CPUs
    if not ncpu:
        try:
            import multiprocessing

            ncpu = multiprocessing.cpu_count()
        except (RuntimeError, NotImplementedError):
            ncpu = 1

    # Build task list (only send contour objects; big arrays go to worker initializer)
    tasks = [contour for contour in contours if len(contour) > 2]

    # Process with progress bar. Use initializer to set shared large objects once per worker.
    source_list = []
    with Pool(
        processes=ncpu, initializer=init_worker, initargs=(image_data, fitsinfo, noise_out, source_fitting)
    ) as pool:
        for catalog_out in tqdm(pool.imap_unordered(_worker_task, tasks), total=len(tasks), desc="Finding sources"):
            if catalog_out[0] > -np.inf:
                source_list.append(catalog_out)

    ra_sorted_list = sorted(source_list, key=itemgetter(0))
    return ra_sorted_list
