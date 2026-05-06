import logging
from multiprocessing import Process, Queue
from operator import itemgetter

import numpy as np
from photutils import centroids
from regions import PolygonSkyRegion

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


def process_contour(contour, image_data, fitsinfo, noise_out, source_fitting="centroid"):
    use_max = 0
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
        peak_flux = 0.0
    if total_flux:
        total_peak_ratio = np.abs((total_flux - peak_flux) / total_flux)
        # Flux density error estimation
        ten_pc_error = 0.1 * total_flux  # a 10% error term as an additional conservative estimate
        beam_error = np.sqrt(source_beams) * noise_out
        flux_density_error = np.sqrt(ten_pc_error**2 + beam_error**2)  # combined error

        # Calculate weighted centroid using selected method
        centroid_method = get_centroid_method(source_fitting)
        try:
            _centroids = centroid_method(data)
        except Exception as e:
            # Fallback to center of mass if any fitting method fails
            logger.warning(
                f"Centroid fitting with {source_fitting} failed: {e}. Falling back to center of mass."
            )
            _centroids = centroid_com(data)

        centroid_x, centroid_y = _centroids
        ra, dec = wcs.all_pix2world(centroid_x, centroid_y, 0)
        # Ensure RA is positive
        if ra < 0:
            ra += 360
        source_flux = (round(total_flux, 5), round(flux_density_error, 5))
        source_size = get_source_size(
            contour, pix_size, mean_beam, image_data, total_peak_ratio, _centroids
        )
        # source_pos = format_source_coordinates(ra, dec)
        source = (ra, dec) + source_flux + source_size
        catalog_out = " ".join(str(src_prop) for src_prop in source)
    else:
        # Dummy source to be eliminated

        catalog_out = ""
    return (ra, catalog_out, use_max)


def multiprocess_contours(
    contours, image_data, fitsinfo, noise_out, ncpu=None, source_fitting="centroid"
):

    def contour_worker(input, output):
        for func, args in iter(input.get, "STOP"):
            result = func(*args)
            output.put(result)

    source_list = []
    # Start worker processes
    if not ncpu:
        try:
            import multiprocessing

            ncpu = multiprocessing.cpu_count()
        except (RuntimeError, NotImplementedError):
            pass
    TASKS = []
    for i in range(len(contours)):
        contour = contours[i]
        if len(contour) > 2:
            x = []
            y = []
            for j in range(len(contour)):
                x.append(contour[j][0])
                y.append(contour[j][1])
            TASKS.append(
                (process_contour, (contour, image_data, fitsinfo, noise_out, source_fitting))
            )
    task_queue = Queue()
    done_queue = Queue()
    # Submit tasks
    for task in TASKS:
        task_queue.put(task)
    for i in range(ncpu):
        Process(target=contour_worker, args=(task_queue, done_queue)).start()

    num_max = 0
    # Get the results from parallel processing
    for i in range(len(TASKS)):
        catalog_out = done_queue.get(timeout=1800)
        if catalog_out[0] > -np.inf:
            source_list.append(catalog_out)
            num_max += catalog_out[2]
    # Tell child processes to stop
    for i in range(ncpu):
        task_queue.put("STOP")

    ra_sorted_list = sorted(source_list, key=itemgetter(0))

    return ra_sorted_list
