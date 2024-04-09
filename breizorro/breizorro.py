import sys
import numpy
import shutil
import logging
import argparse
import os.path
import re
import numpy as np
import traceback
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from reproject import reproject_interp

import regions
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor, as_completed

from scipy.ndimage.morphology import binary_dilation, binary_fill_holes
from scipy.ndimage.measurements import label, find_objects
from scipy.ndimage import maximum_filter
import scipy.special
import scipy.ndimage

from bokeh.models import BoxEditTool, ColumnDataSource, FreehandDrawTool
from bokeh.plotting import figure, output_file, show
from bokeh.themes import built_in_themes
from bokeh.io import curdoc

def create_logger():
    """Create a console logger"""
    log = logging.getLogger(__name__)
    cfmt = logging.Formatter(('%(name)s - %(asctime)s %(levelname)s - %(message)s'))
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.INFO)
    console.setFormatter(cfmt)
    log.addHandler(console)
    return log

LOGGER = create_logger()
NCPU = 1

def get_image(fitsfile):
    """
    Reads FITS file, returns tuple of image_array, header
    """
    LOGGER.info(f"Reading {fitsfile} data")
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
        image = numpy.array(input_hdu.data[:, :])
    elif len(input_hdu.data.shape) == 3:
        image = numpy.array(input_hdu.data[0, :, :])
    else:
        image = numpy.array(input_hdu.data[0, 0, :, :])
    return image, input_hdu.header


def get_image_header(fitsfile):
    LOGGER.info(f"Reading {fitsfile} header")
    input_hdu = fits.open(fitsfile)[0]
    return input_hdu.header


def flush_fits(newimage, fitsfile, header=None):
    LOGGER.info(f"Writing {fitsfile}")
    f = fits.open(fitsfile, mode='update')
    input_hdu = f[0]
    if len(input_hdu.data.shape) == 2:
        input_hdu.data[:, :] = newimage
    elif len(input_hdu.data.shape) == 3:
        input_hdu.data[0, :, :] = newimage
    else:
        input_hdu.data[0, 0, :, :] = newimage
    if header:
        input_hdu.header = header
    f.flush()


def make_noise_map(restored_image, boxsize, quiet=False):
    # Cyril's magic minimum filter
    # Plundered from the depths of https://github.com/cyriltasse/DDFacet/blob/master/SkyModel/MakeMask.py
    if not quiet:
        LOGGER.info("Generating noise map")
    box = (boxsize, boxsize)
    n = boxsize**2.0
    x = numpy.linspace(-10, 10, 1000)
    f = 0.5 * (1.0 + scipy.special.erf(x / numpy.sqrt(2.0)))
    F = 1.0 - (1.0 - f)**n
    ratio = numpy.abs(numpy.interp(0.5, F, x))
    noise = -scipy.ndimage.filters.minimum_filter(restored_image, box) / ratio
    negative_mask = noise < 0.0
    noise[negative_mask] = 1.0e-10
    median_noise = numpy.median(noise)
    median_mask = noise < median_noise
    noise[median_mask] = median_noise
    if not quiet:
        LOGGER.info(f"Median noise value is {median_noise}")
    return noise

def process_residual_cube(cubefile, baseline_image, baseline_flux_threshold,
                          boxsize, threshold, reject_maskfile, reject_threshold, 
                          reject_dilate,
                          max_filter_size=8, 
                          maxplanes=None):

    LOGGER.info(f"Loading residual cube {cubefile}")
    basename = os.path.splitext(cubefile)[0]
    ff = fits.open(cubefile, memmap=True)

    cube = ff[0].data
    header = ff[0].header

    boost = baseline_flux_image = make_boost_cube = None
    if baseline_image:
        boost_cube = f"{basename}.boost.fits"
        if os.path.exists(boost_cube):
            LOGGER.info(f"Loading existing boost cube from {boost_cube}")
            boost = fits.open(boost_cube, memmap=True)[0].data
        else:
            rff = fits.open(baseline_image)[0]
            name = os.path.splitext(baseline_image)[0]
            if rff.shape[-2:] != cube.shape[-2:]:
                LOGGER.info(f"Reprojecting baseline image {baseline_image}")
                
                # form WCSs, dropping non-spatial axes
                target_wcs = WCS(header).dropaxis(2)
                rff_wcs = WCS(rff.header).dropaxis(3).dropaxis(2)
                rff_data = rff.data[0, 0, :, :]

                # Perform the reprojection
                baseline_flux_image, _ = reproject_interp((rff.data[0,0,...], rff_wcs), target_wcs, shape_out=cube.shape[1:])

                outname = f"{name}.reprojected.fits"
                LOGGER.info(f"Saving reprojected baseline image to {outname}")
                fits.PrimaryHDU(baseline_flux_image, header).writeto(outname, overwrite=True)
            else:
                LOGGER.info(f"Loading baseline image from {baseline_image}")
                baseline_flux_image = rff.data
            if max_filter_size:
                LOGGER.info("applying max filter")
                baseline_flux_image = maximum_filter(baseline_flux_image, 
                                                        size=(max_filter_size, max_filter_size), mode='constant')
                outname = f"{name}.baselineflux.fits"
                LOGGER.info(f"saving baseline flux image {outname}")
                fits.PrimaryHDU(baseline_flux_image, header).writeto(outname, overwrite=True)
            boost = np.zeros_like(cube, np.float32)
            make_boost_cube = True

    if reject_maskfile:
        rff = fits.open(reject_maskfile)[0]
        if rff.shape[-2:] != cube.shape[-2:]:
            LOGGER.info(f"Reprojecting external rejection mask from {reject_maskfile}")
            
            # form WCSs, dropping non-spatial axes
            target_wcs = WCS(header).dropaxis(2)
            rff_wcs = WCS(rff.header).dropaxis(3).dropaxis(2)
            rff_data = rff.data[0, 0, :, :]

            # Perform the reprojection
            reprojected_data, _ = reproject_interp((rff.data[0,0,...], rff_wcs), target_wcs, shape_out=cube.shape[1:])

            reject_external_mask = (reprojected_data != 0)

            name = os.path.splitext(reject_maskfile)[0]
            outname = f"{name}.reprojected.fits"
            LOGGER.info(f"Saving reprojected external rejection to {outname}")
            fits.PrimaryHDU(reject_external_mask.astype(np.int8), header).writeto(outname, overwrite=True)
        else:
            LOGGER.info(f"Loading external rejection mask from {reject_maskfile}")
            reject_external_mask = rff.data != 0

    max_sigma_bin = 200
    histogram_bins = np.arange(0, max_sigma_bin+1)
    histogram = np.zeros(len(histogram_bins)-1, int)
    max_sigma = 0

    LOGGER.info(f"Processing mean cube (reject threshold {reject_threshold}, dilate {reject_dilate})")

    nt, ny, nx = cube.shape

    avgcube = cube.mean(0)
    noise_image = make_noise_map(avgcube, boxsize)
    reject_mask = avgcube > reject_threshold * noise_image

    if reject_maskfile:
        reject_mask |= reject_external_mask

    if reject_dilate:
        r = np.arange(-reject_dilate, reject_dilate+1)
        struct = np.sqrt(r[:, np.newaxis]**2 + r[np.newaxis,:]**2) <= reject_dilate
        reject_mask = binary_dilation(input=reject_mask, structure=struct)

    LOGGER.info(f"Saving rejection mask")
    fits.PrimaryHDU(reject_mask.astype(np.int8), header).writeto(f"{basename}.reject.fits", overwrite=True)

    accept_mask = ~reject_mask
    del reject_mask

    detections = {}

    def process_plane(i):
        try:
            image_mask = accept_mask
            if boost is not None:
                # make the boost factor cube if needed
                if make_boost_cube:
                    boost[i] = cube[i] / baseline_flux_image
                    boost[i][baseline_flux_image < 2e-5] = 0
                # filter based on flux boost factor
                image_mask = image_mask & ((boost[i] == 0) | (boost[i] >= baseline_flux_threshold))
            # convert cube to sigmas
            cube[i] /= make_noise_map(cube[i], boxsize, quiet=True)
            img = cube[i]
            if i==0:
                fits.PrimaryHDU(img, header).writeto(f"snr.{i}.fits", overwrite=True)
                fits.PrimaryHDU(image_mask.astype(np.int8), header).writeto(f"mask.{i}.fits", overwrite=True)
            # compute histogram of sigma values
            hist, _ = np.histogram(img[image_mask], histogram_bins)
            maxsig = img[image_mask].max()
            # find stuff above threshold
            img[~image_mask] = 0
            image_mask = (img >= threshold)
            yy, xx = np.nonzero(image_mask)
            if not len(xx):
                return i, [], [], [], [], maxsig, hist
            # LOGGER.info(f"Plane {i} max sigma is {maxsig}")
            return i, xx, yy, \
                    [img[y, x] for x, y in zip(xx, yy)], \
                    [boost[i, y, x] for x, y in zip(xx, yy)] if boost is not None else [1]*len(xx), \
                    maxsig, hist
        except Exception as exc:
            print(f"Error processing plane {i}")
            traceback.print_exc()

    with ThreadPoolExecutor(max_workers=NCPU) as pool:
        if not maxplanes:
            maxplanes = cube.shape[0]
        futures = [pool.submit(process_plane, i) for i in range(maxplanes)] 
        for f in as_completed(futures):
            i, xx, yy, sigmas, boosts, maxsig, hist = f.result()
            for x, y, sigma, b in zip(xx, yy, sigmas, boosts):
                detections.setdefault((x,y),[]).append((i, sigma, b))
            if len(xx) or not i%100:
                LOGGER.info(f"Plane {i} had {len(xx)} detection(s)")
            max_sigma = max(max_sigma, maxsig)
            histogram[...] += hist

    LOGGER.info("saving sigma cube")
    fits.PrimaryHDU(cube, header).writeto(f"{basename}.sigma.fits", overwrite=True)

    LOGGER.info("saving detections")
    sorted_detections = sorted([len(tt), x, y, sorted(tt)] for (x,y), tt in detections.items())[::-1]
    detection_mask = np.zeros_like(cube[0], np.int16)
    for i, (n, x, y, tsb) in enumerate(sorted_detections):
        dets = [f"{t}:{s:.1f}" for t, s, b in tsb]
        maxsig = max([s for t, s, b in tsb])
        print(tsb)
        maxboost = max([b for t, s, b in tsb])
        LOGGER.info(f"Detection {i} at {x},{y} with max sigma {maxsig:.1f} and max boost {maxboost:.2f}: {n} repeats")
        xslice = slice(max(0, x-2), min(nx, x+2))
        yslice = slice(max(0, y-2), min(ny, y+2))
        detection_mask[yslice,xslice] = tsb[0][0]

    fits.PrimaryHDU(detection_mask, header).writeto(f"{basename}.detection.fits", overwrite=True)
    
    print(f"Max sigma is {max_sigma}. Pixels per sigma bin:")
    maxbin = int(max_sigma) + 2
    for bin, num in enumerate(histogram[:maxbin]):
        print(f"{histogram_bins[bin]}-{histogram_bins[bin+1]}: {num}") 

    if make_boost_cube:
        LOGGER.info(f"saving boost cube {boost_cube}")
        fits.PrimaryHDU(boost, header).writeto(boost_cube, overwrite=True)


def resolve_island(isl_spec, mask_image, wcs, ignore_missing=False):
    if re.match("^\d+$", isl_spec):
        return int(isl_spec)
    elif ',' not in isl_spec:
        raise ValueError(f"invalid island specification: {isl_spec}")
    c = SkyCoord(*isl_spec.split(',', 1))
    x, y = wcs.world_to_pixel(c)
    x = round(float(x))
    y = round(float(y))
    value = mask_image[y, x]
    LOGGER.info(f"coordinates {c} correspond to pixel {x}, {y} with value {value}")
    if not value:
        if ignore_missing:
            LOGGER.warning("no island at specified coordinates, ignoring")
        else:
            raise ValueError(f"coordinates {c} do not select a valid island")
    return value

def add_regions(mask_image, regs, wcs):
    for reg in regs:
        if hasattr(reg, 'to_pixel'):
            reg = reg.to_pixel(wcs)
        mask_image += reg.to_mask().to_image(mask_image.shape)

def remove_regions(mask_image, regs, wcs):
    for reg in regs:
        if hasattr(reg, 'to_pixel'):
            reg = reg.to_pixel(wcs)
        mask_image[reg.to_mask().to_image(mask_image.shape) != 0] = 0


def main():
    LOGGER.info("Welcome to breizorro")
    # Get version
    from pkg_resources import get_distribution
    _version = get_distribution('breizorro').version
    LOGGER.info(f"Version: {_version}")

    LOGGER.info("Usage: breizorro --help")
    parser = ArgumentParser(description='breizorro [options] --restored-image restored_image')
    parser.add_argument('-r', '--restored-image', dest="imagename", metavar="IMAGE", 
                        help="Restored image file from which to build mask")
    parser.add_argument('-m', '--mask-image', dest="maskname", metavar="MASK",
                        help="Input mask file(s). Either --restored-image or --mask-image must be specfied.")
    parser.add_argument('-t', '--threshold', dest='threshold', default=6.5,
                        help='Sigma threshold for masking (default = 6.5)')
    parser.add_argument('-b', '--boxsize', dest='boxsize', default=50,
                        help='Box size over which to compute stats (default = 50)')
    parser.add_argument('--savenoise', dest='savenoise', action='store_true', default=False,
                        help='Enable to export noise image as FITS file (default=do not save noise image)')

    parser.add_argument('--merge', dest='merge', metavar="MASK(s)|REG(s)", nargs='+',
                        help='Merge in one or more masks or region files')
    parser.add_argument('--subtract', dest='subtract', metavar="MASK(s)|REG(s)", nargs='+',
                        help='Subract one or more masks or region files')

    parser.add_argument('--number-islands', dest='islands', action='store_true', default=False,
                        help='Number the islands detected (default=do not number islands)')
    parser.add_argument('--remove-islands', dest='remove_isl', metavar='N|COORD', type=str, nargs='+',
                         help='List of islands to remove from input mask. e.g. --remove-islands 1 18 20 20h10m13s,14d15m20s')
    parser.add_argument('--ignore-missing-islands', action='store_true', 
                         help='If an island specified by coordinates does not exist, do not throw an error')
    parser.add_argument('--extract-islands', dest='extract_isl', metavar='N|COORD', type=str, nargs='+',
                         help='List of islands to extract from input mask. e.g. --extract-islands 1 18 20 20h10m13s,14d15m20s')
    parser.add_argument('--make-binary', action="store_true",
                         help='Replace all island numbers with 1')
    parser.add_argument('--invert', action="store_true",
                         help='Invert the mask')
    
    parser.add_argument('-c', '--cube', type=str, 
                        help="Process residual cube instead of main image")
    parser.add_argument('--max-cube-planes', dest='max_cube_planes', metavar='N', type=int, 
                        help='Proecess only the first N cube planes. Useful for debugging.')
    parser.add_argument('--reject-threshold', dest='reject_threshold', type=float, default=5,
                        help='Rejection threshold for mean cube (default %(default)s)')
    parser.add_argument('--reject-dilate', dest='reject_dilate', type=int, default=3,
                        help='Dilate rejection mask by this value (default %(default)s)')
    parser.add_argument('--reject-boxsize', dest='reject_boxsize', type=int, default=50,
                        help='Box size for rejection mask')
    parser.add_argument('--reject-mask', dest='reject_mask', type=str,
                        help='Additional rejection mask')
    parser.add_argument('--max-filter-size', type=int, default=8,
                        help='Dilate baseline flux image by this max-filter (default %(default)s)')
    parser.add_argument('--baseline-flux-threshold', metavar="FACTOR", dest='baseline_flux_threshold', type=float,
                        help='Reject if flux excursion is below baseline restored image times this factor')

    parser.add_argument('--dilate', dest='dilate', metavar="R", type=int, default=0,
                        help='Apply dilation with a radius of R pixels')
    parser.add_argument('--fill-holes', dest='fill_holes', action='store_true', 
                        help='Fill holes (i.e. entirely closed regions) in mask')

    parser.add_argument('--sum-peak', dest='sum_peak', default=None,
                        help='Sum to peak ratio of flux islands to mask in original image.'
                             'e.g. --sum-peak 100 will mask everything with a ratio above 100')

    parser.add_argument('-o', '--outfile', dest='outfile', default='',
                        help='Suffix for mask image (default based on input name')

    parser.add_argument('--gui', dest='gui', action='store_true', default=False,
                         help='Open mask in gui.')
    
    parser.add_argument('-j', '--ncpu', dest='ncpu', default=1, type=int, 
                         help='Number of CPUs to use, for operations that support parallelism.')

    args = parser.parse_args()
    threshold = float(args.threshold)
    boxsize = int(args.boxsize)
    dilate = int(args.dilate)
    savenoise = args.savenoise
    outfile = args.outfile
    global NCPU
    NCPU = args.ncpu

    if args.imagename and args.maskname:
        parser.error("Either --restored-image or --mask-image must be specified, but not both")

    # define input file, and get its name and extension
    input_file = args.imagename or args.maskname or args.cube
#    name, ext = os.path.split(input_file)
    name = '.'.join(input_file.split('.')[:-1])
    ext = input_file.split('.')[-1]


    if args.cube:
        mask_image = process_residual_cube(args.cube,
                                           boxsize=boxsize, threshold=threshold, 
                                           baseline_image=args.imagename,
                                           baseline_flux_threshold=args.baseline_flux_threshold,
                                           reject_maskfile=args.reject_mask,
                                           reject_threshold=args.reject_threshold, 
                                           reject_dilate=args.reject_dilate,
                                           max_filter_size=args.max_filter_size,
                                           maxplanes=args.max_cube_planes)
        sys.exit(0)

    # first, load or generate mask

    if args.imagename:
        input_image, input_header = get_image(input_file)
        LOGGER.info(f"Generating mask using threshold {threshold}")

        noise_image = make_noise_map(input_image, boxsize)
        if savenoise:
            noise_fits = f"{name}.noise.fits"
            shutil.copyfile(input_file, noise_fits)
            flush_fits(noise_image, noise_fits)

        mask_image = input_image > threshold * noise_image

        mask_image[:, -1]=0
        mask_image[:, 0]=0
        mask_image[0, :]=0
        mask_image[-1, :]=0

        mask_header = input_header
        mask_header['BUNIT'] = 'mask'

        out_mask_fits = args.outfile or f"{name}.mask.fits"

    elif args.maskname:
        mask_image, mask_header = get_image(args.maskname)
        LOGGER.info(f"Input mask loaded")

        out_mask_fits = args.outfile or f"{name}.out.{ext}"
    else:
        parser.error("Either --restored-image or --mask-image must be specified")
        sys.exit(1)

    wcs = WCS(mask_header)
    while len(wcs.array_shape) > 2:
        wcs = wcs.dropaxis(len(wcs.array_shape) - 1)

    # next, merge and/or subtract
    def load_fits_or_region(filename):
        ff = regs = None
        # read as FITS or region
        try:
            ff = get_image(filename)
        except OSError:
            try:
                regs = regions.Regions.read(filename)
            except:
                LOGGER.error(f"{merge} is neither a FITS file not a regions file")
                raise
        return ff, regs


    if args.merge:
        for merge in args.merge:
            ff, regs = load_fits_or_region(merge)
            if ff:
                LOGGER.info(f"Treating {merge} as a FITS mask")
                mask_image += ff[0]
                LOGGER.info("Merged into mask")
            else:
                LOGGER.info(f"Merging in {len(regs)} regions from {merge}")
                add_regions(mask_image, regs, wcs)
        mask_image = mask_image != 0
        mask_header['BUNIT'] = 'mask'

    if args.subtract:
        for subtract in args.subtract:
            ff, regs = load_fits_or_region(subtract)
            if fits:
                LOGGER.info(f"treating {subtract} as a FITS mask")
                mask_image[ff[0] != 0] = 0
                LOGGER.info("Subtracted from mask")
            else:
                LOGGER.info(f"Subtracting {len(regs)} regions from {subtract}")
                remove_regions(mask_image, regs, wcs)

    if args.islands:
        LOGGER.info(f"(Re)numbering islands")
        mask_image = mask_image != 0
        # mask_image = mask_image.byteswap().newbyteorder()
        mask_image, num_features = label(mask_image)
        mask_header['BUNIT'] = 'Source_ID'
        LOGGER.info(f"Number of islands: {num_features}")
    
    if args.remove_isl:
        LOGGER.info(f"Removing islands: {args.remove_isl}")
        for isl_spec in args.remove_isl:
            isl = resolve_island(isl_spec, mask_image, wcs, ignore_missing=args.ignore_missing_islands)
            if isl != None:
                mask_image[mask_image == isl] = 0

    if args.extract_isl:
        LOGGER.info(f"Extracting islands: {args.extract_isl}")
        new_mask_image = np.zeros_like(mask_image)
        for isl_spec in args.extract_isl:
            isl = resolve_island(isl_spec, mask_image, wcs)
            new_mask_image[mask_image == isl] = isl
        mask_image = new_mask_image

    if args.make_binary:
        LOGGER.info(f"Converting mask to binary")
        mask_image = mask_image!=0
        mask_header['BUNIT'] = 'mask'

    if args.invert:
        LOGGER.info(f"Inverting mask")
        mask_image = mask_image==0

    if args.dilate:
        LOGGER.info(f"Dilating mask using a ball of R={args.dilate}pix")
        R = args.dilate
        r = np.arange(-R, R+1)
        struct = np.sqrt(r[:, np.newaxis]**2 + r[np.newaxis,:]**2) <= R
        mask_image = binary_dilation(input=mask_image, structure=struct)
        
    if args.fill_holes:
        LOGGER.info(f"Filling closed regions")
        mask_image = binary_fill_holes(mask_image)

    if args.sum_peak:
        # This mainly to produce an image that mask out super extended sources (via sum-to-peak flux ratio)
        # This is useful to allow source finder to detect mainly point-like sources for cross-matching purposes only.
        LOGGER.info(f"Including only flux islands with a sum-peak ratio below: {args.sum_peak}")
        extended_islands = []
        mask_image_label, num_features = label(mask_image)
        island_objects = find_objects(mask_image_label.astype(int))
        for island in island_objects:
            isl_sum = (input_image[island] * mask_image[island]).sum()
            isl_peak = (input_image[island] * mask_image[island]).max()
            isl_sum_peak = isl_sum / isl_peak
            if isl_sum_peak > float(args.sum_peak):
                extended_islands.append(island)
        new_mask_image = np.zeros_like(mask_image)
        new_mask_image = new_mask_image == 0
        for ext_isl in extended_islands:
            isl_slice = mask_image[ext_isl] == 0
            new_mask_image[ext_isl] = isl_slice
        mask_header['BUNIT'] = 'Jy/beam'
        mask_image = input_image * new_mask_image
        LOGGER.info(f"Number of extended islands found: {len(extended_islands)}")

    if args.gui:
        curdoc().theme = 'caliber'

        LOGGER.info("Loading Gui ...")
        d = mask_image
        p = figure(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
        p.x_range.range_padding = p.y_range.range_padding = 0
        p.title.text = out_mask_fits

        # must give a vector of image data for image parameter
        p.image(image=[d], x=0, y=0, dw=10, dh=10, palette="Greys256", level="image")
        p.grid.grid_line_width = 0.5
        src1 = ColumnDataSource({'x': [], 'y': [], 'width': [], 'height': [], 'alpha': []})
        src2 = ColumnDataSource({'xs': [], 'ys': [], 'alpha': []})
        renderer1 = p.rect('x', 'y', 'width', 'height', source=src1, alpha='alpha')
        renderer2 = p.multi_line('xs', 'ys', source=src2, alpha='alpha')
        draw_tool1 = BoxEditTool(renderers=[renderer1], empty_value=1)
        draw_tool2 = FreehandDrawTool(renderers=[renderer2])
        p.add_tools(draw_tool1)
        p.add_tools(draw_tool2)
        p.toolbar.active_drag = draw_tool1
        output_file("breizorro.html", title="Mask Editor")
        show(p)

        LOGGER.info(f"Enforcing that mask to binary")
        mask_image = mask_image!=0
        mask_header['BUNIT'] = 'mask'

    shutil.copyfile(input_file, out_mask_fits)  # to provide a template
    flush_fits(mask_image, out_mask_fits, mask_header)
    LOGGER.info("Done")
