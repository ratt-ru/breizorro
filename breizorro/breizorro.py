import sys
import numpy
import shutil
import logging
import argparse
import os.path
import re
import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord

import regions
from regions import PixCoord
from regions import PolygonSkyRegion, PolygonPixelRegion

from argparse import ArgumentParser

from scipy.ndimage.morphology import binary_dilation, binary_erosion, binary_fill_holes
from scipy.ndimage.measurements import label, find_objects
import scipy.special
import scipy.ndimage

from skimage.draw import polygon as skimage_polygon
from skimage.measure import find_contours

from breizorro.utils import get_source_size, format_source_coordinates, deg2ra, deg2dec
from breizorro.utils import get_image_data, fitsInfo, calculate_beam_area


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


def flush_fits(newimage, fitsfile, header=None):
    LOGGER.info(f"Writing {fitsfile}")
    with fits.open(fitsfile, mode='update') as f:
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


def make_noise_map(restored_image, boxsize):
    # Cyril's magic minimum filter
    # Plundered from the depths of https://github.com/cyriltasse/DDFacet/blob/master/SkyModel/MakeMask.py
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
    LOGGER.info(f"Median noise value is {median_noise}")
    return noise


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
    parser.add_argument('--minimum-size', dest='minsize', type=int,
                        help='Remove islands that have areas fewer than or equal to the specified number of pixels')
    parser.add_argument('--make-binary', action="store_true",
                         help='Replace all island numbers with 1')
    parser.add_argument('--invert', action="store_true",
                         help='Invert the mask')
    
    parser.add_argument('--dilate', dest='dilate', metavar="R", type=int, default=0,
                        help='Apply dilation with a radius of R pixels')
    parser.add_argument('--erode', dest='erode', metavar="N", type=int, default=0,
                        help='Apply N iterations of erosion')
    parser.add_argument('--fill-holes', dest='fill_holes', action='store_true', 
                        help='Fill holes (i.e. entirely closed regions) in mask')

    parser.add_argument('--sum-peak', dest='sum_peak', default=None,
                        help='Sum to peak ratio of flux islands to mask in original image.'
                             'e.g. --sum-peak 100 will mask everything with a ratio above 100')
    parser.add_argument('-ncpu', '--ncpu', dest='ncpu', default=None, type=int,
                        help='Number of processors to use for cataloging.')
    parser.add_argument('-beam', '--beam-size', dest='beam', default=None,
                        help='Average beam size in arcesc incase beam info is missing in image header.')

    parser.add_argument('-o', '--outfile', dest='outfile', default='',
                        help='Suffix for mask image (default based on input name')
    parser.add_argument('--save-catalog', dest='outcatalog', default='',
                        help='Generate catalog based on region mask')
    parser.add_argument('--save-regions', dest='outregion', default='',
                         help='Generate polygon regions from the mask')


    parser.add_argument('--gui', dest='gui', action='store_true', default=False,
                         help='Open mask in gui.')
    args = parser.parse_args()
    threshold = float(args.threshold)
    boxsize = int(args.boxsize)
    dilate = int(args.dilate)
    savenoise = args.savenoise
    outfile = args.outfile

    if args.imagename and args.maskname:
        parser.error("Either --restored-image or --mask-image must be specified, but not both")
    elif not args.imagename and not args.maskname:
        parser.error("Either --restored-image or --mask-image must be specified")

    # define input file, and get its name and extension
    input_file = args.imagename or args.maskname
    name = '.'.join(input_file.split('.')[:-1])
    ext = input_file.split('.')[-1]

    # first, load or generate mask
    if args.imagename:
        input_image, input_header = get_image_data(input_file)
        LOGGER.info(f"Generating mask using threshold {threshold}")

        noise_image = make_noise_map(input_image, boxsize)
        if savenoise:
            noise_fits = f"{name}.noise.fits"
            shutil.copyfile(input_file, noise_fits)
            flush_fits(noise_image, noise_fits)

        mask_image = (input_image > threshold * noise_image).astype('float64')

        mask_image[:, -1]=0
        mask_image[:, 0]=0
        mask_image[0, :]=0
        mask_image[-1, :]=0

        mask_header = input_header
        mask_header['BUNIT'] = 'mask'

        out_mask_fits = args.outfile or f"{name}.mask.fits"

    elif args.maskname:
        mask_image, mask_header = get_image_data(args.maskname)
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
        fits = regs = None
        # read as FITS or region
        try:
            fits = get_image_data(filename)
        except OSError:
            try:
                regs = regions.Regions.read(filename)
            except:
                LOGGER.error(f"{merge} is neither a FITS file not a regions file")
                raise
        return fits, regs


    if args.merge:
        for merge in args.merge:
            fits, regs = load_fits_or_region(merge)
            if fits:
                LOGGER.info(f"Treating {merge} as a FITS mask")
                mask_image += fits[0]
                LOGGER.info("Merged into mask")
            else:
                LOGGER.info(f"Merging in {len(regs)} regions from {merge}")
                add_regions(mask_image, regs, wcs)
        mask_image = mask_image != 0
        mask_header['BUNIT'] = 'mask'

    if args.subtract:
        for subtract in args.subtract:
            fits, regs = load_fits_or_region(subtract)
            if fits:
                LOGGER.info(f"treating {subtract} as a FITS mask")
                mask_image[fits[0] != 0] = 0
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

    if args.minsize:
        LOGGER.info(f"Removing islands that occupy fewer than or equal to {args.minsize} pixels")
        mask_image = mask_image != 0
        island_labels, num_features = label(mask_image)
        island_areas = numpy.array(scipy.ndimage.sum(mask_image,island_labels, numpy.arange(island_labels.max()+1)))
        min_mask = island_areas >= args.minsize
        mask_image = min_mask[island_labels.ravel()].reshape(island_labels.shape)

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

    if args.erode:
        LOGGER.info(f"Applying {args.erode} iteration(s) of erosion")
        N = args.erode
        mask_image = binary_erosion(input=mask_image, iterations=N)
        
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
        shutil.copyfile(input_file, out_mask_fits)  # to provide a template
        flush_fits(mask_image, out_mask_fits, mask_header)
        LOGGER.info("Done")
        sys.exit(1)

    if args.outcatalog or args.outregion:
        contours = find_contours(mask_image, 0.5)
        polygon_regions = []
        for contour in contours:
            # Convert the contour points to pixel coordinates
            contour_pixels = contour
            # Convert the pixel coordinates to Sky coordinates
            contour_sky = wcs.pixel_to_world(contour_pixels[:, 1], contour_pixels[:, 0])
            # Create a Polygon region from the Sky coordinates
            polygon_region = PolygonSkyRegion(vertices=contour_sky, meta={'label': 'Region'})
            # Add the polygon region to the list
            polygon_regions.append(polygon_region)
        LOGGER.info(f"Number of regions found: {len(polygon_regions)}")
        if args.outregion:
            regions.Regions(polygon_regions).write(args.outregion, format='ds9')
            LOGGER.info(f"Saving regions in {args.outregion}")

    if args.outcatalog and args.imagename:
        try:
            import warnings
            # Suppress FittingWarnings from Astropy
            # WARNING: The fit may be unsuccessful; check fit_info['message'] for more information. [astropy.modeling.fitting]
            warnings.resetwarnings()
            warnings.filterwarnings('ignore', category=UserWarning, append=True)
            from breizorro.catalog import multiprocess_contours
        except ModuleNotFoundError:
            LOGGER.error("Running breizorro source detector requires optional dependencies, please re-install with: pip install breizorro[catalog]")
            raise('Missing cataloguing dependencies')
        source_list = []
        image_data, hdu_header = get_image_data(args.imagename)
        fitsinfo = fitsInfo(args.imagename)
        mean_beam = None # Use beam info from the image header by default
        if args.beam:
            mean_beam = float(args.beam)
        if mean_beam:
            LOGGER.info(f'Using user provided size: {mean_beam}')
        elif fitsinfo['b_size']:
            bmaj,bmin,_ = np.array(fitsinfo['b_size']) * 3600.0
            mean_beam = 0.5 * (bmaj + bmin)
        else:
            raise('No beam information found. Specify mean beam in arcsec: --beam-size 6.5')

        noise = np.median(noise_image)
        f = open(args.outcatalog, 'w')
        catalog_out = f'# processing fits image: {args.imagename}  \n'
        f.write(catalog_out)
        image_dimensions = fitsinfo['naxis']
        pixel_size = fitsinfo['ddec'] * 3600.0
        catalog_out = f'# mean beam size (arcsec): {round(mean_beam,2)} \n' 
        f.write(catalog_out)
        catalog_out = f'# original image peak flux (Jy/beam): {image_data.max()} \n'
        f.write(catalog_out)
        catalog_out = f'# noise out (µJy/beam): {round(noise*1000000,2)} \n'
        f.write(catalog_out)
        limiting_flux = noise * threshold
        catalog_out = f'# cutt-off flux  (mJy/beam): {round(limiting_flux*1000,2)} \n'
        f.write(catalog_out)
        LOGGER.info('Submitting distributed tasks. This might take a while...') 
        source_list = multiprocess_contours(contours, image_data, fitsinfo, noise, args.ncpu)
        catalog_out = f"# freq0 (Hz): {fitsinfo['freq0']} \n"
        f.write(catalog_out)
        catalog_out = f'# number of sources detected: {len(source_list)} \n'
        f.write(catalog_out)
        catalog_out = '#\n#format: name ra_d dec_d i i_err emaj_s emin_s pa_d\n'
        f.write(catalog_out)
        for i in range(len(source_list)):
            output = 'src' + str(i) + ' ' + source_list[i][1] + '\n'
            f.write(output)
        f.close()
        LOGGER.info(f'Source catalog saved: {args.outcatalog}')

    if args.gui:
        try:
            from breizorro.gui import display
        except ModuleNotFoundError:
            LOGGER.error("Running breizorro gui requires optional dependencies, please re-install with: pip install breizorro[gui]")
            raise('Missing GUI dependencies')

        LOGGER.info("Loading Gui ...")
        display(args.imagename or args.maskname, mask_image, args.outcatalog, source_list)

    LOGGER.info(f"Enforcing that mask to binary")
    mask_image = mask_image!=0
    mask_header['BUNIT'] = 'mask'

    shutil.copyfile(input_file, out_mask_fits)  # to provide a template
    flush_fits(mask_image, out_mask_fits, mask_header)
    LOGGER.info("Done")
