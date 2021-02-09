import sys
import numpy
import shutil
import logging
import argparse
import os.path
import numpy as np
from astropy.io import fits
from argparse import ArgumentParser

from scipy.ndimage.morphology import binary_dilation, binary_fill_holes
from scipy.ndimage.measurements import label
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
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(cfmt)
    log.addHandler(console)
    return log

LOGGER = create_logger()

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


def merge_mask(finalmaskdata, maskfiles, dilate=0):
    """Merge a list of mask file and perfom dilation"""
    for maskfile in maskfiles:
        finalmaskdata |= get_image(maskfile)[0]
    if dilate:
        finalmaskdata = binary_dilation(finalmaskdata, iterations=dilate)
    return finalmaskdata


def subtract_mask(maskfiles, dilate=0):
    """Subtract a list of masks from the index 0 maskfile and perfom dilation"""
    init_mask = maskfiles[0]
    finalmaskdata = get_image(init_mask)
    for maskfile in maskfiles[1:]:
        finalmaskdata -= get_image(maskfile)
    finalmaskdata[finalmaskdata<0] = 0
    if dilate:
        finalmaskdata = binary_dilation(finalmaskdata, iterations=dilate)
    return finalmaskdata


def main():
    LOGGER.info("Welcome to breizorro")
    # Get version
    from pkg_resources import get_distribution
    _version = get_distribution('breizorro').version
    LOGGER.info(f"Version: {_version}")

    LOGGER.info("Usage: breizorro --help")
    parser = ArgumentParser(description='breizorro [options] --image restored_image')
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

    parser.add_argument('--merge', dest='merge', metavar="MASK(s)", nargs='+',
                        help='Merge in one or more masks')
    parser.add_argument('--subtract', dest='subtract', metavar="MASK(s)", nargs='+',
                        help='Subract one or more masks')

    parser.add_argument('--number-islands', dest='islands', action='store_true', default=False,
                        help='Number the islands detected (default=do not number islands)')
    parser.add_argument('--remove-islands', dest='remove_isl', metavar='N', type=int, nargs='+',
                         help='List of islands to remove from input mask. e.g. --remove-islands 1 18 20')
    parser.add_argument('--extract-islands', dest='extract_isl', metavar='N', type=int, nargs='+',
                         help='List of islands to extract from input mask. e.g. --extract-islands 1 18 20')
    parser.add_argument('--make-binary', action="store_true",
                         help='Replace all island numbers with 1')
    
    parser.add_argument('--dilate', dest='dilate', metavar="R", type=int, default=0,
                        help='Apply dilation with a radius of R pixels')
    parser.add_argument('--fill-holes', dest='fill_holes', action='store_true', 
                        help='Fill holes (i.e. entirely closed regions) in mask')

    parser.add_argument('-o', '--outfile', dest='outfile', default='',
                        help='Suffix for mask image (default based on input name')

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

    # define input file, and get its name and extension
    input_file = args.imagename or args.maskname
    name, ext = os.path.split(input_file)

    # first, load or generate mask

    if args.imagename:
        input_image, input_header = get_image(args.input_file)
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

    # next, merge and/or subtract
    if args.merge:
        for merge in args.merge:
            mask_image += get_image(merge)[0]
            LOGGER.info("Merged into mask")
        mask_image = mask_image != 0
        mask_header['BUNIT'] = 'mask'

    if args.subtract:
        for subtract in args.subtract:
            mask_image[get_image(subtract)[0] != 0] = 0
            LOGGER.info("Subtracted from mask")

    if args.islands:
        LOGGER.info(f"(Re)numbering islands")
        mask_image = mask_image != 0
        # mask_image = mask_image.byteswap().newbyteorder()
        mask_image, num_features = label(mask_image)
        mask_header['BUNIT'] = 'Source_ID'
        LOGGER.info(f"Number of islands: {num_features}")
    
    if args.remove_isl:
        LOGGER.info(f"Removing islands: {args.remove_isl}")
        for isl in args.remove_isl:
            mask_image[mask_image == isl] = 0

    if args.extract_isl:
        LOGGER.info(f"Extracting islands: {args.extract_isl}")
        new_mask_image = np.zeros_like(mask_image)
        for isl in args.extract_isl:
            new_mask_image[mask_image == isl] = isl
        mask_image = new_mask_image

    if args.make_binary:
        LOGGER.info(f"Converting mask to binary")
        mask_image = mask_image!=0
        mask_header['BUNIT'] = 'mask'

    if args.dilate:
        LOGGER.info(f"Dilating mask using a ball of R={args.dilate}pix")
        R = args.dilate
        r = np.arange(-R, R+1)
        struct = np.sqrt(r[:, np.newaxis]**2 + r[np.newaxis,:]**2) <= R
        # print(struct)
        mask_image = binary_dilation(input=mask_image, structure=struct)
        
    if args.fill_holes:
        LOGGER.info(f"Filling closed regions")
        binary_fill_holes(mask_image, output=mask_image)

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
