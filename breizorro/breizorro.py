import sys
import numpy
import shutil
import logging
import argparse
import scipy.special
import scipy.ndimage
from astropy.io import fits
from argparse import ArgumentParser
from scipy.ndimage.measurements import label
from scipy.ndimage.morphology import binary_dilation

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
    LOGGER.info(f"Reading {fitsfile} data")
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
        image = numpy.array(input_hdu.data[:, :])
    elif len(input_hdu.data.shape) == 3:
        image = numpy.array(input_hdu.data[0, :, :])
    else:
        image = numpy.array(input_hdu.data[0, 0, :, :])
    return image


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

def merge_masks(maskfiles, dilate=0):
    """Merge a list of mask file and perfom dilation"""
    shutil.copyfile(modelfits[0], outfile)
    init_mask = maskfile[0]
    finalmaskdata = getImage(init_mask)
    for maskfile in maskfiles[1:]:
        finalmaskdata += getImage(maskfile)
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
    parser.add_argument('--restored-image', dest="imagename", required=False,
                        help="Restored image file")
    parser.add_argument('--mask-image', dest="masknames", nargs='+', type=argparse.FileType('r'),
                        help="Input mask file")
    parser.add_argument('--merge', dest='merge', action='store_true', default=False,
                        help='Merge a list of masks')
    parser.add_argument('--threshold', dest='threshold', default=6.5,
                        help='Sigma threshold for masking (default = 6.5)')
    parser.add_argument('--boxsize', dest='boxsize', default=50,
                        help='Box size over which to compute stats (default = 50)')
    parser.add_argument('--dilate', dest='dilate', default=0,
                        help='Number of iterations of binary dilation (default=0)')
    parser.add_argument('--number-islands', dest='islands', action='store_true', default=False,
                        help='Number the islands detected (default=do not number islands)')
    parser.add_argument('--remove-islands', dest='remove_isl', metavar='N', type=int, nargs='+',
                         help='List of islands to remove from input mask. e.g. --remove-islands 1,18,20')
    parser.add_argument('--savenoise', dest='savenoise', action='store_true', default=False,
                        help='Enable to export noise image as FITS file (default=do not save noise image)')
    parser.add_argument('--outfile', dest='outfile', default='',
                        help='Suffix for mask image (default=restored_image.replace(".fits",".mask.fits"))')
    args = parser.parse_args()
    threshold = float(args.threshold)
    boxsize = int(args.boxsize)
    dilate = int(args.dilate)
    savenoise = args.savenoise
    outfile = args.outfile

    if not args.imagename and not args.maskname:
        LOGGER.info("Please specify a FITS file")
        sys.exit()


    if args.imagename:
        input_fits = args.imagename.rstrip('/')
        input_image = get_image(input_fits)

        noise_image = make_noise_map(input_image, boxsize)
        if savenoise:
            noise_fits = input_fits.replace('.fits', '.noise.fits')
            shutil.copyfile(input_fits, noise_fits)
            flush_fits(noise_image, noise_fits)

        mask_image = input_image > threshold * noise_image

        mask_image[:, -1]=0
        mask_image[:, 0]=0
        mask_image[0, :]=0
        mask_image[-1, :]=0

        if dilate != 0:
            LOGGER.info(f"Dilation mask, {dilate} iteration(s)")
            dilated = binary_dilation(input=mask_image, iterations=dilate)
            mask_image = dilated

        if outfile == '':
            out_mask_fits = input_fits.replace('.fits', '.mask.fits')
        else:
            out_mask_fits = outfile

        shutil.copyfile(input_fits, out_mask_fits)
        flush_fits(mask_image, out_mask_fits)

    if args.masknames:
        import IPython; IPython.embed()
        input_fits = args.masknames[0].rstrip('/')
        input_mask_image = get_image(args.maskname[0])
        mask_header = get_image_header(args.maskname[0])
        if args.islands:
            input_mask_image = input_mask_image.byteswap().newbyteorder()
            labeled_mask, num_features = label(input_mask_image)
            input_mask_image = labeled_mask
            mask_header['BUNIT'] = 'source_ID'
            LOGGER.info(f"Number of islands: {num_features}")
            out_mask_fits = input_fits.replace('.fits', '_isl.fits')
        elif args.remove_isl:
            LOGGER.info(f"Removing islands: {args.remove_isl}")
            for isl in args.remove_isl:
                input_mask_image = numpy.where(input_mask_image==isl, 0, input_mask_image)
        elif args.merge:
            LOGGER.info(f"Merging the following masks: {args.masknames}")
            input_mask_image = merge_mask(args.masknames, args.dilate)
        else:
            LOGGER.info(f"All islands are converted to 1")
            input_mask_image[input_mask_image>=1] = 1
            mask_header['BUNIT'] = 'Jy/beam'

        mask_image=input_mask_image
        if outfile == '':
            out_mask_fits = input_fits.replace('.fits', '_isl.fits')
        else:
            out_mask_fits = outfile

        shutil.copyfile(input_fits, out_mask_fits)
        flush_fits(mask_image, out_mask_fits, mask_header)


    LOGGER.info("Done")
