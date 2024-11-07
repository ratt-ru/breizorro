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

from regions import PixCoord
from regions import PolygonSkyRegion, PolygonPixelRegion

def deg_to_hms(ra_deg):
    ra_hours = ra_deg / 15  # 360 degrees = 24 hours
    hours = int(ra_hours)
    minutes = int((ra_hours - hours) * 60)
    seconds = (ra_hours - hours - minutes / 60) * 3600
    return hours, minutes, seconds


def deg_to_dms(dec_deg):
    degrees = int(dec_deg)
    dec_minutes = abs(dec_deg - degrees) * 60
    minutes = int(dec_minutes)  
    seconds = (dec_minutes - minutes) * 60
    return degrees, minutes, seconds


def format_source_coordinates(ra_deg, dec_deg):
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    ra_str = coord.ra.to_string(unit=u.hour, sep=':', precision=2)
    dec_str = coord.dec.to_string(unit=u.deg, sep=':', precision=2)
    return ra_str, dec_str

def format_source_coordinates0(coord_ra_deg, coord_dec_deg):
    h,m,s = deg_to_hms(coord_ra_deg)
    if h < 10:
        h  = '0' + str(h)
    else:
        h = str(h)
    if m < 10: 
        m = '0' + str(m)
    else:
        m = str(m)
    s = round(s,2)
    if s < 10:
        s = '0' + str(s)
    else:
        s = str(s)
    if len(s) < 5:
        s = s + '0'
    h_m_s = h + ':' + m + ':' + s

    d,m,s = deg_to_dms(coord_dec_deg)
    if d >= 0 and d < 10:
        d = '0' + str(d)
    elif d < 0 and abs(d) < 10: 
        d = '-0' + str(d)
    else:
        d = str(d)
    if m < 10:
        m = '0' + str(m)
    else:
        m = str(m)
    s = round(s,2)
    if s < 10:
        s = '0' + str(s)
    else:
        s = str(s)
    if len(s) < 5:
        s = s + '0'
    d_m_s = d + ':' + m + ':' + s

    src_pos = (h_m_s,  d_m_s)
    return src_pos

def deg2dec(dec_deg, deci=2):
    """Converts declination in degrees to dms coordinates

    Parameters
    ----------
    dec_deg : float
      Declination in degrees
    dec: int
      Decimal places in float format

    Returns
    -------
    dms : str
      Declination in degrees:arcmin:arcsec format
    """
    DD          = int(dec_deg)
    dec_deg_abs = np.abs(dec_deg)
    DD_abs      = np.abs(DD)
    MM          = int((dec_deg_abs - DD_abs)*60)
    SS          = round((((dec_deg_abs - DD_abs)*60)-MM), deci)
    return "%s:%s:%s"%(DD,MM,SS)


def deg2ra(ra_deg, deci=2):
    """Converts right ascension in hms coordinates to degrees

    Parameters
    ----------
    ra_deg : float
    ra in degrees format

    Returns
    -------
    HH:MM:SS : str

    """
    if ra_deg < 0:
       ra_deg = 360 + ra_deg
    HH     = int((ra_deg*24)/360.)
    MM     = int((((ra_deg*24)/360.)-HH)*60)
    SS     = round(((((((ra_deg*24)/360.)-HH)*60)-MM)*60), deci)
    return "%s:%s:%s"%(HH,MM,SS)


def calculate_area(bmaj, bmin, pix_size):
    """
    Calculate the area of an ellipse represented by its major and minor axes,
    given the pixel size.

    Parameters:
        bmaj (float): Major axis raduis of the ellipse in arcseconds.
        bmin (float): Minor axis radius of the ellipse in arcseconds.
        pix_size (float): Pixel size in arcseconds.

    Returns:
        float: Calculated area of the ellipse in square pixels.
    """
    # Calculate the semi-major and semi-minor axes in pixels
    a_pixels = bmaj / pix_size
    b_pixels = bmin / pix_size

    # Calculate the area of the ellipse using the formula: π * a * b
    area = np.pi * a_pixels * b_pixels

    return area


def fitsInfo(fitsname=None):
    """Get fits header info.

    Parameters
    ----------
    fitsname : fits file
        Restored image (cube)

    Returns
    -------
    fitsinfo : dict
        Dictionary of fits information
        e.g. {'wcs': wcs, 'ra': ra, 'dec': dec,
        'dra': dra, 'ddec': ddec, 'raPix': raPix,
        'decPix': decPix,  'b_size': beam_size,
        'numPix': numPix, 'centre': centre,
        'skyArea': skyArea, 'naxis': naxis}

    """
    with fits.open(fitsname) as hdu:
        hdr = hdu[0].header
    ra = hdr['CRVAL1']
    dra = abs(hdr['CDELT1'])
    raPix = hdr['CRPIX1']
    dec = hdr['CRVAL2']
    ddec = abs(hdr['CDELT2'])
    decPix = hdr['CRPIX2']
    wcs = WCS(hdr)
    numPix = hdr['NAXIS1']
    naxis = hdr['NAXIS']
    try:
        beam_size = (hdr['BMAJ'], hdr['BMIN'], hdr['BPA'])
    except:
        beam_size = None
    try:
        centre = (hdr['CRVAL1'], hdr['CRVAL2'])
    except:
        centre = None
    try:
        freq0=None
        for i in range(1, hdr['NAXIS']+1):
            if hdr['CTYPE{0:d}'.format(i)].startswith('FREQ'):
                freq0 = hdr['CRVAL{0:d}'.format(i)]
    except:
        freq0=None

    skyArea = (numPix * ddec) ** 2
    fitsinfo = {'wcs': wcs, 'ra': ra, 'dec': dec, 'naxis': naxis,
                'dra': dra, 'ddec': ddec, 'raPix': raPix,
                'decPix': decPix, 'b_size': beam_size,
                'numPix': numPix, 'centre': centre,
                'skyArea': skyArea, 'freq0': freq0}
    return fitsinfo


def get_image_data(fitsfile):
    """
    Reads a FITS file and returns a tuple of the image array and the header.
    
    Parameters:
    fitsfile (str): The path to the FITS file.
    
    Returns:
    tuple: A tuple containing the image array and the header.
    """
    with fits.open(fitsfile) as input_hdu:
        # Check the dimensionality of the data and extract the appropriate slice
        if len(input_hdu[0].data.shape) == 2:
            image = np.array(input_hdu[0].data[:, :])
        elif len(input_hdu[0].data.shape) == 3:
            image = np.array(input_hdu[0].data[0, :, :])
        else:
            image = np.array(input_hdu[0].data[0, 0, :, :])
        
        header = input_hdu[0].header
    
    return image, header


def maxDist(contour, pixel_size, x_centroid, y_centroid):
    """Calculate maximum extent and position angle of a contour.

    Parameters:
    contour : list of [x, y] pairs
        List of coordinates defining the contour.
    pixel_size : float
        Size of a pixel in the image (e.g., arcseconds per pixel).
    x_centroid: int
        X-coordinate of the centroid.
    y_centroid: int
        Y-coordinate of the centroid.

    Returns:
    e_maj : float
        Major axis of the contour in angular units (e.g., arcseconds).
    e_min : float
        Minor axis of the contour in angular units (e.g., arcseconds).
    pos_angle : float
        Position angle of the contour (in degrees).
    """
    contour_array = np.array(contour)
    distances = np.linalg.norm(contour_array - [y_centroid, x_centroid], axis=1)

    # Find the indices of the points with maximum and minimum distances
    max_idx = np.argmax(distances)
    min_idx = np.argmin(distances)

    # Calculate major and minor axes
    e_maj = np.max(distances) * pixel_size
    e_min = np.min(distances) * pixel_size

    # Calculate position angle
    dx, dy = contour_array[max_idx] - [x_centroid, y_centroid]
    pos_angle = np.degrees(np.arctan2(dy, dx))

    return e_maj, e_min, pos_angle

def calculate_beam_area(bmaj, bmin, pix_size):
    """
    Calculate the area of an ellipse represented by its major and minor axes,
    given the pixel size.

    Parameters:
        bmaj (float): Major axis of the ellipse in arcseconds.
        bmin (float): Minor axis of the ellipse in arcseconds.
        pix_size (float): Pixel size in arcseconds.

    Returns:
        area (float): Calculated area of the ellipse in square pixels.
    """
    # Calculate the semi-major and semi-minor axes in pixels
    a_pixels = bmaj / pix_size
    b_pixels = bmin / pix_size

    # Calculate the area of the ellipse using the formula: π * a * b
    area = np.pi * a_pixels * b_pixels 

    return area


def get_source_size(contour, pixel_size, mean_beam, image, int_peak_ratio, centroids):
    result = maxDist(contour,pixel_size,*centroids)
    src_angle = result[:-1]
    pos_angle = result[-1]
    contour_pixels = PixCoord([c[0] for c in contour], [c[1] for c in contour])
    p = PolygonPixelRegion(vertices=contour_pixels, meta={'label': 'Region'})
    source_beam_ratio =  p.area / mean_beam
    # first test for point source
    point_source = False
    if (int_peak_ratio <= 0.2) or (src_angle[0] <= mean_beam):
        point_source = True
    if source_beam_ratio <=  1.0:
        point_source = True
    if point_source:
        src_size = (0.0, 0.0, 0.0)
    else:
        emaj = round(src_angle[0],2)
        emin = round(src_angle[-1],2)
        pa = round(pos_angle,2)
        src_size = (emaj, emin, pa)
    return src_size
