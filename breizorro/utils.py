import sys
import numpy
import shutil
import logging
import argparse
import os.path
import re
import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord


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


def format_source_coordinates(coord_ra_deg, coord_dec_deg):
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


def calculate_area(bmaj, bmin, pix_size):
    """
    Calculate the area of an ellipse represented by its major and minor axes,
    given the pixel size.

    Parameters:
        bmaj (float): Major axis of the ellipse in arcseconds.
        bmin (float): Minor axis of the ellipse in arcseconds.
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

def calculate_weighted_centroid(x, y, flux_values):
    # Calculate the total flux within the region  
    total_flux = np.sum(flux_values)
    # Initialize variables for weighted sums
    weighted_sum_x = 0
    weighted_sum_y = 0
    # Loop through all pixels within the region
    for xi, yi, flux in zip(x, y, flux_values):
        # Add the weighted contribution of each pixel to the centroid
        weighted_sum_x += xi * flux
        weighted_sum_y += yi * flux
    # Calculate the centroid coordinates
    centroid_x = weighted_sum_x / total_flux
    centroid_y = weighted_sum_y / total_flux
    return centroid_x, centroid_y


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
    hdu = fits.open(fitsname)   
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

def maxDist(contour, pixel_size):
    """Calculate maximum extent and position angle of a contour.

    Parameters:
    contour : list of [x, y] pairs
        List of coordinates defining the contour.
    pixel_size : float
        Size of a pixel in the image (e.g., arcseconds per pixel).

    Returns:
    ang_size : float
        Maximum extent of the contour in angular units (e.g., arcseconds).
    pos_angle : float
        Position angle of the contour (in degrees).
    """
    src_size = 0
    pos_angle = None

    # Convert the contour to a numpy array for easier calculations
    contour_array = np.array(contour)

    # Calculate pairwise distances between all points in the contour
    for i in range(len(contour_array)):
        for j in range(i+1, len(contour_array)):
            # Calculate Euclidean distance between points i and j
            distance = np.linalg.norm(contour_array[i] - contour_array[j]) * pixel_size

            # Calculate positional angle between points i and j
            dx, dy = contour_array[j] - contour_array[i]
            angle = np.degrees(np.arctan2(dy, dx))

            # Update max_distance, max_points, and pos_angle if the calculated distance is greater
            if distance > src_size:
                src_size = distance
                pos_angle = angle

    return src_size, pos_angle


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
    a_pixels = b_major / pixel_size
    b_pixels = b_minor / pixel_size

    # Calculate the area of the ellipse using the formula: π * a * b
    area = np.pi * a_pixels * b_pixels

    return area


def get_source_size(contour, pixel_size, mean_beam, image, int_peak_ratio):
    result = maxDist(contour,pixel_size)
    src_angle = result[0]
    pos_angle = result[1]
    contour_pixels = PixCoord([c[0] for c in contour], [c[1] for c in contour])
    p = PolygonPixelRegion(vertices=contour_pixels, meta={'label': 'Region'})
    source_beam_ratio =  p.area / mean_beam
    # first test for point source
    point_source = False
    if (int_peak_ratio <= 0.2) or (src_angle <= mean_beam):
        point_source = True
    if source_beam_ratio <=  1.0:
        point_source = True
    if point_source:
        src_size = (0.0, 0.0)
        print(f"Point source because {int_peak_ratio} <= 0.2 and {src_angle} <= {mean_beam}")
    else:
        ang = round(src_angle,2)
        pa = round(pos_angle,2)
        src_size = (ang, pa)
    return src_size
