import numpy as np
import scipy as sp
import scipy.ndimage
from astropy.io import fits
import sys
import numpy.ma as ma
import math
import pidly
import os
import glob

def badpixelmap(med_flat, objname, flatheader, max_std = 3.0):
    '''
    badpixelmap
    ----------
    Routine to create a map of the bad pixels from a median-combined flatfield. Determines which pixels have values 
    above or below some number of standard deviations of the flatfield median, and sets them to NaNs.
        
    Inputs
    ----------
    med_flat 	 		: (array) median-combined, dark-subtracted flatfield (output from flat_combine)
    objname				: (string) name of target, e.g., "HIP123"
    flatheader			: (FITS header) header object taken from an individual flatfield header
    max_std				: (float) number of standard deviations above which a pixel is considered bad. Default = 3
    
    Returns
    ----------
    badflat       		: (array) output image of the bad pixel map
    badpixname 			: (FITS file) output FITS file of the bad pixel map

    Dependents
    ----------
    None
    '''
    date = flatheader['DATE']
    date = date[0:10]

    # calculate the standard deviation of the median stacked flat, not including the edges of the array
    std = np.std(med_flat[100:924,100:1000])

    # some trickery so now that within each pixel is the number of deviations that pixel is away from the median value
    tmp = abs(med_flat - np.median(med_flat))/std  

    # create a list of array indices where the standard deviation is above the max_std.
    bad_locations = tmp > max_std

    badflat = np.copy(med_flat)
    badflat[bad_locations]= np.nan # set to NaNs instead of zero (Jan 2016)
    badpixname = objname + '_' + date +'_badflat.fits'
    fits.writeto(badpixname, badflat, flatheader, overwrite = True, output_verify='silentfix')

    return badflat
