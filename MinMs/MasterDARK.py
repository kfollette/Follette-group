
import numpy as np
import scipy as sp
import scipy.ndimage
from astropy.io import fits
import sys
import numpy.ma as ma
import math
import os
import glob
import fnmatch
import shutil

def dark_combine(path_to_raw_darks, sci_exptime, imsize, objname):
	'''
	dark_combine
    ----------
    Routine to find corresponding darks to the science exposure times and median stack them. 
    Requires darks to exist in folder in upper level directory.
        
    Inputs
    ----------
    path_to_raw_darks	: (string) path to location of dark directory for given target; usu. './DARK/'
    sci_exptime    		: (float) exposure time of the science files to be calibrated
    imsize		   		: (int) dimension of FITS image array
    objname             : (string) name of target directory 
    
    Returns
    ----------
    med_dark       		: (array) median-combined dark frame
    
    darkname       		: (FITS file) output fits file of master dark frame and header taken from first raw dark.

    Dependents
    ----------
    sci_exptime    		: (float)
    
    History
    ----------
    (v1) 2018-06-19 - KWD: Updated from original script for Python 3 compatibility and astropy updates.
    (v2) 2018-06-20 - SGC: Added objname parameter.
    '''
	darklist = glob.glob(path_to_raw_darks + '*fits')
	n = len(darklist)

	im = fits.open(darklist[0], ignore_missing_end=True)
	darkheader = im[0].header

	date = darkheader['DATE']
	date = date[0:10]


	# Make a new list of darks that match the provided exptime:
	matching_darks = []

	# Select only those darks that match the given exposure time:
	for each_dark in darklist:
		header = fits.getheader(each_dark)
		if header['EXPTIME'] == sci_exptime:
			matching_darks.append(each_dark)
		else:
			pass

	number_of_darks = len(matching_darks)

	print(f'Found {number_of_darks} darks with exposure times of {sci_exptime}.')

	if number_of_darks == 0:
		raise ValueError('There are no matching dark frames.')

	darkarray = np.zeros([imsize, imsize, number_of_darks])

	for ii in range(0, number_of_darks):
		header = fits.getheader(matching_darks[ii])

		# check to add only those those darks with correct exptime
		if header['EXPTIME'] == sci_exptime:
			im = fits.getdata(matching_darks[ii]) 
			if len(im.shape) == 3: # check for data cubes
				assert not np.any(np.isnan(im))
				im = np.median(im, axis=0) # if data cube, then median frames
			darkarray[:,:,ii] = im

	med_dark = np.median(darkarray, axis=2)
	darkname = objname + '_' + date + '_masterdark_' + str(sci_exptime) + 's.fits'
	fits.writeto(darkname, med_dark, darkheader, overwrite=True, output_verify='silentfix')
	return med_dark