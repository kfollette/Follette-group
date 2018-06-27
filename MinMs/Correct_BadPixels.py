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


def correct_bad_pixels(sciarray, badflat, width = 5, sd_cut = 4.0):
    '''
    correct_bad_pixels
    ----------
    For a given array of FITS files, interpolates over bad pixels, and identifies and cleans cosmic rays. 
        
    Inputs
    ----------
    sciarray 	 		: (array) stacked (x by y by n) array of science image data (stack of science images)
    badflat 			: (array) image of the bad pixel map
    width				: (int) size in pixels of the median filter to apply. Default = 5 pixels
    sd_cut				: (float) standard deviation of pixel threshold for hot pixel/cosmic ray identification. Default = 4.0
    
    Returns
    ----------
    sciarray_corr       : (array) stacked array of science frames with bad pixels corrected 
    
    Dependents
    ----------
    None
    '''
    # initialize a new empty array with the same dimensions as sciarray
    sciarray_corr = np.zeros(sciarray.shape)
    
    # make a copy of sciarray which we can edit as needed
    sciarray_copy = np.copy(sciarray)
    
    for ii in range(0,n):

        # cleans the bad pixels identified from the flatfield
        med_im = sp.ndimage.median_filter(sciarray_copy[:,:,ii], width)
        tmp = sciarray_copy[:,:,ii]
        tmp[np.isnan(badflat)] = med_im[np.isnan(badflat)]
        sciarray_corr[:,:,ii] = tmp

		# This part is working out the local variance, and cleans cosmic rays
        med_im = sp.ndimage.median_filter(sciarray_corr[:,:,ii],width)
        av_im = sp.ndimage.filters.uniform_filter(sciarray_corr[:,:,ii],width)
        avsq_im = sp.ndimage.filters.uniform_filter(pow(sciarray_corr[:,:,ii],2.0),width)
        var_im = avsq_im - (pow(av_im,2.0))
        sd_im = np.sqrt(var_im)

        ind = np.where(abs((sciarray_corr[:,:,ii] - av_im)/sd_im) > sd_cut)
        if np.sum(ind) > 0:
            tmp = sciarray_corr[:,:,ii]
            tmp[ind] = med_im[ind]
            sciarray_corr[:,:,ii] = tmp

    del med_im
    del av_im
    del avsq_im
    del var_im
    del sd_im

    return sciarray_corr
