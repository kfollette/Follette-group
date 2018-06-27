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

def process_flats(path_to_raw_flats, path_to_raw_darks, imsize, flatcategory, objname):
    '''
    process_flats
    ----------
    Routine to combine individual flat frames into a master flat. Depending on type of flat, 
    will: (0 - SKY) find corresponding darks to the flat exposure times, median stack the darks, and 
    subtract the master dark from each flatfield, then stack the flats into a master flat.
    Or, will (1 - LAMP) determine which of the flats are "on" and "off" and combine each of the on 
    and off sequences into a master "on" and master "off", then subtract them.
        
    Inputs
    ----------
    path_to_raw_flats       : (string) path to location of flatfield directory for given target; usu. './FLAT/'
    path_to_raw_darks	    : (string) path to location of dark directory for given target; usu. './DARK/'
    imsize                  : (int) dimension of FITS image array
    flatcategory            : (boolean) describe whether flats are lamps (1) or sky/twilight (0)
    objname                 : (string) name of target, e.g., "HIP123"
        
    Returns
    ----------
    med_flat                : (array) median-combined flat frame
    master_flat             : (array) normalized, combined flatframe
    flatname                : (FITS file) output fits file of median flatfield (can be lamp or sky)
    masterflatname          : (FITS file) output fits file of master flatfield

    Dependents
    ----------
    flat_exptime            : (float)
    dark_combine            : (function)
    '''

    flatlist = glob.glob(path_to_raw_flats + '*fits')
    n = len(flatlist)

    flatheader = fits.getheader(flatlist[0], ignore_missing_end=True)

    date = flatheader['DATE']
    date = date[0:10]

    # if flats are lamp flats:
    if flatcategory == 1:

        # check exposure times: lamp flats should all be the same
        flattimes = [fits.getheader(im)['EXPTIME'] for im in flatlist]
        if all(x == flattimes[0] for x in flattimes):
            print("Exposure time: "+str(flattimes[0]))
        else:
            raise Exception("Exposure times for given list of lamp flats do not match.")


        flatarray = np.zeros([imsize, imsize, n])

        for ii in range(0,n):
            flatdata = fits.getdata(flatlist[ii], ignore_missing_end=True)
            if len(flatdata.shape) == 3: #check for data cubes
                assert not np.any(np.isnan(flatdata))
                flatdata = np.median(flatdata, axis=0) #if data cube, then median frames
            flatarray[:,:,ii] = flatdata

        tmp_medval = np.average(flatarray)
        flaton = np.zeros([imsize,imsize,n])
        flatoff = np.zeros([imsize,imsize,n])

        counton = 0
        countoff = 0    

        for ii in range(0, n):
            if np.median(flatarray[:,:,ii]) < tmp_medval:
                flatoff[:,:,countoff] = flatarray[:,:,ii]
                countoff += 1
            else:
                flaton[:,:,counton] = flatarray[:,:,ii]
                counton += 1
        
        medflatoff = np.median(flatoff, axis = 2)
        medflaton = np.median(flaton,axis = 2)    
        med_flat = medflaton - medflatoff
        flatname = objname + '_' + date + '_medianlampflat.fits'
        fits.writeto(flatname, med_flat, flatheader, overwrite = True, output_verify='silentfix')            



    # if flats are twilight/skyflats:
    if flatcategory == 0:
        
        # determine if a master dark with appropriate exposure time exists in main directory:
        exptime = flatheader['EXPTIME']
        masterdark = glob.glob('*masterdark*'+str(exptime)+'*fits')
        
        if len(masterdark) == 0:
            print('Creating new master dark for flats with '+str(exptime)+'s exposure.')
            med_dark = dark_combine(path_to_raw_darks, exptime, imsize, objname)

        elif len(masterdark) > 1:
            print('Error: Too many matching darks for given exptime.')

        else:
            med_dark = fits.getdata(masterdark[0])

        # subtract off the median dark frame from each of the twiflats, IF median dark is the same exposure length as the twilight flats. 
        for i in range (0,n):
            flatarray[:,:,i] -= med_dark
        
        median_flat1 = np.median(flatarray[:,:,0])
        
        for i in range(0,n):
            flatmedian = np.median(flatarray[:,:,i])
            flatarray[:,:,i] *= (median_flat1/flatmedian)

        med_flat = np.median(flatarray,axis=2)
        flatname = objname + '_' + date + '_medianskyflat.fits'
        fits.writeto(flatname, med_flat, flatheader, overwrite = True, output_verify='silentfix')

    # some housekeeping to free up memory (at some point)
    del flatarray
    del flaton
    del flatoff
        
    ind = np.where(med_flat == 0)
    if np.sum(ind) > 0:    
        med_flat[ind] = 0.001


    # Create the master flat: scale median flat by its median value, to get pixel values near 1.0
    flatscalefactor = np.nanmedian(med_flat) # instead of np.median (Jan 2016)

    master_flat = med_flat/flatscalefactor

    if flatcategory == 1:
        masterflatname = objname + '_' + date + '_masterlampflat.fits'
    elif flatcategory == 0:
        masterflatname = objname + '_' + date + '_masterskyflat.fits'

    fits.writeto(masterflatname, master_flat, flatheader, overwrite = True, output_verify='silentfix')
    

    return med_flat, master_flat, flatheader
