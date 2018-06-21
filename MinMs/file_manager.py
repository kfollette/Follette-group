
# coding: utf-8

# In[2]:

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

"""
fits_manager.py: Script to organize fits files of a given HIP directory according to their calibration type
Usage: python fits_manager.py.py

History:
(v1) 2018-06-13 - SGC
"""


# In[10]:

HIP_list = glob.glob('*.fits')

def file_manager(HIP_list):
    '''
    fileinfo
    ----------
    Organizes fits files acording to their calibration (DARK, OBJECT, SKY, FLAT_LAMP) for a given HIP target.
        
    Inputs
    ----------
    HIP_list      : (list) list of fits filenames, each in string format
    
    Returns
    ----------
    (print output) : prints out the HIP directory that it is currently in. Then prints out the name of the type of calibration for each fits files to check if code is running as it should.

    Dependents
    ----------
    None
    '''
    print(HIP_list)
    for files in HIP_list:
    
        header = fits.getheader(files)
    
        if header['HIERARCH ESO DPR TYPE'] == 'DARK':
            print("dark")
            os.rename(files, 'DARK/' + files)
        
        elif header['HIERARCH ESO DPR TYPE'] == 'OBJECT':
            print("object")
            os.rename(files, 'OBJECT/' + files)
        
        elif header['HIERARCH ESO DPR TYPE'] == 'SKY':
            print("sky")
            os.replace(files, 'SKY/' + files)
        
        elif header['HIERARCH ESO DPR TYPE'] == 'FLAT,LAMP':
            print("flat,lamp")
            os.rename(files, 'FLAT_LAMP/' + files)
        else:
            print('error')
            

#file_manager(HIP_list)

