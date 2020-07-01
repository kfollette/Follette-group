import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
import glob
from astropy import units as u
from astropy import constants as const
import time
import pandas as pandas 
import sys

"""
Script to generate formatted txt file for estimating stellar parameters based on spectral type. 
Input should be the name of the object and its corresponding spectral type.

Last updated: 
    2020-06-30, v1: (KWD) Began script from original excel file (TempSpTyRelation.xlsx)
"""

def numericalspty(targname, spty, spty_err):
    """
    numericalspty
    ---------------
    Converts a list of spectral types into numerical values for interpolation purposes.
    The conversion is anchored at M0 = 50, with each +/- subclass corresponding to +/- 1, including 
        decimals for fractional spectral types.
    Code cannot currently handle additional spty information beyond class+subclass (e.g., 'M5Ve'
        for an emission line star) and ignores any characters beyond two decimal places.
    
    
    inputs
    ---------------
    targname                 : (array, str) array of target names formatted as strings
    spty                     : (array, str) array of spty values with class and subclass
    spty_err                 : (array, float) default = False; array of spty errors
    
    
    returns
    ---------------
    spty_num                 : (array) numerical spectral type
    spty_num_err             : (array) y-shift in pixels
    {filename}_SpTy.txt      : (text file) output of the covariance estimation                      

    
    
    dependencies
    ---------------
    none           
    
    
    todo
    ---------------
    - more flexible approach to 'unusual' spectral types?
    """
    if len(targname) != len(spty):
        raise V
    for target in targname:


    return spty_num, spty_num_err


    