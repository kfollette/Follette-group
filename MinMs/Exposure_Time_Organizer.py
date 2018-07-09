import numpy as np
import scipy as sp
from astropy.io import fits
import sys
import os
import glob
import shutil

DARKlist= glob.glob('./DARK/*fits')
OBJlist = glob.glob('./OBJECT/*fits')
def exposure_time_check(DARKlist, OBJlist):
    
    print('DARK fits exp times:')
    for i in DARKlist:
        darkheader = fits.getheader(i)
        print(darkheader['EXPTIME'])
        
        if darkheader['EXPTIME'] == 0.5:
            shutil.move(i, './DARK/SHORT/')
        elif darkheader['EXPTIME'] > 0.5:
            shutil.move(i, './DARK/LONG/')
            
     
    print('OBJ fits exp times:')
    for j in OBJlist:
        objheader = fits.getheader(j)
        print(objheader['EXPTIME'])
        
        if objheader['EXPTIME'] == 0.5:
            shutil.move(j, './OBJECT/SHORT/')
        elif objheader['EXPTIME'] > 0.5:
            shutil.move(j, './OBJECT/LONG/')
        
########################################################
  


    
    
    