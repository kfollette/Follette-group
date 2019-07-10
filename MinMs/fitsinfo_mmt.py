import numpy as np
import scipy as sp
import scipy.ndimage
from astropy.io import fits
import sys
import numpy.ma as ma
import math
import os
import glob

fitslist = glob.glob('*fits')

n = len(fitslist)

for i in range(0,n):
        header=fits.getheader(fitslist[i])
        #fitstype = header['HIERARCH ESO DPR TYPE']
        fitsname = header['OBJECT']
        fitstime = header['EXPTIME']
        fitsrot = header['ROT']
        parang = header['PA']
        filt = header['FILTER']
        #print fitsname, fitstype, fitstime
        #print 'Object Name, Exptime, Rotation, Parallactic Angle, Filter'
        print(fitsname,',', fitstime,',', fitsrot,',', parang,',', filt)
