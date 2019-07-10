import numpy as np
import sys
import glob
from astropy.io import fits
import os
from File_Display_Sorter import *

cwd = os.getcwd()

allfits = glob.glob('*.fits')

for eachfile in allfits:
    header = fits.getheader(eachfile, ignore_missing_end=True, silentfix=True)
    objname = header['OBJECT']
    if os.path.isdir("./%s" % objname) == False:
        os.mkdir("./%s" % objname)
        print(f"\nNow creating a new directory for {objname}") 
        os.rename(eachfile,"./%s/%s" % (objname, eachfile))
    elif os.path.isdir("./%s" % objname) == True:
        os.rename(eachfile, "./%s/%s" % (objname, eachfile))
    else:
        print("Something has gone wrong!")

# Sort into subfolders based on filter
alldir = glob.glob('HIP*')

for dir in alldir:
    os.chdir('./' + dir)
    fitslist = glob.glob('*.fits')
    
    filters = []

    for idx, im in enumerate(fitslist):
        hdr = fits.getheader(im)
        filters.append(hdr['FILTER'])

    for filt in np.unique(filters):
        sorter('*fits', filt, 'FILTER')
        
        filt_new = filt.split(' ')[0] + filt.split(' ')[1]        

        os.chdir('./' + filt_new)
        
        newfits = glob.glob('*.fits')
        
        exptimes = []

        for idx, im in enumerate(newfits):
            hdr = fits.getheader(im)
            exptimes.append(hdr['EXPTIME'])

        for exp in np.unique(exptimes):
            exptime('*fits', exp, str(exp))

        os.chdir('../')
    
    os.chdir('../')   
