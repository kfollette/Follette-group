"""
Takes in 2 cubes, both with the third dimension being a time series, 
one of which is in H-alpha and one of which is in the continuum wavelength.
Outputs multiple 2-layer cubes, one for each H-alpha/continuum pair

ex: 2 1000-layer cubes input, 1000 2-layer cubes output.

Original cube files are command-line inputs, the first being H-alpha and the second being continuum.

Last modified: 6/30/2017
"""

import sys
from astropy.io import fits
import numpy as np
import os

def slice(linefile, contfile, rotoff):
    
    ###INPUTS###
    hdulist = fits.open(linefile)
    linecube = hdulist[0].data
    header_data = fits.getheader(linefile, 0) #will be the same as the cont header
    hdulist.close()
    hdulist = fits.open(contfile)
    contcube = hdulist[0].data
    hdulist.close()
    hdulist = fits.open(rotoff)
    r = hdulist[0].data
    hdulist.close()
    
    ###OUTPUT PATH###
    if not os.path.exists("sliced"):
        os.makedirs("sliced")
    os.chdir("sliced")
    
    ###MAIN CODE###
    size = linecube.shape[0]
    size2 = contcube.shape[0]
    size3 = r.shape[0]
    if(not size == size2 and size == size3):
        sys.exit("Line, Continuum, and rotoff data are not all the same size.")
        
    for z in range(size):
        SL = linecube[z]
        SC = contcube[z]
        NewSlice = np.vstack(([SL],[SC]))
        hdu = fits.PrimaryHDU(NewSlice, header_data)
        hdulist = fits.HDUList([hdu])
        prihdr = hdulist[0].header
        prihdr.set('rotoff', str(r[z]))
        hdulist.writeto("Cube_" + str(z+1) + ".fits", overwrite=True)
        print("Sliced layer " + str(z+1))
        
    print("Done slicing.")