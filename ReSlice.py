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

###INPUTS###
hdulist = fits.open(sys.argv[1])
Line_data = hdulist[0].data
hdulist.close()
hdulist = fits.open(sys.argv[2])
Cont_data = hdulist[0].data
hdulist.close()

###MAIN CODE###
size = Line_data.shape[0]
size2 = Cont_data.shape[0]
if(not size == size2):
    sys.exit("Line and Continuum data cubes are not the same size.")

count = 1
for z in range(size):
    SL = Line_data[z]
    SC = Cont_data[z]
    NewSlice = np.vstack(([SL],[SC]))
    hdu = fits.PrimaryHDU(NewSlice)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("Cube_" + str(count) + ".fits", overwrite=True)
    count += 1
