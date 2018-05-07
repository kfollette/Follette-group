"""
Takes in final line reduction, final continuum reduction, and list of scale factors
Creates final SDI reduction by subtracting cont*avg_scale from line
"""

from astropy.io import fits
import numpy as np

def run(line,cont,scales,KL):
    scale = np.mean(fits.getdata(scales))
    line = fits.getdata(line)[KL] 
    cont = fits.getdata(cont)[KL]
    SDI = line-(cont*scale)
    print("scale factor is ", scale)
    fits.writeto("Final_SDI_result.fits", SDI, overwrite=True)
