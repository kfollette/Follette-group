from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

def compare(pyscale,idlscale):
    hdulist = fits.open(pyscale)
    ps = hdulist[0].data
    hdulist.close()
    hdulist = fits.open(idlscale)
    idls = hdulist[0].data
    hdulist.close()
    
    hdulist = fits.open('cosmic_arr.fits')
    cosmics = hdulist[0].data
    hdulist.close()
    count=0
    for a in range (idls.shape[0]):
        if not a in cosmics:
            idls = np.delete(idls,a-count)
            count = count+1

    plt.figure(1, figsize=(10,10))
    plt.plot(ps,idls,'co', markersize=3)
    plt.xlabel('python scale factor')
    plt.ylabel('IDL scale factor')
    plt.axis([0.98,1.14,0.98,1.14])
    plt.savefig('scale_comparison.png', overwrite=True)
    plt.close()