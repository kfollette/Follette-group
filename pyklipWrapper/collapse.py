import numpy as np
import sys
from astropy.io import fits

filename = sys.argv[1]
kl = list(map(int, sys.argv[2].split(",")))
file = fits.getdata(filename)[:,0,:,:]
all_kl = list(map(int, fits.open(filename)[0].header['KLMODES'][1:-1].split(",")))

collapsed = np.zeros((1, file.shape[1], file.shape[2]))
for i in kl:
    collapsed += file[all_kl.index(i)]

collapsed = collapsed/len(kl)
hdu = fits.PrimaryHDU(collapsed)
hdulist = fits.HDUList([hdu])
hdulist.writeto(filename[:-5] + "_collapsed.fits", overwrite=True)
