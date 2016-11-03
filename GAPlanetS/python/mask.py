import numpy as np
from astropy.io import fits
import glob

imlist = glob.glob('sliced*.fits')
mask = fits.getdata('ctrlrad_mask_27to42.fits')
nims = len(imlist)

for index in range(nims):
    im, head = fits.getdata(imlist[index], header=True)
    im = im * mask
    fits.writeto(imlist[index], im, clobber=True, header=head)
