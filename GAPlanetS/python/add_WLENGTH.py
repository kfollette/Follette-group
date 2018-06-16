from astropy.io import fits
import numpy as np
import glob

def addwl(dir, wl):
    filelist=glob.glob(dir+'/*.fits')
    print('length of file list: ', len(filelist))
    for i in np.arange(len(filelist)):
        hdulist=fits.open(filelist[i])
        header=hdulist[0].header
        header['WLENGTH']=wl
        hdulist[0].header=header
        hdulist.writeto(filelist[i],clobber=True)
