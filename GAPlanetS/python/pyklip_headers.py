from astropy.io import fits
import numpy as np
import glob
import gaussfitter


def addwl(dir, wl):
    filelist=glob.glob(dir+'/*.fits')
    print('length of file list: ', len(filelist))
    for i in np.arange(len(filelist)):
        hdulist=fits.open(filelist[i])
        header=hdulist[0].header
        header['WLENGTH']=wl
        hdulist[0].header=header
        hdulist.writeto(filelist[i],clobber=True)

def addstarpeak(dir,amp, fwhm):
    filelist = glob.glob(dir + '/*.fits')
    dummy_im = fits.getdata(filelist[0])
    size = dummy_im.shape
    xcen = int((size[0]-1)/2)
    ycen = int((size[1]-1)/2)
    # Returns (background height, amplitude, x, y, width_x, width_y, rotation angle)
    input_parameters = [0, amp, 25, 25, fwhm, fwhm, 0]

    for i in np.arange(len(filelist)):
        im = fits.getdata(filelist[i])
        head = fits.getheader(filelist[i])
        p = gaussfitter.moments(im, circle=False, rotate=False, vheight=False)
        head['STARPEAK']=p[0]
        fits.writeto(filelist[i], im, header=head, clobber=True)

