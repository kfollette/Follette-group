from astropy.io import fits
import numpy as np
import glob
import gaussfitter


def addwl(dir, wl):
    filelist=glob.glob(dir+'/*.fits')
    print('length of file list: ', len(filelist))
    for i in np.arange(len(filelist)):
        hdulist=fits.open(filelist[i])
        header=hdulist[0].header
        header['WLENGTH']=wl
        hdulist[0].header=header
        hdulist.writeto(filelist[i],overwrite=True)

def addstarpeak(dir,amp, fwhm, debug=False):
    filelist = glob.glob(dir + '/*.fits')
    filelist = sorted(filelist, key=len)
    dummy_im = fits.getdata(filelist[0])
    size = dummy_im.shape
    xcen = int((size[0]-1)/2)
    ycen = int((size[1]-1)/2)
    # Returns (background height, amplitude, x, y, width_x, width_y, rotation angle)
    input_parameters = [0, amp, xcen, ycen, fwhm, fwhm, 0]

    dummyim=fits.getdata(filelist[0])
    dim=dummyim.shape[1]
    cen = int((dim-1)/2)
    diff=np.zeros(len(filelist))

    for i in np.arange(len(filelist)):
        im = fits.getdata(filelist[i])
        head = fits.getheader(filelist[i])
        p = gaussfitter.moments(im, circle=False, rotate=False, vheight=False, estimator=numpy.ma.median)
        head['STARPEAK']=p[0]
        fits.writeto(filelist[i], im, header=head, overwrite=True)

        if debug==True:
            #print(filelist[i])
            #print('fit peak is:', p[0], '. max pixel is: ', np.nanmax(im[cen-10:cen+10,cen-10:cen+10]))
            diff[i] = p[0]-np.nanmax(im[cen-10:cen+10,cen-10:cen+10])

    if debug==True:
        print('standard deviation of differences is: ', np.std(diff))
        print('max difference is:', np.max(abs(diff)))
        print('median difference is:', np.median(diff))

def addradec(dir,ra, dec):
    filelist = glob.glob(dir + '/*.fits')
    for i in np.arange(len(filelist)):
        im = fits.getdata(filelist[i])
        head = fits.getheader(filelist[i])
        head['RA'] = ra
        head['DEC'] = dec
        fits.writeto(filelist[i], im, header=head, overwrite=True)