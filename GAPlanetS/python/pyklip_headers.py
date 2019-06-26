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

def addstarpeak(dir,amp, fwhm, debug=False, mask=False):
    filelist = glob.glob(dir + '/*.fits')
    ##sort sequentially
    filelist.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
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
	peaks=[]

    for i in np.arange(len(filelist)):
        im = fits.getdata(filelist[i])
        head = fits.getheader(filelist[i])
        if mask==True:
        	#replace nans in image with zeros, just for the purposes of fitting (median should be robust)
        	maskedim=np.copy(im)
        	maskedim[np.isnan(im)==True]=0.
        	#fits.writeto('test.fits',maskedim,overwrite=True)
        	p = gaussfitter.moments(maskedim, circle=False, rotate=False, vheight=False, estimator=np.nanmedian)
        else:
        	p = gaussfitter.moments(im, circle=False, rotate=False, vheight=False)
        #print(p)
        head['STARPEAK']=p[0]
        peaks.append(p[0])
        fits.writeto(filelist[i], im, header=head, overwrite=True)

        if debug==True:
            #print(filelist[i])
            #print('fit peak is:', p[0], '. max pixel is: ', np.nanmax(im[cen-10:cen+10,cen-10:cen+10]))
            diff[i] = p[0]-np.nanmax(im[cen-10:cen+10,cen-10:cen+10])
	
	fits.writeto('peaks.fits', peaks, overwrite=True)
	
    if debug==True:
        print('standard deviation of difference between fit peak and max pixel is: ', np.std(diff))
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