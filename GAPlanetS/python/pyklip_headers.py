from astropy.io import fits
import numpy as np
import glob
from astropy.modeling import models, fitting, functional_models
import pdb


def addstarpeak(dir, debug=False, mask=False, ghost=False, wl='Line'):
    """
    fits star (or ghost) and adds STARPEAK to header. needed for PyKLIP.
    :param dir: directory
    :param amp: amplitude guess
    :param fwhm: fwhm
    :param debug: will print comparison of peak pixel values to fit values if set
    :param mask: will mask NaNs with zeros for fitting if set
    :param ghost: will fit ghost in lieu of star and return ghost peak and estimate for star peak in header.
    :param wl: set if ghost is set so it knows which scale factor to pull. values are "Line" or "Cont"
    :return: list of star (or ghost) peaks
    """
    filelist = glob.glob(dir + '/*.fits')
    ##sort sequentially
    filelist.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
    dummy_im = fits.getdata(filelist[0])
    size = dummy_im.shape
    xcen = int((size[0] - 1) / 2)
    ycen = int((size[1] - 1) / 2)
    # guesses for fit parameters. used as starting point
    # input_parameters = [0, amp, xcen, ycen, fwhm, fwhm, 0]
    if ghost == True:
        xcen = int(xcen + 157.5)
        ycen = int(ycen - 7)
        if wl == "Line":
            ghost_scale = 184.7
        if wl == "Cont":
            ghost_scale = 193.6
        if wl == "False":
            "wl keyword must be set for ghost calibration"
            return
    diff = np.zeros(len(filelist))
    peaks = []
    fwhmlist = []
    #size of image stamp
    width = 31
    stmpsz = int((width - 1) / 2)

    for i in np.arange(len(filelist)):
        im = fits.getdata(filelist[i])
        # make a copy - only invoked in case of ghost=True, but can't think of a cleverer way
        imcopy = np.copy(im)
        head = fits.getheader(filelist[i])
        #crop around star (or ghost) for fitting
        imcopy = imcopy[ycen - stmpsz-1:ycen + stmpsz, xcen - stmpsz-1:xcen + stmpsz]

        #set up fit
        y,x = np.mgrid[:width,:width]
        #Moffat PSF model
        g_init=models.Moffat2D(np.nanmax(imcopy), stmpsz, stmpsz, 6, 1)
        fit_g=fitting.LevMarLSQFitter()

        #if ghost == True:
            # CROP AROUND GHOST
            #im = im[ycen - 50:ycen + 50 + 1, xcen - 50:xcen + 50 + 1]

        # do fit
        p=fit_g(g_init,x,y,imcopy)

        # populate headers for each individual image
        if ghost == True:
            head['GSTPEAK'] = p.amplitude.value
            head['STARPEAK'] = p.amplitude.value * ghost_scale
            if p.fwhm is not np.nan:
                head['FWHM'] = p.fwhm
            else:
                head['FWHM'] = np.median(fwhmlist)
        else:
            head['STARPEAK'] = p.amplitude.value
            if p.fwhm is not np.nan:
                head['FWHM'] = p.fwhm
            else:
                head['FWHM'] = np.median(fwhmlist)

        # record peak
        peaks.append(p.amplitude.value)

        # record fwhm
        fwhmlist.append(p.fwhm)


        # print a warning if any peak values are unphysical
        if (p.amplitude.value < 0) or (p.amplitude.value > 17000) or (np.isnan(p.amplitude.value)) == True:
            print("warning: unphysical peak value of", p.amplitude.value, 'for image', i + 1)

        # write out file with peak info in header
        if ghost == True:
            fits.writeto(filelist[i], im, header=head, overwrite=True)
        else:
            fits.writeto(filelist[i], im, header=head, overwrite=True)

        if debug == True:
            # print(filelist[i])
            # print('fit peak is:', p[0], '. max pixel is: ', np.nanmax(im[cen-10:cen+10,cen-10:cen+10]))
            imsz = im.shape[1]
            imcen = int((imsz - 1) / 2.)
            diff[i] = p.amplitude.value - np.nanmax(im[imcen - 10:imcen + 10, imcen - 10:imcen + 10])

    # write out list of peaks one directory up so KLIP doesn't try to pull it
    if ghost == True:
        fits.writeto(dir + '../' + str(wl) + 'ghostpeaks.fits', np.array(peaks), overwrite=True)
        fits.writeto(dir + '../' + str(wl) + 'fwhmlist.fits', np.array(fwhmlist), overwrite=True)
    else:
        fits.writeto(dir + '../' + str(wl) + 'starpeaks.fits', np.array(peaks), overwrite=True)
        fits.writeto(dir + '../' + str(wl) + 'fwhmlist.fits', np.array(fwhmlist), overwrite=True)

    if debug == True:
        print('standard deviation of difference between fit peak and max pixel is: ', np.std(diff))
        print('max difference is:', np.max(abs(diff)))
        print('median difference is:', np.median(diff))
    return (peaks)


def addradec(dir, ra, dec):
    """
    adds ra and dec info to headers in cases where wasn't propagated through pipeline
    should be obsolete with pipeline bug fixes, but maintaining for posterity
    """
    filelist = glob.glob(dir + '/*.fits')
    for i in np.arange(len(filelist)):
        im = fits.getdata(filelist[i])
        head = fits.getheader(filelist[i])
        head['RA'] = ra
        head['DEC'] = dec
        fits.writeto(filelist[i], im, header=head, overwrite=True)


def addwl(dir, wl):
    """
    adds wavelength in microns to headers in cases where wasn't propagated through pipeline
    should be obsolete with pipeline bug fixes, but maintaining for posterity
    """
    filelist = glob.glob(dir + '/*.fits')
    print('length of file list: ', len(filelist))
    for i in np.arange(len(filelist)):
        hdulist = fits.open(filelist[i])
        header = hdulist[0].header
        header['WLENGTH'] = wl
        hdulist[0].header = header
        hdulist.writeto(filelist[i], overwrite=True)
