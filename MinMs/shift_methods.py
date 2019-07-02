# numerical python
import numpy as np

# library for various data calculations
import scipy.signal

# for shifting
from scipy.ndimage.interpolation import shift






def centroid(data_arr,xcen,ycen,nhalf=5,derivshift=1.):
    '''
    centroid
    -----------------
    based on dimension-indepdendent line minimization algorithms implemented in IDL cntrd.pro 
    
    inputs
    ----------------
    data_arr      : (matrix of floats) input image
    xcen          : (int) input x-center guess
    ycen          : (int) input y-center guess
    nhalf         : (int, default=5) the excised box of pixels to use
                    recommended to be ~(2/3) FWHM (e.g. only include star pixels).
    derivshift    : (int, default=1) degree of shift used to calculate derivative. 
                     larger values can find shallower slopes more efficiently
    
    
    
    outputs
    ---------------
    xcenf         : the centroided x value
    ycenf         : the centroided y value
    
    dependencies
    ---------------
    numpy         : imported as np
    
    
    
    
    also see another implementation here:
    https://github.com/djones1040/PythonPhot/blob/master/PythonPhot/cntrd.py
    
    '''
    # input image requires the transpose to 
    #
    # find the maximum value near the given point
    data = data_arr[int(ycen-nhalf):int(ycen+nhalf+1),int(xcen-nhalf):int(xcen+nhalf+1)]


    yadjust = nhalf - np.where(data == np.max(data))[0][0]
    xadjust = nhalf - np.where(data == np.max(data))[1][0]
    
    xcen -= xadjust
    ycen -= yadjust
    
    #
    # now use the adjusted centers to find a better square
    data = data_arr[int(ycen-nhalf):int(ycen+nhalf+1),int(xcen-nhalf):int(xcen+nhalf+1)]

    #
    # make a weighting function
    ir = (nhalf-1) > 1 
    
    # sampling abscissa: centers of bins along each of X and Y axes
    nbox = 2*nhalf + 1
    dd = np.arange(nbox-1).astype(int) + 0.5 - nhalf
    
    #Weighting factor W unity in center, 0.5 at end, and linear in between 
    w = 1. - 0.5*(np.abs(dd)-0.5)/(nhalf-0.5) 
    sumc   = np.sum(w)
    
    #
    # fancy comp sci part to find the local maximum
    #
    # this uses line minimization using derivatives
    # (see text such as Press' Numerical Recipes Chapter 10), 
    # treating X and Y dimensions as indepdendent (generally safe for stars). 
    # In this sense the method can be thought of as a two-step gradient descent.
    
    
    #
    # find X centroid
    
    # shift in Y and subtract to get derivative
    deriv = np.roll(data,-1,axis=1) - data.astype(float)
    deriv = deriv[nhalf-ir:nhalf+ir+1,0:nbox-1]
    deriv = np.sum( deriv, 0 )                    #    ;Sum X derivatives over Y direction

    sumd   = np.sum( w*deriv )
    sumxd  = np.sum( w*dd*deriv )
    sumxsq = np.sum( w*dd**2 )
    
    dx = sumxsq*sumd/(sumc*sumxd)
    
    xcenf = xcen - dx
    
    #
    # find Y centroid
    
    # shift in X and subtract to get derivative
    deriv = np.roll(data,-1,axis=0) - data.astype(float)    # Shift in X & subtract to get derivative
    deriv = deriv[0:nbox-1,nhalf-ir:nhalf+ir+1]
    deriv = np.sum( deriv,1 )               #    ;Sum X derivatives over Y direction

    sumd   = np.sum( w*deriv )
    sumxd  = np.sum( w*dd*deriv )
    sumxsq = np.sum( w*dd**2 )
    
    dy = sumxsq*sumd/(sumc*sumxd)
    
    ycenf = ycen - dy
    
    return xcenf,ycenf





def cross_image(im1, im2, **kwargs):
    '''
    cross_image
    ---------------
    calcuate cross-correlation of two images in order to find shifts
    
    
    inputs
    ---------------
    im1                      : (matrix of floats)  first input image
    im2                      : (matrix of floats) second input image
    boxsize                  : (integer, optional) subregion of image to cross-correlate
    
    
    returns
    ---------------
    xshift                   : (float) x-shift in pixels
    yshift                   : (float) y-shift in pixels
    
    
    dependencies
    ---------------
    scipy.signal.fftconvolve : two-dimensional fourier convolution
    centroid                 : a centroiding algorithm of your choosing or defintion
    numpy                    : imported as np
    
    todo
    ---------------
    -add more **kwargs capabilities for centroid argument
    
    '''
    
    # The type cast into 'float' is to avoid overflows:
    im1_gray = im1.astype('float')
    im2_gray = im2.astype('float')

    # Enable a trimming capability using keyword argument option.
    if 'boxsize' in kwargs:
        im1_gray = im1_gray[0:kwargs['boxsize'],0:kwargs['boxsize']]
        im2_gray = im2_gray[0:kwargs['boxsize'],0:kwargs['boxsize']]

    # Subtract the averages (means) of im1_gray and im2_gray from their respective arrays     
    im1_gray -= np.nanmean(im1_gray)
    im2_gray -= np.nanmean(im2_gray)
    
    # guard against extra nan values
    im1_gray[np.isnan(im1_gray)] = np.nanmedian(im1_gray)
    im2_gray[np.isnan(im2_gray)] = np.nanmedian(im2_gray)


    # Calculate the correlation image using fast Fourrier Transform (FFT)
    # Note the flipping of one of the images (the [::-1]) to act as a high-pass filter
    corr_image = scipy.signal.fftconvolve(im1_gray, im2_gray[::-1,::-1], mode='same')

    # Find the peak signal position in the cross-correlation, which gives the shift between the images
    corr_tuple = np.unravel_index(np.nanargmax(corr_image), corr_image.shape)
    
    try: # try to use a centroiding algoritm to find a better peak
        xcenc,ycenc = centroid(corr_image.T,corr_tuple[0],corr_tuple[1],nhalf=10,derivshift=1.)

    except: # if centroiding algorithm fails, use just peak pixel
        xcenc,ycenc = corr_tuple
        
    # Calculate shifts (distance from central pixel of cross-correlated image)
    xshift = xcenc - corr_image.shape[0]/2.
    yshift = ycenc - corr_image.shape[1]/2.

    return xshift,yshift


def shift_image(image,xshift,yshift):
    '''
    shift_image
    -------------
    wrapper for scipy's implementation that shifts images according to values from cross_image
    
    inputs
    ------------
    image           : (matrix of floats) image to be shifted
    xshift          : (float) x-shift in pixels
    yshift          : (float) y-shift in pixels
    
    outputs
    ------------
    shifted image   : shifted, interpolated image. 
                      same shape as input image, with zeros filled where the image is rolled over
    
    
    '''
    return scipy.ndimage.interpolation.shift(image,(xshift,yshift))


