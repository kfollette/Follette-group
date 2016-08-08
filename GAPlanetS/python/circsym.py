from astropy.io import fits
import numpy as np
import math
import scipy

from scipy.ndimage.filters import gaussian_filter
import scipy.ndimage.interpolation as sci
import scipy.optimize as optimize


"""
NAME: visao_circlesym

PURPOSE: Finds center of circular symmetry of median combinations of registered images, and shifts the image cube
        to this center

DESCRIPTION: Opens the image cubes created by visao_reg and takes a median of the cube. If mask is specified, creates
            a mask of ones and zeros to block out saturated pixels near the center of the star. Center_circlesym is
            given the median image, one dimensional x and y arrays populated with values starting at half the
            dimensions of the median image, the pixel radius for which to search, and if specified a mask blocking
            out the saturated images. An x and y shift are returned and each image in the image cube is shifted by
            those values. Then the image cube is written out. This is done for both the line and the continuum images
            and then if sdi is specified, each image in the sdi cube is shifted by the save amount as the continuum.

OPTIONAL *ARGS:
    'flat':     string, specifies if flat is included in the name of the image cube that will be opened
    'sdi':      string, specified if a spectral differential image cube has been created

OPTIONAL KWARGS:
    clip:       integer, dimensions of each image in the image cube
    rmax:       integer, maximum radius at which to attempt to center the star
    msk:        integer, radius of the central mask to be constructed to block out the star in saturated images

OUTPUTS:    Does not return anything
    'Line_???_reg_circsym.fits':    image cube of the line shifted and centered images
    'Cont_???_reg_circsym.fits':    image cube of the continuum shifted and centered images

HISTORY:
    Created: 2016 by Kate Follette
        Modified: 2016-06-30 by Follette: does line and continuum separately to save memory
    Py Trans:   2016-08-05 by Wyatt Mullen, wmullen1@stanford.edu

"""
def visao_circlesym(*args, **kwargs):
    if 'flat' in args: namestr = '_flat_'
    else: namestr = '_'

    if 'clip' in kwargs:
        clip = kwargs['clip']
        Line = fits.open('Line_clip' + str('%03d' % (clip,)) + '_reg.fits')[0]
    else:
        Line = fits.open('Line_reg.fits')[0]
    Linemed = np.nanmedian(Line.data, axis=0) #median of the image cube
    nims = Line.shape[0]

    dimy = Linemed.shape[0]
    dimx = Linemed.shape[1]
    if 'msk' in kwargs:
        kwargs['msk'] = msk
        print('Constructing mask with radius of ' + str(msk) + ' pixels.')
        mask = mkmask(dimx, dimy, msk)

    if 'rmax' in kwargs: lim = kwargs['rmax']
    else: lim = dimx / 2

    xr = np.arange(dimx // 2 + 1.0) - dimx // 4
    yr = np.arange(dimy // 2 + 1.0) - dimy // 4
    print('Calculating center of circular symmetry for median line image.')

    #main circle sym function call, returns center of star
    if 'msk' in kwargs:
        Line_xc, Line_yc = center_circlesym(Linemed, xr, yr, lim, mask=mask)
    else:
        Line_xc, Line_yc = center_circlesym(Linemed, xr, yr, lim)
    Line_shift = [(dimy - 1) / 2 - Line_yc, (dimx - 1) / 2 - Line_xc]

    print('Center of circular symmetry for median Line image is (' + str(Line_xc) + ',' + str(Line_yc) + ')')
    print('Shifting all Line images by ' + str(Line_shift))
    Line_cent = np.zeros(Line.shape)
    for i in range(0,nims):
        print('Shifting Line image ' + str(i + 1) + ' of ' + str(nims))
        sci.shift(Line.data[i,:,:], Line_shift, order=1, output=Line_cent[i,:,:]) # shifting the line images

    print('Writing centered line image cube')
    if 'clip' in kwargs:
        fits.writeto(('Line_clip' + str('%03d' % (clip,)) + '_reg_circsym.fits'), Line_cent, clobber=True)
    else:
        fits.writeto('Line_reg_circsym.fits', Line_cent, clobber=True)
    del Line_cent
    del Linemed

    ###################################### Start Over With Continuum ############################################

    if 'clip' in kwargs:
        clip = kwargs['clip']
        Cont=fits.open('Cont_clip' + str('%03d' % (clip,)) + '_reg.fits')[0] #opening image cube
    else:
        Cont = fits.open('Cont_reg.fits')[0]

    Contmed = np.nanmedian(Cont.data, axis=0) #takes median of image cube
    if 'msk' in kwargs:
        kwargs['msk'] = msk
        print('Constructing mask with radius of ' + str(msk) + ' pixels.')
        mask = mkmask(dimx, dimy, msk)

    xr = np.arange(dimx // 2 + 1.0) - dimx / 4 #initializing xr and yr with values
    yr = np.arange(dimy // 2 + 1.0) - dimy / 4
    print('Calculating center of circular symmetry for median line image.')

    #main center_circlesym function call, returns the center of the image
    if 'msk' in kwargs:
        Cont_xc, Cont_yc = center_circlesym(Contmed, xr, yr, lim, mask=mask)
    else:
        Cont_xc, Cont_yc = center_circlesym(Contmed, xr, yr, lim)
    Cont_shift = [(dimy - 1) / 2 - Cont_yc, (dimx - 1) / 2 - Cont_xc]
    print('Center of circular symmetry for median Continuum image is (' +  str(Cont_xc) + ',' + str(Cont_yc) + ')')
    print('Shifting all continuum images by' + str(Cont_shift)) #the shift reads as y,x

    Cont_cent = np.zeros(Cont.shape)
    for i in range(0,nims):
        print('Shifting continuum image ' + str(i + 1) + ' of ' + str(nims))
        sci.shift(Cont.data[i,:,:], Cont_shift, order=1, output=Cont_cent[i,:,:]) #shift each image

    print('Writing centered continuum image cube.')
    if 'clip' in kwargs:
        fits.writeto(('Cont_clip' + str('%03d' % (clip,)) + '_reg_circsym.fits'), Cont_cent, clobber=True)
    else:
        fits.writeto('Cont_reg_circsym.fits', Cont_cent, clobber=True)
    del Cont_cent #clear memory for use
    del Contmed

    ################################# Shift SDI ###########################################

    if 'sdi' in args:
        if 'clip' in kwargs:
            SDI1 = fits.open('SDI_sc' + str('%.2f' % clip) + '_clip' + str('%03d' % (clip,))
                             + str(namestr) + 'reg.fits')[0]
        else:
            SDI1 = fits.open('SDI_sc' + str('%.2f' % clip) + str(namestr) + 'reg.fits')[0]
        SDI1_cent = np.zeros((nims, dimy, dimx))

        for i in range(0,nims):
            print('Shifting SDI image ' + str(i+1) + ' of ' + str(nims))
            sci.shift(SDI1[i,:,:], Cont_shift, order=1, output=SDI1_cent[i,:,:])

        print('Writing centered SDI image cube')
        if 'clip' in kwargs:
            fits.writeto('SDI_sc' + str('%.2f' % clip) + '_clip' + str('%03d' % (clip,)) + 'reg_circsym.fits',
                         SDI1_cent, clobber=True)
        else:
            fits.writeto('SDI_sc' + str('%.2f' % clip) + 'reg_circsym.fits', SDI1_cent, clobber=True)



"""
NAME: center_circlesym

PURPOSE: Finds the center of a star assuming circular symmetry

DESCRIPTION: Iterates through the vectors of x and y offsets from the image and calculates the distance from the
            center for each pixel using rarr. Then a grid is gradually built up where pixels near the center of
            the star have the fewest pixel counts and pixels far away have the highest. This occurs in a radius
            around the star specified by rmax. The location of the minimum pixel is found in the grid and gcntrd
            is called to find the center of the grid based on a gaussian fit. X and Y coordinates are returned,
            which then are correlated with the actual coordinates of the image to provide the center which are
            returned.

INPUTS:
    im:     2-d grid from the median fits image from visao_circlesym
    xr/yr:  vectors of x and y offsets from the center of the image
    rmax:   maximum radius to consider when aligning pixels/calculating the center

OPTIONAL KWARGS:
    mask:   optional 1/0 (or maybe true/false) mask, with 0 (false) specifying pixels to ignore

OUTPUTS:
    xc/yc:     The calculated centroid of the star

HISTORY:
    Created:    ???
    Py Trans:   2016-08-05 by Wyatt Mullen, wmullen1@stanford.edu
"""
def center_circlesym(im, xr, yr, rmax, **kwargs):
    print('Starting center_circlesym')
    dimy = im.shape[0]
    dimx = im.shape[1]
    grid = np.zeros((len(yr),len(xr)))

    print('Out of ' + str(len(xr)) + ' rows, this many have finished: ')
    for i in range(0,len(xr)):

        print(str(i) + ', ', end='', flush=True) #print statement to show progress
        for j in range(0,len(yr)):

            #calling rarr to figure out distance from center
            r, x, y = rarr(dimx, dimy, 'pixels', xcen=xr[i], ycen=yr[j])
            for k in range(0, rmax+1):
                if ('mask' in kwargs) and (kwargs['mask'].shape == im.shape):
                    idx = np.where((r >= k) & (r < k+1) & (kwargs['mask'] > 0))
                else:
                    idx = np.where((r >= k) & (r < k+1))

                if len(idx[0]) != 0:
                    sd = np.nanstd(im[idx])
                    if not np.isfinite(sd): sd = 0
                    grid[j][i] = grid[j][i] + sd / np.fabs(np.nanmedian(im[idx]))

    minPos = np.argmin(grid)
    pos = np.unravel_index(minPos, grid.shape) #finding pos of min element of grid

    #fits.writeto('Circsym_grid_cont.fits', grid, clobber=True)
    xcc, ycc = gcntrd(-1 * grid, pos[1], pos[0], 0.5 * len(xr)) #calculating centroid for small grid

    #xcc = 112.490338
    #ycc = 112.51506884
    xc = xr[0] + xcc * (xr[1] - xr[0]) + 0.5 * (dimx - 1) #converting position from small grid to actual image
    yc = yr[0] + ycc * (yr[1] - yr[0]) + 0.5 * (dimy - 1)
    print('Finishing center_circlesym')
    return (xc, yc)

"""
NAME: gcntrd

PURPOSE: Compute the stellar centroid by Gaussian fits to marginal X,Y, sums

DESCRIPTION: GCNTRD uses the DAOPHOT "FIND" centroid algorithm by fitting Gaussians to the marginal X,Y distributions.
            User can specify bad pixels (either by using the MAXGOOD keyword or setting them to NaN) to be ignored in
            the fit. Pixel values are weighted toward the center to avoid contamination by neighboring stars.

INPUTS:
    img:    two dimensional image array
    x,y:    scalar or vector integers giving approximate stellar center
    fwhm:   floating scalar full width half max to compute centroid, centroid computed withing box of half width
            equal to 1.5 sigma = 0.637 * FWHM

OUTPUTS:
    xcen:   computed x centroid position, same number of points as x
    ycen:   computed y centroid position
            (values for xcen and ycen will not be computed if the computed centroid falls outside of the box, or if
            there are too many bad pixels, or if the best-fit Gaussian has a negative height. If the centroid cannot be
            computed, then a message is displayed (unless 'silent' is set) and xcen and ycen are set to -1.

OPTIONAL ARGS:
    silent:         gcntrd will not print error message if centroid cannot be found
    keepcenter:

OPTIONAL KWARGS:
    maxgood - int:  only pixels with values leess than maxgood are used to determine centroid

HISTORY:
    Created: 2004-06, W. Landsman  following algorithm used by P. Stetson in DAOPHOT2.
        Modified:   2008-03, W. Landsman to allow shifts of more than 1 pixel from initial guess
                    2009-01, W. Landsman to perform Gaussian convolution first before finding max pixel to smooth noise
    Py Trans: 2016-08 by Wyatt Mullen, wmullen1@stanford.edu
"""
def gcntrd(img, x, y, fwhm, *args, **kwargs):
    print('Starting gcntrd')
    sz_image = img.shape
    if len(sz_image) != 2: raise TypeError('Invalid dimensions - image array must be 2 dimensional')
    xsize = sz_image[1]
    ysize = sz_image[0]

    if isinstance(x, int) or isinstance(x, float):
        npts = 1
    else:
        npts = len(x)

    maxbox = 13 #why this value?
    radius = max(0.637 * fwhm, 2.001)
    radsq = radius**2
    sigsq = (fwhm / 2.35482)**2
    nhalf = min(int(radius), (maxbox - 1) // 2)
    nbox = 2 * nhalf + 1 # of pix in side of convolution box

    ix = np.array([round(x)])
    iy = np.array([round(y)])

    g = np.zeros((nbox, nbox))
    row2 = (np.arange(nbox, dtype=float) - nhalf)**2
    g[nhalf,:] = row2
    for i in range(1, nhalf+1):
        temp = row2 + i**2
        g[nhalf-i,:] = temp
        g[nhalf+i,:] = temp

    mask = g <= radsq
    good = np.where(mask)
    pixels = len(good[0])
    g = math.e**(-0.5*g/sigsq)

    """ In fitting Gaussians to the marginal sums, pixels will arbitrarily be
    assigned weights ranging from unity at the corners of the box to
    NHALF^2 at the center (e.g. if NBOX = 5 or 7, the weights will be

                                     1   2   3   4   3   2   1
          1   2   3   2   1          2   4   6   8   6   4   2
          2   4   6   4   2          3   6   9  12   9   6   3
          3   6   9   6   3          4   8  12  16  12   8   4
          2   4   6   4   2          3   6   9  12   9   6   3
          1   2   3   2   1          2   4   6   8   6   4   2
                                     1   2   3   4   3   2   1

    respectively).  This is done to desensitize the derived parameters to
    possible neighboring, brighter stars. """

    x_wt = np.zeros((nbox, nbox))
    wt = nhalf - np.fabs(np.arange(nbox) - nhalf) + 1
    for i in range(0, nbox):
        x_wt[i,:] = wt
    y_wt = np.transpose(x_wt)
    pos = str(x) + ' ' + str(y)

    if 'keepcenter' not in args:
        c = g*mask
        sumc = np.sum(c)
        sumcsq = np.sum(c**2) - sumc**2 / pixels
        sumc = sumc / pixels
        c[good] = (c[good] - sumc) / sumcsq

    xcen, ycen = [], []
    for i in range(0,npts):
        if 'keepcenter' not in args:
            if (ix[i] < nhalf) or ((ix[i] + nhalf) > xsize - 1) or (iy[i] < nhalf) or ((iy[i] + nhalf) > ysize-1):
                if 'silent' not in args:
                    raise RuntimeError('Position ' + str(pos[i]) + ' is too near edge of image.')
            x1 = max(ix[i] - nbox, 0)
            x2 = min(ix[i] + nbox, xsize - 1)
            y1 = max(iy[i] - nbox, 0)
            y2 = min(iy[i] + nbox, ysize - 1)
            h = img[y1:y2 + 1, x1:x2 + 1]

            h = scipy.ndimage.convolve(h, c)
            h = h[nbox - nhalf: nbox + nhalf + 1, nbox - nhalf: nbox + nhalf + 1]
            d = img[iy[i] - nhalf: iy[i] + nhalf + 1, ix[i] - nhalf: ix[i] + nhalf + 1]

            if 'maxgood' in kwargs:
                ig = np.where(d < maxgood)
                mx = np.nanmax(d[ig])

            mx = np.nanmax(h) #max pix val in bigbox
            mx_pos = np.where(h == mx) #num pix w/max val
            idx = mx_pos[1] % nbox #x coord of max pix
            idy = mx_pos[0] % nbox #y coord of max pix
            if len(mx_pos[0]) > 1:
                idx = round(np.sum(idx) / len(idx))
                idy = round(np.sum(idy) / len(idy))
            else:
                idx = idx
                idy = idy

            xmax = ix[i] - nhalf + idx #x coord in original image array
            ymax = iy[i] - nhalf + idy #y coord in original image array
        else: #if keepcenter is specified
            xmax = ix[i]
            ymax = iy[i]

        ########################################################################################
        #check *new* center location for range
        #added by Hogg

        if (xmax < nhalf) or ((xmax + nhalf) > xsize-1) or (ymax < nhalf) or ((ymax + nhalf) > ysize-1):
            if 'silent' not in args:
                raise RuntimeError('Position ' + str(pos[i]) + ' is too near edge of image.')
        #########################################################################################

        d = img[ymax - nhalf: ymax + nhalf + 1, xmax - nhalf: xmax + nhalf + 1]
        #extract subimage centered on max pixel, skipping debugging
        if 'maxgood' in kwargs:
            mask = (d < maxgood) #if we need values for this we should use np.amin(d, maxgood)
        else: # isinstance(img[0][0], float):
            mask = np.isfinite(d)
            #mask = np.zeros(d.shape)
        #else:
            #mask = np.ones(nbox, nbox)
        maskx = np.sum(mask, 1) > 0
        masky = np.sum(mask, 0) > 0

        #at least 3 points are needed in the partial sum to compute gaussian
        if np.sum(maskx) < 3 or np.sum(masky) < 3:
            if 'silent' not in args:
                raise RuntimeError('Position ' + str(pos[i]) + ' has insufficient good points')

        ywt = y_wt * mask
        xwt = x_wt * mask
        wt1 = wt * maskx
        wt2 = wt * masky
        """ Centroid computation:   The centroid computation was modified in Mar 2008 and now differs from DAOPHOT
        which multiplies the correction dx by 1/(1+abs(dx)). The DAOPHOT method is more robust (e.g. two different
        sources will not merge) especially in a package where the centroid will be subsequently be redetermined
        using PSF fitting.   However, it is less accurate, and introduces biases in the centroid histogram.
        The change here is the same made in the IRAF DAOFIND routine
        (see http://iraf.net/article.php?story=7211&query=daofind)"""

        #computation for x centroid
        sd = np.nansum(d * ywt, 0)
        sg = np.nansum(g * ywt, 0)
        sumg = np.nansum(wt1 * sg)
        sumgsq = np.nansum(wt1 * sg * sg)
        sumgd = np.nansum(wt1 * sg * sd)
        sumd = np.nansum(wt1 * sd)
        p = np.nansum(wt1)
        xvec = nhalf - np.arange(nbox)
        dgdx = sg * xvec
        sdgdxs = np.nansum(wt1 * dgdx ** 2)
        sdgdx = np.nansum(wt1 * dgdx)
        sddgdx = np.nansum(wt1 * sd * dgdx)
        sgdgdx = np.nansum(wt1 * sg * dgdx)

        #height of the best-fitting marginal Gaussian. If this is not positive then centroid will not make sense
        hx = (sumgd - sumg * sumd / p) / (sumgsq -sumg**2 / p)
        if hx <= 0:
            if 'silent' not in args:
                raise RuntimeError('Position ' + str(pos[i]) + ' cannot be fit by a Gaussian')

        skylvl = (sumd - hx * sumg) / p
        dx = (sgdgdx - (sddgdx - sdgdx * (hx * sumg + skylvl * p))) / (hx * sdgdxs / sigsq)
        if math.fabs(dx) >= nhalf:
            if 'silent' not in args:
                raise RuntimeError('Position ' + str(pos[i]) + ' is too far from initial guess')

        xcen.append(xmax + dx) #x centroid in original array

        # computation for y centroid
        sd = np.nansum(d * xwt, 1)
        sg = np.nansum(g * xwt, 1)
        sumg = np.nansum(wt2 * sg)
        sumgsq = np.nansum(wt2 * sg * sg)
        sumgd = np.nansum(wt2 * sg * sd)
        sumd = np.nansum(wt2 * sd)
        p = np.nansum(wt2)
        yvec = nhalf - np.arange(nbox)
        dgdy = sg * yvec
        sdgdys = np.nansum(wt1 * dgdy ** 2)
        sdgdy = np.nansum(wt1 * dgdy)
        sddgdy = np.nansum(wt1 * sd * dgdy)
        sgdgdy = np.nansum(wt1 * sg * dgdy)

        hy = (sumgd - sumg * sumd / p) / (sumgsq - sumg ** 2 / p)
        if hy <= 0:
            if 'silent' not in args:
                raise RuntimeError('Position ' + str(pos[i]) + ' cannot be fit by a Gaussian')

        skylvl = (sumd - hy * sumg) / p
        dy = (sgdgdy - (sddgdy - sdgdy * (hy * sumg+ skylvl * p))) / (hy * sdgdys / sigsq)
        if math.fabs(dy) >= nhalf:
            if 'silent' not in args:
                raise RuntimeError('Position ' + str(pos[i]) + ' is too far from initial guess')

        ycen.append(ymax + dy) # y centroid in original array

    print('Finishing gcntrd')
    return (xcen, ycen)


"""
NAME: rarr

PURPOSE: Create a 2D array, dimx*dimy, with value at each point being the radius of that pixel from the center.
        The default center is x=0.5*(dimx-1), y=0.5*(dimy-1), but this can be changed. The radius is specified
        relative to dimx/2, that is a pixel a distance dimx/2 from the center will have radius 1.0.  You can
        have the radius in pixels instead by setting the /pixels keyword.  If desired, will also calculate the
        angle at each pixel.

INPUTS:
    dimx:   the x dimension
    dimy:   the y dimension

OPTIONAL MAP KEYS:
    xcen:   the desired x center offset, in units of pixels relative to 0.5*(dimx-1)
    ycen:   the desired y center offset, in units of pixels relative to 0.5*(dimy-1)
    pixels: output r, x, and y in units of pixels

OUTPUTS:
    r:  the returned radius array (goes from 0 to 1 along x-axis, unless pixels is set)
    x:  the x coordinate of each pixel (goes from -1 to 1 unless pixels is set
    y:  the y coordinate of each pixel (same as x, unless rectangular)
    theta:  (optional) set to the angle in radians of each pixel

HISTORY:
    Created:    Jared Males, jrmales@email.arizona.edu
        Updated:    2012-12-27 by Jared Males (updated/documented)
                    2013-03-14 by Jared Males (added pixels keyword)
    Py Trans:   2016-07-15 by Wyatt Mullen, wmullen1@stanford.edu
"""
def rarr(dimx, dimy, *args, **kwargs):
    x = np.array([(0.5 + i) for i in range(0,dimx)])
    x = x / dimx * 2.0 - 1.0
    x = np.array([x for i in range(0,dimy)])
    y = np.array([(0.5 + i) for i in range(0,dimy)])
    y = (y / dimy * 2.0 - 1.0) * (dimy / dimx)
    y = np.transpose(np.array([y for i in range(0,dimx)]))

    if 'xcen' in kwargs:
        xcen = kwargs['xcen']
        xc = 2 * xcen / dimx
        x = x - xc
    if 'ycen' in kwargs:
        ycen = kwargs['ycen']
        yc = 2 * ycen / dimy
        y = y - yc

    r = np.sqrt(x**2 + y**2)
    if 'theta' in args: theta = math.atan(y,x)

    if 'pixels' in args:
        r = 0.5 * r * dimx
        x = 0.5 * x * dimx
        y = 0.5 * y * dimx #shouldn't this be dimy????

    if 'theta' in args:
        return (r, x, y, theta)
    else:
        return (r, x, y)



"""
NAME: mkmask

PURPOSE: Creates an image mask of ones and zeros or ones and nans

DESCRIPTION: Takes image dimensions (xdim, ydim) and creates a 2 dimensional grid of zeros and ones used for a mask of
            a specified radius. The mask is circular and centralized (unless specified) with pixels of zeros under the
            mask and ones outside. If specified the mask can be made of nans or the mask can be inverted.

INPUTS:
    xdim/ydim:  dimensions of the 2-d mask array to be created
    r:          radius of the circular mask in pixels

    OPTIONAL KEYWORDS:
        'fits':     writes mask to fits file
        'reverse':  switches ones and zeros so that mask is outside of circle
        'NaN':      makes mask out of NaNs instead of zeros
    OPTIONAL MAP KEYS:
        'cen' = (xpos,ypos):    coordinates of alternate center

OUTPUTS:
    mask:      2d image mask array

HISTORY:
    CREATED: 2015 by Kate Follette, kbf@stanford.edu
        modified 2016-03-29 to allow off center mask and reversed options
        modified 2016-06-07 to add nan option
    PY TRANS: 2016-07-18 by Wyatt Mullen, wmullen1@stanford.edu

"""
def mkmask(xdim, ydim, r, *args, **kwargs):
    print('Starting mkmask')
    mask = np.ones((ydim, xdim))
    for x in range(0,ydim):
        for y in range(0,xdim):
            if 'cen' in kwargs:
                cen = kwargs['cen']
                if math.sqrt((cen[0] - x)**2 + (cen[1] - y)**2) < r: mask[y,x] = 0
            else:
                if math.sqrt(((xdim - 1) / 2 - x)**2 + ((ydim - 1) / 2 - y)**2) < r: mask[y,x] = 0

    #inverting mask
    if 'reverse' in args:
        ones = np.where(mask != 0)
        zeros = np.where(mask == 0)
        mask[ones] = 1
        mask[zeros] = 0
    if 'NaN' in args:
        zeros = np.where(mask == 0)
        mask[zeros] = np.NaN
    if 'fits' in args:
        fits.writeto('mask.fits', mask, clobber=True)

    print('Finishing mkmask')
    return mask
