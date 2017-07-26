"""
Precisely centers an image cube based on finding the center of circular symmetry.

To call: new_circlesym.circsym(filename, maximum radius)
Will output a fits file containing the centered cube, as well as returning the cube as an array

Last modified: 7/24/17
"""
from astropy.io import fits
import numpy as np
import image_registration as ir
import math
import scipy
from astropy.convolution import convolve
from astropy.modeling import models, fitting
from matplotlib import pyplot as plt

def circsym(file,rmax):
    """
    main function for created a centered cube, uses the rest of the functions in this file
    INPUTS:
        file - a string containing the name of the fits file that has been registered and clipped.
        rmax - the maximum radius to look at when determining how 'good' a given point is as the center.
    OUTPUTS:
        final_shifted_cube - an array containing data for the centered cube.
    FILE OUTPUTS:
        (Cont or Line)_clip451_flat_reg_circsym.fits
        NOTE: default is 451 right now, not sure what to do about that
    """
    if(file[0:4]=='Cont'):
        namestring = "Cont_"
    elif(file[0:4]=='Line'):
        namestring = "Line_"
    else:
        namestring = "_"
    
    hdulist = fits.open(file)
    clip_cube = hdulist[0].data
    header_data = fits.getheader(file,0)
    hdulist.close()
    center = clip_cube.shape[1]/2. - 0.5
    print("read in cube")
  
    med = np.nanmedian(clip_cube, axis=0)
    print("constructed median image")
    
    centerX, centerY = sym_center(med, rmax)
 
    for z in range(clip_cube.shape[0]):
        temp = ir.fft_tools.shift2d(clip_cube[z], center-centerX, center-centerY)
        clip_cube[z] = temp
    final_shifted_cube = clip_cube
    print("shifted cube")
    
    hdu = fits.PrimaryHDU(final_shifted_cube, header_data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(str(namestring)+"clip451_flat_reg_circsym.fits", overwrite=True)
    
    return final_shifted_cube
    
def sym_center(image,rmax):
    """
    determines the center of circular symmetry of the median of the cube
    INPUTS:
        image - median image of cube, in 2D array form
        rmax - the maximum radius to look at when determining how 'good' a given point is as the center.
    """ 
    size = image.shape[0] #will be square images
    r = np.arange(0,21)
    grid = np.zeros((r.shape[0],r.shape[0]))
    print("completed setup")
    
    for y in range(r.shape[0]):
        #print("starting row " + str(y))
        for x in range(r.shape[0]):
            yc = r[y] + (size/2.-0.5) - (r.shape[0]/2.-0.5)
            xc = r[x] + (size/2.-0.5) - (r.shape[0]/2.-0.5)
            radii = make_radial_profile(image, xc, yc)
            for k in range(rmax+1):
                sd = np.nanstd(image[(radii>=k) & (radii<k+1)])
                grid[y][x] = grid[y][x] + sd/abs(np.nanmedian(image[(radii>=k) & (radii<k+1)]))
    print("created grid")
    
    minPos = np.argmin(grid)
    pos = np.unravel_index(minPos, grid.shape) #finding pos of min element of grid

    #xcc, ycc = gcntrd(-1 * grid, pos[1], pos[0], 0.5 * len(r)) #calculating centroid for small grid
    xcc, ycc = fit_gaussian(-1*grid, pos[1], pos[0])
    xcc = xcc[0]
    ycc = ycc[0]
    print("grid centered at " + str(xcc)+","+str(ycc))
    
    #converting position from small grid to actual image:
    xc = xcc + (size/2.-0.5) - (r.shape[0]/2.-0.5)
    yc = ycc + (size/2.-0.5) - (r.shape[0]/2.-0.5)
    print('Finishing sym_center; center of image is ' + str(xc)+","+str(yc))
    return (xc, yc)

def make_radial_profile(image, centerX, centerY):
    rArray = np.zeros((image.shape[0],image.shape[1]))
    for y in range(image.shape[0]):
        for x in range(image.shape[1]):
            rArray[y][x] = np.sqrt((x - centerX)**2 + (y - centerY)**2).astype(np.double)
    return rArray

"""
ALTERNATIVE TO GCNTRD THAT I FOUND:
"""
def fit_gaussian(image, xcen, ycen):
    """
    fits a gaussian function to the grid and finds the center of that function
    INPUTS:
        image - 2D array
        xcen, ycen - approximate center
    OUTPUTS:
        xcc, ycc - center of fitted gaussian
    """
    #setup, sets the most negative value (grid is always inverted) as 0
    #so that the grid is positive but still inverted
    min_guess = np.min(image)
    image = image-min_guess
    
    max_guess = image[xcen][ycen]
    g_init = models.Gaussian2D(amplitude=max_guess, x_mean=xcen, y_mean=ycen, x_stddev=5, y_stddev=5)
    fit_g = fitting.LevMarLSQFitter()
    y,x = np.mgrid[:image.shape[0],:image.shape[1]]
    z = image
    g = fit_g(g_init, x, y, z)
    
    xcc = g.x_mean
    ycc = g.y_mean
    
    '''
    #code for making visualizations of the grid, gaussian model, and gaussian fit
    plt.figure(figsize=(10, 3))
    plt.subplot(1, 3, 1)
    plt.imshow(z, origin='upper', interpolation='nearest')
    plt.title("Grid")
    plt.colorbar()
    plt.subplot(1, 3, 2)
    plt.imshow(g(x, y), origin='upper', interpolation='nearest')
    plt.title("Model_fitted")
    plt.colorbar()
    plt.subplot(1, 3, 3)
    plt.imshow(g_init(x, y), origin='upper', interpolation='nearest')
    plt.title("Model_orig")
    plt.colorbar()
    plt.show()
    
    #print(fit_g.fit_info['message'])
    '''
    return (xcc,ycc)
