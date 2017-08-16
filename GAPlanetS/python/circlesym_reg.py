'''
Registers an image cube by finding the center of circular symmetry in each layer 
of the cube (to sub-pixel precision) and shifting the layer so that is the center of the image.
Also crops the cube to a given size.

To call: circlesym_reg.register(filename,rmax,clip size)

last modified: 8/16/17
'''
import image_registration as ir
from astropy.io import fits
import numpy as np
#from astropy.convolution import convolve
#from astropy.modeling import models, fitting
#from matplotlib import pyplot as plt
import math
import astropy.convolution as conv
from scipy import optimize

def register(file,rmax,clip):
    '''
    reads in an image cube from a fits file, finds and shifts to the center of circular symmetry for each layer,
    clips the cube to a given size, and writes the result to a fits file. 
    This is the main function; it calls the rest of the functions in this module.
    
    INPUTS:
        file - the input filename. Ideally is either Cont_flat_preproc.fits or Line_flat_preproc.fits,
        if you're following the image processing pipeline.
        rmax - the maximum radius you want it to check from a given center pixel for standard deviations 
        when determining if that pixel is a good center or not.
        clip - desired size of the resulting image cube.
    OUTPUTS:
        clip_cube - a 3D array containing the resulting cube.
    FITS FILE OUTPUTS:
        (Line or Cont)_clip(clip size)_flat_reg_circsym.fits
    '''
    if(file[0:4]=='Cont'):
        namestring = "Cont_"
    elif(file[0:4]=='Line'):
        namestring = "Line_"
    else:
        namestring = "_"
    
    global centerX
    global centerY
    
    hdulist = fits.open(file)
    original = hdulist[0].data
    #original = original[212:215]
    header_data = fits.getheader(file,0)
    hdulist.close()
    centerX = int(original.shape[2]/2. - 0.5)
    centerY = int(original.shape[1]/2. - 0.5)
    print("read in cube")
    
    size = original.shape[0]
    
    for z in range(size):
        print(z)
        xc, yc = sym_center(original[z],rmax)
        temp = ir.fft_tools.shift2d(original[z], centerX-xc, centerY-yc)
        original[z] = temp
    shifted_cube = original
    print("shifted cube")
    
    clip_cube = clip_all(shifted_cube, centerX, centerY, clip)
    
    print("writing to fits file")
    hdu = fits.PrimaryHDU(clip_cube, header_data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(str(namestring)+"clip"+str(clip)+"_flat_reg_circsym.fits", overwrite=True)
    
    return clip_cube
    
def clip_all(cube, centerX, centerY, clip):
    """
    clips data cube to desired square size.
    INPUTS: 
        cube - data cube, in 3 dimensional array form.
        centerX, centerY - coordinates in original cube that will be the center of the new cube.
        clip - size to clip to.
    OUTPUTS:
        clipped_cube - clipped data cube in 3 dimensional array form.
    """
    clipped_cube = np.zeros((cube.shape[0],clip,clip))
    half = math.floor(clip/2)
    for z in range(cube.shape[0]):
        for y in range(-half,half+1):
            for x in range(-half,half+1):
                '''
                check if you're trying ot pull values from outside vertical bounds of image.
                if so, that area will be 0, separated by a line of nans.
                doesn't check x since the original cubes are much wider in x.
                (will not be an issue in this module, but was useful when used in previous modules)
                '''
                if centerY+y==cube.shape[1]:
                    clipped_cube[z][y+half][x+half] = np.nan
                elif centerY+y>cube.shape[1]:
                    clipped_cube[z][y+half][x+half] = 0
                else:
                    clipped_cube[z][y+half][x+half] = cube[z][centerY+y][centerX+x]
    return clipped_cube
        
def sym_center(image,rmax):
    '''
    finds the center of circular symmetry of a single image.
    INPUTS: 
        image - image in the form of a 2D array
        rmax - the maximum radius you want it to check from a given center pixel for standard deviations 
        when determining if that pixel is a good center or not.
    OUTPUTS:
        xc,yc - center of circular symmetry of the image
    '''
    global xMax
    global yMax
    
    #smoothing before finding approximate center
    #helps avoid issues with cosmic rays
    gauss = conv.Gaussian2DKernel(stddev=6)
    new_im = conv.convolve(image, gauss)

    #find approximate center of image
    xMax = 0
    yMax = 0
    max_val = 0
    for y in range(centerY-200,centerY+200):
        for x in range(centerX-200,centerX+200):
            if new_im[y][x] > max_val:
                xMax = x
                yMax = y
                max_val = new_im[y][x]
    
    #create 21 by 21 grid of values showing whether each possible pixel in the grid is a good center
    #lower values mean better center            
    grid = np.zeros((21,21))          
    for y in range(21):
        print("starting row " + str(y))
        for x in range(21):
            yc = yMax - (10) + y
            xc = xMax - (10) + x
            radii = make_radial_profile(image, xc, yc)
            for k in range(rmax+1):
                ids = np.where((radii>=k) & (radii<k+1))
                for a in range(ids[0].shape[0]):
                    ids[0][a] = ids[0][a] + yMax-50
                for b in range(ids[1].shape[0]):
                    ids[1][b] = ids[1][b] + xMax-50
                sd = np.nanstd(image[ids])
                grid[y][x] = grid[y][x] + sd/abs(np.nanmedian(image[ids]))
    print("created grid")

    #fitting a 2D gaussian to the grid to determine actual 'best' center
    maxval = np.nanmax(grid)
    data = grid*-1
    data = data+maxval
    p = fitgaussian(data)
    xcc = p[2]
    ycc = p[1]
    print("grid centered at " + str(xcc)+","+str(ycc))
    
    #converting position from small grid to actual image:
    xc = xMax - (21/2.-0.5) + xcc
    yc = yMax - (21/2.-0.5) + ycc
    print('Finishing sym_center; center of image is ' + str(xc)+","+str(yc))
    return (xc, yc)
    
def make_radial_profile(image, xc, yc):
    '''
    creates a 101 by 101 grid containing each pixel's radial distance from a given image center
    INPUTS:
        image - image in form of a 2D array
        xc,yc - chosen center coordinates (in regard to image, not radial grid)
    OUTPUTS:
        rArray - radial profile grid in form of a 2D array
    '''
    #xc and yc are coordinates in original image
    rArray = np.zeros((101,101))
    for y in range(101):
        for x in range(101):
            xi = x + xMax-50 #x coordinate in original image
            yi = y + yMax-50 #y coordinate in original image
            rArray[y][x] = np.sqrt((xi - xc)**2 + (yi - yc)**2).astype(np.double)
    return rArray

#######################################
####### FUNCTIONS FOR GAUSSIAN ########
#######################################

'''all of these taken from here: https://stackoverflow.com/a/30791245'''

def gaussian(p, x, y):
    height, center_x, center_y, width_x, width_y = p
    return height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    total = np.nansum(data)
    X, Y = np.indices(data.shape)
    center_x = np.nansum(X*data)/total
    center_y = np.nansum(Y*data)/total
    row = data[int(center_x), :]
    col = data[:, int(center_y)]
    width_x = np.nansum(np.sqrt(abs((np.arange(col.size)-center_y)**2*col))
                        /np.nansum(col))
    width_y = np.nansum(np.sqrt(abs((np.arange(row.size)-center_x)**2*row))
                        /np.nansum(row))
    height = np.nanmax(data)
    return height, center_x, center_y, width_x, width_y

def errorfunction(p, x, y, data):
    return gaussian(p, x, y) - data

def fitgaussian(data):
    params = moments(data)
    X, Y = np.indices(data.shape)
    mask = ~np.isnan(data)
    x = X[mask]
    y = Y[mask]
    data = data[mask]
    p, success = optimize.leastsq(errorfunction, params, args=(x, y, data))
    return p

