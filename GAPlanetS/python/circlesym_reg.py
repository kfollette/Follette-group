import image_registration as ir
from astropy.io import fits
import numpy as np
from astropy.convolution import convolve
from astropy.modeling import models, fitting
from matplotlib import pyplot as plt
import math

def register(file,rmax,clip):
    
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
    #original = original[0:10]
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
    
    print("writing to fits file")
    hdu = fits.PrimaryHDU(shifted_cube, header_data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("temp_cube.fits", overwrite=True)
    
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
                '''
                if centerY+y==cube.shape[1]:
                    clipped_cube[z][y+half][x+half] = np.nan
                elif centerY+y>cube.shape[1]:
                    clipped_cube[z][y+half][x+half] = 0
                else:
                    clipped_cube[z][y+half][x+half] = cube[z][centerY+y][centerX+x]
    return clipped_cube
        
def sym_center(image,rmax):
    global xMax
    global yMax
    xMax = 0
    yMax = 0
    max_val = 0
    for y in range(centerY-200,centerY+200):
        for x in range(centerX-200,centerX+200):
            if image[y][x] > max_val:
                xMax = x
                yMax = y
                max_val = image[y][x]
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
    
    minPos = np.argmin(grid)
    pos = np.unravel_index(minPos, grid.shape) #finding pos of min element of grid
    
    xcc, ycc = fit_gaussian(-1*grid, pos[1], pos[0])
    xcc = xcc[0]
    ycc = ycc[0]
    print("grid centered at " + str(xcc)+","+str(ycc))
    
    #converting position from small grid to actual image:
    xc = xMax - (21/2.-0.5) + xcc
    yc = yMax - (21/2.-0.5) + ycc
    print('Finishing sym_center; center of image is ' + str(xc)+","+str(yc))
    return (xc, yc)
    
def make_radial_profile(image, xc, yc):
    #xc and yc are coordinates in image
    rArray = np.zeros((101,101))
    for y in range(101):
        for x in range(101):
            xi = x + xMax-50 #x coordinate in image
            yi = y + yMax-50 #y coordinate in image
            rArray[y][x] = np.sqrt((xi - xc)**2 + (yi - yc)**2).astype(np.double)
    return rArray

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
