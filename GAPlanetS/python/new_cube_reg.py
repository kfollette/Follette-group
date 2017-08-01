"""
Aligns images within a data cube and approximately centers them; 
run new_circlesym as the next step to precisely center the whole cube

To call: new_cube_reg.make_reg_cube(filename, number of reference frame, clip number)
Will output a fits file of the registered and clipped cube, 
  as well as returning the clipped cube as an array to be used if needed.

Last modified: 7/21/17
"""
from astropy.io import fits
import numpy as np
import image_registration as ir
import sys
import math
from scipy.ndimage import gaussian_filter as gauss
'''
Uses image_registration, a package developed to align data cubes based on extended emission.
http://image-registration.readthedocs.io/en/latest/
'''

def make_reg_cube(file, ref, clip):
    """
    main function for created a registered and clipped cube, uses the rest of the functions in this file
    INPUTS: 
        file - a string containing the name of a fits file that has been preprocessed.
        ref - number of the reference frame that the images will be registered to.
        clip - the size of the clipped cube that will be output.
    OUTPUTS:   
        clip_cube - an array containing data for the clipped and registered cube
    FILE OUTPUTS:
        (Cont or Line)_clip(clip)_flat_reg.fits
    EXAMPLE:
        make_reg_cube("Line_flat_preproc.fits", 818, 451, line)
        outputs a file named "Line_clip451_flat_reg.fits"
        make_reg_cube("Line_flat_preproc.fits", 818, 451)
        outputs a file named "_clip451_flat_reg.fits"
    """
    global original 
    global reg_cube
    
    if(file[0:4]=='Cont'):
        namestring = "Cont_"
    elif(file[0:4]=='Line'):
        namestring = "Line_"
    else:
        namestring = "_"
    
    #open file
    hdulist = fits.open(file)
    original = hdulist[0].data
    header_data = fits.getheader(file,0)
    hdulist.close()
    print("read in cube")
    
    reg_cube = register(original, ref)
    print("registered cube")
    
    center = find_max(reg_cube[ref])
    print(center)
    
    clip_cube = clip_all(reg_cube, center[0], center[1], clip)
    
    print("writing to fits file")
    hdu = fits.PrimaryHDU(clip_cube, header_data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(str(namestring)+"clip"+str(clip)+"_flat_reg.fits", overwrite=True)
    
    return clip_cube

def register(cube, ref):
    """
    runs through each layer of a data cube and determines how far to shift it
    to match a given reference layer, and then shifts it.
    INPUTS:
        cube - data cube in 3 dimensional array form.
        ref - number of the reference layer.
    OUTPUTS:
        cube - 3 dimensional array in which each layer is shifted to match 
             the reference. Will have same dimensions as the original.
    """
    size = cube.shape[0]
    for z in range(size):
        #determines how far to shift to match ref image
        shift = ir.chi2_shift(cube[ref],cube[z])
        shiftX, shiftY = shift[0], shift[1]
        
        #actually shifts
        shifted_image = ir.fft_tools.shift2d(cube[z],-shiftX,-shiftY)
        cube[z] = shifted_image
    return cube

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

def find_max(image):
    """
    finds the brightest pixel in the image within a range around the center
    INPUTS:
        image in 2 dimensional array form
    OUTPUTS:
        coordinates of brightest pixel within +/- 100 from the image center
    """
    m = 0
    centerx = 0
    centery = 0
    for y in range(int(image.shape[0]/2)-100,int(image.shape[0]/2)+100):
        for x in range(int(image.shape[1]/2)-100,int(image.shape[1]/2)+100):
            if (image[y][x] > m):
                m = image[y][x]
                centerx = x
                centery = y
    return centerx, centery
