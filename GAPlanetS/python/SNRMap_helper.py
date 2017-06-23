"""
contains all functions necessary to use SNRMap.py
"""

from astropy.io import fits
import numpy as np
import math
import os
import sys

def read_file(filename): 
    """
    This function reads a FITS image cube to memory

    Required Inputs:
    1. String containing path to desired pyklipped image file
    
    Example:
    read_file("med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits")

    Last Modified:
    6/19/2017
    """ 
    hdulist = fits.open(filename)
    indivData = hdulist[0].data
    hdulist.close()
    print("Read " + filename  + " in to memory")
    return indivData

def radial_profile(center, y, x):
    """
    This function calculates every pixel's distance from the center, to begin the annuli construction for noise mapping
    
    Required Inputs:
    1. 2-element array of [x][y] coordinates of center
    2. y-value examined
    3. x-value examined

    Last Modified:
    6/19/2017
    """
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.double)
    return r

def generate_radial_profile(filename, saveOutput):
    """          
    This function calls for radial_profile for every pixel in a 2D image.         
    
    Required Inputs: 
    1. if saveOutput is true it won't matter; otherwise filename is a string containing the name of the original input file, for use when saving output
    2. a boolean that is true if you want to save the output to a fits file
                                                                                                                                                                                                       
    Last Modified: 
    6/23/2017 
    
    """ 
    rArray = np.zeros((450,450))
    center = (225.00, 225.00)
    for y in range(450):
        for x in range(450):
            rArray[y][x] = radial_profile(center, y, x)
    if(saveOutput):
        filename = filename[0:-5]
        hdu = fits.PrimaryHDU(rArray)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(filename + "_radprof.fits", overwrite=True)
        print("Wrote " + filename + "_radprof.fits to " + os.getcwd())
    return rArray

def generate_mask(radial_profile, z):
    """
    This function loops through every pixel in a 2D image and determines whether or not to mask its value

    Required Inputs:
    1. Radial profile image
    2. Radius of annuli to mask

    Last Modified:
    6/21/2017
    """
    mask = np.zeros((450,450))
    for y in range(450):
        for x in range(450):
            if (radial_profile[y][x] >= z-.5 and radial_profile[y][x] <= z+.5):
                mask[y][x] = 1
            else:
                mask[y][x] = np.nan
    sys.stdout.write("Generating mask %d of 225   \r" % (z) )
    sys.stdout.flush()
    return mask

def build_mask_cube(mask_cube, filename, saveOutput):
    """
    This function constructs a FITS image cube containing the mask cube generated
    
    Required inputs:
    1. Mask cube data
    2. if saveOutput is true it won't matter; otherwise filename is a string containing the name of the original input file, for use when saving output
    3. a boolean that is true if you want to save the output to a fits file

    Last Modified:
    6/23/2017
    """
    if(saveOutput):
        filename = filename[0:-5]
        hdu = fits.PrimaryHDU(mask_cube)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(filename + "_mask_cube.fits", overwrite=True)
        print("Wrote " +filename+ "_mask_cube.fits to " + os.getcwd())
    return mask_cube

def multiply_by_noise_mask(mask, data):
    """                                                                                                                                                                                                 
    This function is used specifically to mask the planet candidate from the image before noise is calculated (using wedge mask from optional parameter)
    
    Required Inputs:
    1. Mask 2D data
    2. Original data image   

    Last Modified:
    2/28/2017
    """
    newImage = np.zeros((450,450))
    for y in range(450):
        for x in range(450):
            newImage[y][x] = data[y][x] * mask[y][x]
    return newImage

def multiply_by_mask(mask_cube, z, indiv):
    """
    This function multiplies every pixel in a 2D image with the corresponding pixel in a mask image

    Required Inputs:
    1. Mask cube data
    2. Radius to focus annuli
    3. Original data image

    Last Modified:
    6/21/2017
    """    
    data = indiv
    newImage = np.zeros((450,450))
    for y in range(450):
        for x in range(450):
            newImage[y][x] = data[y][x] * mask_cube[z][y][x]
    sys.stdout.write("Generated multiplied image %d of 225   \r" % (z) )
    sys.stdout.flush()
    return newImage

def build_multiplied_cube(multiplied_cube, filename, saveOutput):
    """
    This function constructs a FITS image cube containing the multiplied cube generated                                              

    Required inputs:
    1. Multiplied cube data 
    2. if saveOutput is true it won't matter; otherwise filename is a string containing the name of the original input file, for use when saving output
    3. a boolean that is true if you want to save the output to a fits file
                                         
    Last Modified:  
    6/23/2017 
    """
    if(saveOutput):
        filename = filename[0:-5]
        hdu = fits.PrimaryHDU(multiplied_cube)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(filename + "_multcube.fits", overwrite=True)
        print("Wrote " + filename + "_multcube.fits to " + os.getcwd())
    return multiplied_cube

def calculate_std(image, z):
    """
    This function calculates standard deviation for every pixel in an annulus of given radius z

    Required Inputs:
    1. A 2D image
    2. A value of radius to examine

    Example:
    calculate_std(multiplied_cube[z], z)

    Last Modified:
    6/21/2017
    """
    values = []
    for y in range(450):
        for x in range(450):
            value = image[y][x]
            if not math.isnan(value):
                values.append(value)
    std = np.nanstd(values)
    sys.stdout.write("Calculated standard deviation %d of 225   \r" % (z) )
    sys.stdout.flush()
    
    return std

def replacePixels(mask_cube, z, stds, reference, indiv):
    """
    This function recursively replaces one annulus of pixels at a time, passing the next iteration its current form to maintain continuity in new image

    Required Inputs:
    1. Mask cube data
    2. Radius of annuli to replace
    3. Array of standard deviation data
    4. Reference image --> the output of the previous iteration
    5. The original image, for comparison 

    Example:
    replacePixels(mask_cube, z, stds, noise_cube[z-1], indiv)

    Last Modified:
    6/21/2017
    """
    data = indiv
    new_values = np.zeros((450,450))
    for y in range(450):
        for x in range(450):    
            value = data[y][x] * mask_cube[z][y][x] 
            if not (math.isnan(value)):
                new_values[y][x] = stds[z]
            else:
                new_values[y][x] = reference[y][x]
    sys.stdout.write("Replaced pixels in %d out of 220 slices   \r" % (z) )
    sys.stdout.flush()
    return new_values

def build_noise_map(noise_cube, filename, saveOutput):
    """ 
    This function constructs a FITS image containing the noise map generated    
                                                                                                                                                                                                         
    Required inputs:
    1. Noise cube data 
    2. if saveOutput is true it won't matter; otherwise filename is a string containing the name of the original input file, for use when saving output
    3. a boolean that is true if you want to save the output to a fits file
    
    Last Modified:
    6/23/2017
    """
    noise_map = noise_cube[219]
    if(saveOutput):
        filename = filename[0:-5]
        hdu = fits.PrimaryHDU(noise_map)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(filename + "_noisemap.fits", overwrite=True)
        print("Wrote " + filename + "_noisemap.fits to " + os.getcwd())
    return noise_map

def create_signal_to_noise_map(noise_map, indiv, filename):
    """
    This function divides the signal (original) image by the noise map at every pixel

    Required Inputs:
    1. Noise map data (2D image)
    2. Original data (2D image)
    3. if data not read in from file it won't matter, otherwise is the original filename, to be used when saving output

    Last Modified:
    6/23/2017
    """
    data = indiv
    snr_map = np.zeros((450,450))
    print("Dividing images")
    for y in range(450):
        for x in range(450):
            snr_map[y][x] = data[y][x] / noise_map[y][x]
    if (isinstance(filename, str)):
        filename = filename[0:-5]
        hdu = fits.PrimaryHDU(snr_map)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(filename + '_SNRMap.fits', overwrite=True)
        print("Wrote "+filename+"_SNRMap.fits to " + os.getcwd())
    return snr_map

def implant_custom_mask(theta1, theta2, r1, r2):
    """
    This function creates a custom wedge mask between two radii r1,r2 and in a range of angles theta1,theta2. Intended to mask planet in calculation of standard deviation

    Required Inputs:
    1. Inner angle (0 to 360 degrees, with 0=to the right, counterclockwise)
    2. Outer angle (0 to 360 degrees, with 0=to the right, counterclockwise) MUST BE BIGGER THAN THETA1
    3. Inner radius
    4. Outer radius

    Example:
    implant_custom_mask(100, 150, 60, 70)

    Last Modified:
    2/28/2017
    """
    theta1x = theta1-360
    theta2x = theta2-360
    use_thetaX = False
    if (theta1 > theta2):
        print("Theta1 must be smaller than Theta2")
        return
    if (theta1 < 90 and theta2 > 90):
        use_thetaX = True
    elif (theta1 >= 90 and theta2 > 90):
        use_thetaX = True
    theta1 = math.radians(theta1) + (math.pi/2)
    theta2 = math.radians(theta2) + (math.pi/2)
    if (use_thetaX == True):
        theta1x = math.radians(theta1x) + (math.pi/2)
        theta2x = math.radians(theta2x) + (math.pi/2)
    mask = np.ones((450,450))
    for y in range(450):
        for x in range(450):
            a = x-225
            if (a == 0):
                a = .01
            b = y-225
            if (b == 0):
                b = .01
            r = math.sqrt(a**2 + b**2)
            theta = np.arctan2(b,a)
            if (use_thetaX == False):
                if (r > r1 and r < r2):
                    if (theta > theta1 and theta < theta2):
                        mask[y][x] = np.nan
            elif (use_thetaX == True):
                if (r > r1 and r < r2):
                    if ((theta > theta1 and theta < theta2) or (theta > theta1x and theta < theta2x)):
                        mask[y][x] = np.nan
    return mask

def print2D(array):
    '''
    function for debugging. will print a 450 by 450 array in 2D grid form but at size 45 by 45; useful for checking masks.
    obviously most values are excluded; this is just to be able to quickly check if something is working.
    
    Required Inputs:
    1. a 450 by 450 2D array
    
    Last modified:
    6/19/17
    '''
    for y in range(450):
        if(y%10==0):
            for x in range(450):
                if(x%10==0):
                    if(math.isnan(array[x][y])):
                        print('n', end=' ')
                    else:
                        print((array[x][y]).astype(np.int), end=' ')
            print('')
    return