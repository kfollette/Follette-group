"""
Written by Elijah Spiro [espiro18@amherst.edu]

This code is to be used as a metric for evaluating the goodness of a klipped image. It takes in a "final product" klipped image output by pyklip and produces a signal:noise ratio image. The higher the value of signal vs. noise, the better we can trust the image. 

This code is meant to be made generic enough that it can be called repeatedly and its results categorized as one step in the chain of automating the process of determining ideal pyklip parameters for a given data set.

Eventual program flow:
1. Master code runs pyklip on a given data set with intial best guess parameters (movement, IWA, annuli, etc). 
2. THIS CODE IS CALLED AND RUN ON THE OUTPUT OF PYKLIP, AND THE SNR MAP PRODUCED IS ANALYZED
3. Master code runs pyklip again, varying parameters slightly in an MCMC algorithm to avoid local maxima and best explore parameter space
4. THIS CODE IS CALLED AGAIN ON THE NEW KLIPPED IMAGE, OUTPUT IS AGAIN ANALYZED, AND WE DETERMINE IF THE CHANGES IMPROVED OR WORSENED THE QUALITY OF THE RESULT
5. This process is repeated indefinitely until ideal parameters are found for the given data set, and the best possible quality pyklipped image has been produced 
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import math
import os


def read_file(filename):
    """
    This function reads a FITS image cube to memory

    Required Inputs:
    1. String containing path to desired pyklipped image file
    
    Example:
    read_file("med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits")

    Last Modified:
    2/26/2017
    """ 
    os.chdir("../HD142527/8Apr14/revamped/")
    hdulist = fits.open(filename)
    indivData = hdulist[0].data
    hdulist.close()
    print("Read " + filename  + " in to memory")
    indiv = indivData
    return indivData

def radial_profile(data, center, y, x):
    """
    This function calculates every pixel's distance from the center, to begin the annuli construction for noise mapping
    
    Required Inputs:
    1. 2D data cube
    2. 2-element array of [x][y] coordinates of center
    3. y-value examined
    4. x-value examined
    
    Example:
    radial_profile(indiv, center, y, x)

    Last Modified:
    2/26/2017
    """
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.double)
    return r

def generate_radial_profile(indiv):
    """                                                                                                                                                                      
    This function calls for radial_profile for every pixel in a 2D image. Once the pixel distances are all determined, it constructs a new FITS file           
    
    Required Inputs:                                                                                                                                                                                     
    1. 2D data cube                                                                                                                                                                                          
    Example:                                                                                                                                                                                          
    radial_profile(data)                                                                                                                                                                     
                                                                                                                                                                                                       
    Last Modified:                                                                                                                                                                                           2/26/2017                                                                                                                                                                                      
    """ 
    rArray = np.zeros((450,450))
    center, radi = (225.00, 225.00), 50
    for y in range(450):
        for x in range(450):
            rArray[y][x] = radial_profile(indiv, center, y, x)
    hdu = fits.PrimaryHDU(rArray)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("radial_profile.fits", overwrite=True)
    print("Wrote radial_profile.fits to " + os.getcwd())

def read_radial_profile():
    """
    This function reads in a FITS image to memory, skipping the step of generating a radial profile if the FITS file is already in place
    
    Required Inputs:
    None

    Example:
    radial_profile = read_radial_profile()

    Last Modified:
    2/26/2017
    """
    hdulist = fits.open("radial_profile.fits")
    radial_profile = hdulist[0].data
    hdulist.close()
    print("Read radial_profile.fits back in to memory")
    return radial_profile



def generate_mask(radial_profile, z):
    """
    This function loops through every pixel in a 2D image and determines whether or not to mask its value

    Required Inputs:
    1. Radial profile image
    2. Radius of annuli to mask

    Example:
    generate_mask(114)

    Last Modified:
    2/28/2017
    """
    mask = np.zeros((450,450))
    for y in range(450):
        for x in range(450):
            if (radial_profile[y][x] >= z-.5 and radial_profile[y][x] <= z+.5):
                mask[y][x] = 1
            else:
                mask[y][x] = np.nan
    print("Generating mask " + str(z) + " of 225")
    return mask

def build_mask_cube(mask_cube):
    """
    This function constructs a FITS image cube containing the mask cube generated
    
    Required inputs:
    1. Mask cube data

    Example:
    build_mask_cube(mask_cube_data)

    Last Modified:
    2/26/2017
    """
    hdu = fits.PrimaryHDU(mask_cube)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("mask_cube.fits", overwrite=True)
    print("Wrote mask_cube.fits to " + os.getcwd())

def read_mask_cube():
    """                                                                                                                                                                                                      This function reads in a FITS image to memory, skipping the step of generating a new mask if the FITS file is already in place                                                                     
                                                                                                                                                                                                         
    Required Inputs:                                                                                                                                                                                   
    None                                                                                                                                                                                                    
    Example:                                                                                                                                                                                             
    mask_cube = read_mask_cube()                                                                                                                                                                  
                                                                                                                                                                                                         
    Last Modified:                                                                                                                                                                                      
    2/2/2017                                                                                                                                                                                                 """   
    hdulist = fits.open("mask_cube.fits")
    mask_cube = hdulist[0].data
    hdulist.close()
    print("Read mask_cube.fits back in to memory")
    return mask_cube

def multiply_by_noise_mask(mask, data):
    """                                                                                                                                                                                                 
    This function is used specifically to mask the planet candidate from the image before noise is calculated (using wedge mask from optional parameter)
    
    Required Inputs:
    1. Mask 2D data
    2. Original data image
    
    Example:                                                                                                                                                                                             
    multiply_by_mask(mask_data, indiv)                                                                                                                                                                                                                                                                                                                                                                                Last Modified:                                                                                                                                                                                     
    2/28/2017                                                                                                                                                                                                """
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

    Example:
    multiply_by_mask(mask_data, 114)

    Last Modified:
    2/26/2017
    """    
    data = indiv
    newImage = np.zeros((450,450))
    for y in range(450):
        for x in range(450):
            newImage[y][x] = data[y][x] * mask_cube[z][y][x]
    print("Generated multiplied image " + str(z) + " of 225")
    return newImage

def build_multiplied_cube(multiplied_cube):
    """                                                                                                                                                                                                      This function constructs a FITS image cube containing the multiplied cube generated                                              

    Required inputs:
    1. Multiplied cube data                                                                                                                                                                                                                                                                                                           
    Example:                                                                                                                                                                                            
    build_multiplied_cube(multiplied_cube_data)                                                                                                                                                                          
    Last Modified:                                                                                                                                                                                           2/26/2017                                                                                                                                                                                                """
    hdu = fits.PrimaryHDU(multiplied_cube)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("multiplied_cube.fits", overwrite=True)
    print("Wrote multiplied_cube.fits to " + os.getcwd())

def read_multiplied_cube():
    """                                                                                                                                                                                                      This function reads in a FITS image to memory, skipping the step of generating a new multiplied cube if the FITS file is already in place                                              
    
    Required Inputs:                                                                                                                                                                                         None                                                                                                                                                                                               
                                                                                                                                                                                                         
    Example:                                                                                                                                                                                            
    multiplied_cube = read_multiplied_cube()                                                                                                                                                                                                                                                                                                                                                                          Last Modified:                                                                                                                                                                                           2/2/2017                                                                                                                                                                                                 """
    hdulist = fits.open("multiplied_cube.fits")
    multiplied_cube = hdulist[0].data
    hdulist.close()
    print("Read multiplied_cube.fits back in to memory")
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
    2/26/2017
    """
    values = []
    for y in range(450):
        for x in range(450):
            value = image[y][x]
            if not math.isnan(value):
                values.append(value)
    std = np.nanstd(values)
    print("Calculated standard deviation " + str(z) + " of 225") 
    return std

def graph_stds(stds):
    """
    This is an optional function to display the graph of standard deviation as a function of radius during calculations
    
    Required Inputs:
    1. An array of standard deviation values

    Example:
    graph_stds(std_data)

    Last Modified:
    2/26/2017
    """
    plt.title("Standard Deviation as a function of Radius")
    plt.plot(stds, "b--")
    plt.ylabel("Standard Deviation")
    plt.xlabel("Radius")
    plt.show()

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
    2/26/2017
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
    print("Replaced pixels in " + str(z+1) + " out of 220 slices")
    return new_values

def build_noise_cube(noise_cube):
    """                                                                                                                                                                                                      This function constructs a FITS image cube containing the noise cube generated                                                                                                                      
                                                                                                                                                                                                         
    Required inputs:                                                                                                                                                                                         1. Noise cube data                                                                                                                                                                                                                                                                                                                                                                               
    Example:                                                                                                                                                                                                 build_noise_cube(noise_cube_data)                                                                                                                                                             \
    
    Last Modified:                                                                                                                                                                                           2/26/2017                                                                                                                                                                                                """
    noise_map = noise_cube[219]
    hdu = fits.PrimaryHDU(noise_map)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("noise_map.fits", overwrite=True)
    print("Wrote noise_map.fits to " + os.getcwd())

def read_noise_map():
    """                                                                                                                                                                                                      This function reads in a FITS image to memory, skipping the step of generating a noise profile if the FITS file is already in place                                                                     
    Required Inputs:                                                                                                                                                                                         None                                                                                                                                                                                                     
    Example:                                                                                                                                                                                                 noise_map = read_noise_map()                                                                                                                                                                   
                                                                                                                                                                                                         
    Last Modified:                                                                                                                                                                                           2/26/2017                                                                                                                                                                                                """
    hdulist = fits.open("noise_map.fits")
    noise_map = hdulist[0].data
    hdulist.close()
    print("Read noise_map.fits back in to memory")
    return noise_map

def create_signal_to_noise_map(noise_map, indiv):
    """
    This function divides the signal (original) image by the noise map at every pixel

    Required Inputs:
    1. Noise map data (2D image)
    2. Original data (2D image)

    Example:
    create_signal_to_noise_map(noise_data, original_data)

    Last Modified:
    2/26/2017
    """
    data = indiv
    snr_map = np.zeros((450,450))
    print("Dividing images")
    for y in range(450):
        for x in range(450):
            snr_map[y][x] = data[y][x] / noise_map[y][x]
    hdu = fits.PrimaryHDU(snr_map)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("snr_map.fits", overwrite=True)
    print("Wrote snr_map.fits to " + os.getcwd())
    return snr_map

def implant_custom_mask(theta1, theta2, r1, r2):
    """
    This function creates a custom wedge mask between two radii r1,r2 and in a range of angles theta1,theta2. Intended to mask planet in calculation of standard deviation

    Required Inputs:
    1. Inner angle (0-360 degrees, with 0=straight up)
    2. Outer angle (0-360 degrees, with 0=straight up) MUST BE BIGGER THAN THETA1
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
            #print("a = " + str(a) + " , b = " + str(b) + " , x = " + str(x) + " , y = " + str(y) + " , r = " + str(r) + " , theta = " + str(theta))
            if (use_thetaX == False):
                if (r > r1 and r < r2):
                    if (theta > theta1 and theta < theta2):
                        mask[y][x] = np.nan
            elif (use_thetaX == True):
                if (r > r1 and r < r2):
                    if ((theta > theta1 and theta < theta2) or (theta > theta1x and theta < theta2x)):
                        mask[y][x] = np.nan
    hdu = fits.PrimaryHDU(mask)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("custom_mask.fits", overwrite=True)
    print("Wrote custom_mask.fits to " + os.getcwd())
    return mask
    

def SNRMap(filename, from_scratch, graph, mask):
    """
    This function ties together all of the other functions in this file. Call it once with proper inputs and the rest of the program will run in the correct order

    Required Inputs:
    1. String containing filename of original klipped image
    2. Boolean (answer true/false) for whether it's the first time running on a given file (some steps may be skipped on multiple runs)

    Optional Inputs:
    1. graph=Boolean (answer true/false) for whether to display the standard deviation vs. radius graph along the way (program will freeze until the graph is closed)
    2. mask=Call to function to generate custom mask, if desired

    Example:
    SNRMap("med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits", False, graph=False, mask=implant_custom_mask(100,120,30,35))  

    Last Modified:
    2/26/2017
    """
    indiv = read_file(filename)
    wedge_masked_original = multiply_by_noise_mask(mask, indiv)
    generate_radial_profile(indiv)
    radial_profile = read_radial_profile()
    if (from_scratch == True):
        mask_cube = np.zeros((225,450,450))
        for z in range(1,225):
            mask_cube[z] = generate_mask(radial_profile, z)
        build_mask_cube(mask_cube)
    mask_cube = read_mask_cube()
    if (from_scratch == True):
        multiplied_cube = np.zeros((225,450,450))
        for z in range(1,225):
            multiplied_cube[z] = multiply_by_mask(mask_cube, z, wedge_masked_original)
        build_multiplied_cube(multiplied_cube)
    multiplied_cube = read_multiplied_cube()
    stds = []
    for z in range(5, 225):
        stds.append(calculate_std(multiplied_cube[z], z))
    if (graph == True):
        graph_stds(stds)
    noise_cube = np.zeros((220,450,450))
    for z in range(220):
        if (z == 0):
           noise_cube[z] = replacePixels(mask_cube, z, stds, indiv, indiv)
        elif not (z == 0):
            noise_cube[z] = replacePixels(mask_cube, z, stds, noise_cube[z-1], indiv)
    if (from_scratch == True):
        build_noise_cube(noise_cube)
    noise_map = read_noise_map()
    snr_map = create_signal_to_noise_map(noise_map, indiv)





#########################################################################################################################################################################################################
###########################################################################       RUN PROGRAM        ####################################################################################################

SNRMap("med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits", True, False, implant_custom_mask(100,160,15,30))
#########################################################################################################################################################################################################
