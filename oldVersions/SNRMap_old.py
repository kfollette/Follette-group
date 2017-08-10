"""
This code is to be used as a metric for evaluating the goodness of a klipped image. 
It takes in a "final product" klipped image output by pyklip and produces a signal:noise ratio image. 
The higher the value of signal vs. noise, the better we can trust the image. 

To call this from another python file: (to be added)

"""

import SNRMap_helper_old as snr
import numpy as np

def create_map(filename, saveOutput, mask=np.ones((450,450))):
    """
    creates the actual signal-to-noise map using functions stored in SNRMap_helper.
    
    Required Inputs:
    1. String containing filename of original klipped image OR object containing data already taken from original klipped image
    2. Boolean for whether you want to save the function outputs to fits files rather than just keeping them as objects.
        DO NOT set saveOutput to True if your input is an object and not a file; you will break it if you do.
    
    Optional Inputs:
    1. mask = custom mask (2D 'image' array)
    If you don't have a wedge mask set up but you need one, first call: 
        SNRMap_helper.implant_custom_mask(theta_start, theta_end, inner_radius, outer_radius)
        the angles measure counterclockwise starting with 0 directly to the right, and are in degrees. 
        the starting angle must be smaller than the ending angle.
        
    file input example:
        SNRMap.create_map("med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits", True)
    object input example, with mask:
        SNRMap.create_map(data, False, mask=custom_mask)

    Last Modified:
    6/23/2017
    """
    if(isinstance(filename, str)):
        indiv = snr.read_file(filename)
    else:
        indiv = filename
    wedge_masked_original = snr.multiply_by_noise_mask(mask, indiv)
    radial_profile = snr.generate_radial_profile(filename, saveOutput)
    mask_cube = np.zeros((225,450,450))
    for z in range(1,225):
        mask_cube[z] = snr.generate_mask(radial_profile, z)
    mask_cube = snr.build_mask_cube(mask_cube, filename, saveOutput)
    multiplied_cube = np.zeros((225,450,450))
    for z in range(1,225):
        multiplied_cube[z] = snr.multiply_by_mask(mask_cube, z, wedge_masked_original)
    multiplied_cube = snr.build_multiplied_cube(multiplied_cube, filename, saveOutput)
    stds = []
    for z in range(5, 225):
        stds.append(snr.calculate_std(multiplied_cube[z], z))
    noise_cube = np.zeros((220,450,450))
    for z in range(220):
        if (z == 0):
            noise_cube[z] = snr.replacePixels(mask_cube, z, stds, indiv, indiv)
        else:
            noise_cube[z] = snr.replacePixels(mask_cube, z, stds, noise_cube[z-1], indiv)
    noise_map = snr.build_noise_map(noise_cube, filename, saveOutput)
    snr_map = snr.create_signal_to_noise_map(noise_map, indiv, filename)
    return snr_map, noise_map, multiplied_cube, mask_cube, radial_profile
