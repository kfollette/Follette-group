import numpy as np
import math
import statistics as stat 
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.convolution as conv
import os
import sys
import time
from copy import deepcopy

def nameOutput(filename, output):
    
    if (output == None):
        
        #if no output name is specified at runtime and input data is a fits file
        if(isinstance(filename, str)):
            outputName = filename[:-5] + "_SNRMap.fits"
        
        #if data type is not a string (file name), names output file after date and time 
        else:
            outputName = "SNRMap_"  +(time.strftime("%d-%m-%Y")) +"_" +(time.strftime("%H-%M-%S"))+".fits" 
        
    else: 
        outputName = output
        
    return str(outputName)
    

def read_file(filename): 
    """
    This function reads a FITS image cube to memory

    Required Inputs:
    1. String containing path to desired pyklipped image file
    
    Example:
    read_file("med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits")
    
    Written by:
    Elijah Spiro

    Last Modified:
    6/19/2017
    """ 
    hdulist = fits.open(filename)
    indivData = hdulist[0].data
    hdulist.close()
    print("Read " + filename  + " in to memory")
    return indivData



def convertAngle(theta):
    """
    This function takes in an angle as input and converts it to be between 0 and 360 if necessary
    
    Reuired Inputs:
    1. Angle to be converted
            
    Example:
    convertAngle(-30)
    
    Written by:
    Clare Leonard
    
    Last Modified:
    6/28/2016
    
    """
    #modifies angle measurement to fit on a scale from 0 to 360 if it doesn't already
    
    if (theta < 0):
        theta = theta+360
        
    elif (theta >360):
        theta = theta -360
        
    return theta



def inWedge(theta, theta1, theta2):
    
    """
    This function takes in three angles values (in degrees) and returns true if the first of these falls within a wedge starting at theta1 and ending at theta2.
    
    Reuired Inputs:
    1. Position angle of point being tested
    2. Start angle of wedge
    3. End angle of wedge
    
    Examples:
    inWedge(100, 70, 80)
        *would return False
    inWedge(100, 80, 70)    
        *would return True
    
    Written by:
    Clare Leonard
    
    Last Modified:
    6/27/2016
    
    """
    #checks to see if designated angle falls within masked region
    if (theta1 > theta2):
        return (theta <= theta2 or theta >= theta1)
    elif (theta2 > theta1):
        return (theta <= theta2 and theta >= theta1)
    elif (theta2 == theta1):
        return (theta == theta1)
    else: 
        return (False)



def isPlanet(radius, theta, planets):
    
    
    """
    This function takes in the polar coordinates of a point to be tested and a tuple containing lists of parameters for planets
    in the data to be masked.
    
    Reuired Inputs:
    1. Integer radius of point to be tested
    2. Angle coordinate of point to be tested
    3. Tuple containing the following lists:
        a. List of radial coordinates of planets in data
        b. List of corresponding position angles of planets in data (must be same length of a)
        c. List containing radial thickness of desired mask on either side of the planet, followed by the desired angular thickness
            
    Example:
    isPlanet(20, 70, planetData)
        where  >>> planetData = [12, 20, 30, 50], [40, 100, 60, 150], [10, 5]
    
    Written by:
    Clare Leonard
    
    Last Modified:
    6/28/2016
    
    """

    #returns False if there are no planets to mask, showing that the pixel of interest does not fall within any masked region
    if (planets == None):
        return False
  
    #stores lists found in 'planets' tuple as separate variables
    rads, PAs, wid = planets
    
    #stores both arguements of 'wid' parameter in separate variables
    r_wid, pa_wid = wid
    
    for x in range (len(rads)):
       
        #checks to see if point falls within masked radii
        if ((radius < rads[x] + r_wid) and (radius > rads[x] - r_wid)):
            
            #converts position angle and upper and lower angle limits t fall between 0 and 360 degrees
            PA = convertAngle(PAs[x])
            theta1 = PA - pa_wid
            theta2 = convertAngle(theta1)
            theta2 = PA + pa_wid
            theta2 = convertAngle(theta2)
            
            #returns true if the point falls within the bounds of the angle limits, as well as within specified radii
            if(inWedge(theta, theta1, theta2)):
                return True
            
    #returns false if point either doesnt fall between masked radii or masked angles        
    return False
    



def toPolar(x, y):
    
    """
    This function takes a set of pixel coordinates and a set of reference coordinates and transforms the pixel coordinates into polar coordinates.
    
    Reuired Inputs:
    1. Integer x index of pixel
    2. Integer y index of pixel
    3. Integer x index of reference (center) pixel
    4. Integer y index of reference (center) pixel
    
    Exmple:
    toPolar(317, 12, 225, 225)
    
    Written by:
    Clare Leonard
    
    Last Modified:
    6/27/2016
    
    """
    
    #defines pixel radius as the distance from said pixel to the center pixel rounded to an integer
    r = int(np.sqrt((x-XCenter)**2+(y-YCenter)**2))
    #if (r-int(r)>=.5):
        #r = int(r)+1
    #elif (r-int(r)<.5): 
        #r = int(r)
    
    #defines pixel angle 'theta' as the arctangent of the y distance from center divided by the x distance from center
    theta = math.degrees(math.atan2((y-YCenter),(x-XCenter)))
    
    #indexing of the image requires reflecting calculated angle accross the y axis
    theta = theta *-1
    
    #makes sure angle is between 0 and 360
    if(theta<0): 
         theta = theta + 360

    #return calculated polar coordinates
    return (r,theta)



def noisemap(indiv, planets, fwhm, method='stdev'):
    
    """
    This function takes a filename and a list of parameters for objects to mask and outputs a dictionary object of
    integer value radii pointing to the standard deviation for pixel values in the image at that radius from the center.
    
    Reuired Inputs:
    1. Numpy array containing all pixel values for an image
    2. Touple containing the following lists:
        a. List of radial coordinates of planets in data
        b. List of corresponding position angles of planets in data (must be same length of a)
        c. List containing radial thickness of deired mask on either side of the planet, followed by the disired angular thickness
    
    Example:
    stdevMap(indiv, planetData)
        where  >>> planetData = [12, 20, 30, 50], [40, 100, 60, 150], [10, 5]
        and indiv is a numpy array of pixel values
    
    Written by:
    Clare Leonard
    
    Last Modified:
    6/28/2016
    
    """
    #creates empty dictionary objects to store unmasked pixel values and radial standard deviations
    noise_={}
    radialProfs = {}
    radialProfs_abs = {}
    
    #finds the size of the image
    xDim, yDim = np.shape(indiv)

    
    #loops through every pixel in the image
    for x in range (xDim): 
        for y in range (yDim):         
            
            #converts pixel values to polar coordinates
            radius, angle = toPolar(x, y)  
            
            #adds pixel values to radial profile dictionary with the radius as key. ignores masked pixels. 
            if(not isPlanet(radius, angle, planets) and not np.isnan(indiv[x][y])):
                
                #appends pixel value to list associated with radius if the key already exists, adds key and starts new list if not
                if (radius in radialProfs):
                    radialProfs[radius].append(indiv[x][y])
                    radialProfs_abs[radius].append(abs(indiv[x][y]))
                else:
                    radialProfs[radius] = [indiv[x][y],]
                    radialProfs_abs[radius] = [abs(indiv[x][y]), ]
        
                     
    
    #loops through each key in radial profile dictionary, and takes standard deviation of list of pixel values
    #adds standard deviation to stdevs_ dictionary with radius as the key
    #ignores data points if there are too few at a certain radius to take a standard deviation. These pixels will eventually become nans
    for r in radialProfs.keys():
        try:
            if (len(radialProfs[r]) > 8):
                #corrects for small number statistics following Mawet et al. 2008 by penalizing radii where circumference/FWHM is small
                if method == 'stdev':
                    noise_[r]= np.nanstd(radialProfs[r])/np.sqrt(1+(fwhm/(2*math.pi*r)))
                if method == 'med':
                    noise_[r]=np.nanmedian(radialProfs_abs[r])/np.sqrt(1+(fwhm/(2*math.pi*r)))
        except: 
            pass
        
    #returns dictionary holding standard deviations
    return noise_




def create_map(filename, fwhm, smooth = False, planets = None, saveOutput = True, outputName = None, method = 'all', checkmask=False, makenoisemap=False):
    """
    creates signal to noise ratio map of image.
    
    Required Input:
    1. String containing filename of original klipped image OR object containing data already taken from original klipped image

    Optional Inputs:
    1. Tuple containing the following lists:
        a. List of radial coordinates of planets in data
        b. List of corresponding position angles of planets in data (must be same length of a)
        c. List containing radial thickness of desired mask on either side of the planet, followed by the desired angular thickness
         *default value: None*
    2. Boolean designating whether or not to save the completed map to disk 
         *default value: True*
    
    file input example, without mask, saving final map to disk:
        SNRMap.create_map("med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits", saveOutput = True)
    object input example, with mask, without saving final map to disk:
        SNRMap.create_map(data, planets = planetData) 
            (where  >>> planetData = [12, 20, 30, 50], [40, 100, 60, 150], [10, 5])
            
    Written by:
    Clare Leonard

    Last Modified:
    Feb 2019 by KBF - added checkmask and makenoisemap keywords, removed default smooth
    Mar 2019 by KBF - returning max pixel under mask, adding loop over 3rd dimension so can generate 3D SNRmaps, return snrs and masked images
    """

    print('this is the REPAIRED SNRMap code')

    #checks data type of 'filename'
    # if 'filename' is a string, assumes it is a filepath and reads in file
    if(isinstance(filename, str)):
        inp = read_file(filename)
        
    #if data type is not a string, reads in python object holding data
    else:
        inp = np.copy(filename)

    #smooth input image by specified amount
    if smooth > 0:
        print("smoothing")
        gauss = conv.Gaussian2DKernel(stddev=smooth)
        inpsm =conv.convolve(inp, gauss, preserve_nan=True)
        inp = inpsm
    
    #creates dictionary holding the standard deviation of pixlel values at each radius 
    #stdMap = stdevMap(inp, planets, fwhm)
  
    #gets size of pixel value array
    try:
        zdim, ydim, xdim = np.shape(inp)
    except:
        ydim, xdim = np.shape(inp)
        zdim = 1

    global XCenter
    global YCenter 
    XCenter = (xdim-1)/2
    YCenter = (ydim-1)/2

    if method == 'all':
        origmethod='all'
        methods = ['stdev', 'med']
        nmethods=2
    else:
        origmethod=method
        methods = [method]
        nmethods=1

    Output = np.zeros((nmethods, zdim, ydim, xdim))

    if checkmask == True:
        msks = np.ones((zdim,ydim,xdim))
        msk = np.ones((ydim, xdim))

    if makenoisemap == True:
        noises = np.ones((nmethods,zdim,ydim,xdim))
        noise = np.ones((ydim, xdim))

    snrs = np.zeros((nmethods,zdim))
    snr_sums = np.zeros((nmethods,zdim))
    planet_pixels = np.ones((ydim, xdim)) * np.nan
    planet_pixels_pos = np.ones((ydim, xdim)) * np.nan

    #initialize method counter for loop
    methodctr = 0

    for method in methods:

        #KL mode loop
        for s in range (zdim):
            try:
                indiv = np.copy(inp[s,:,:])
            except:
                indiv = np.copy(inp)

        #creates dictionary holding the noise value at each radius for this method
            NoiseMap = noisemap(indiv, planets, fwhm, method=method)


        #loops through all pixels in array
            for x in range (xdim):
                for y in range (ydim):

                #converts indices to polar coordinates
                    radius, angle = toPolar(x,y)

                    if checkmask==True:
                        if (isPlanet(radius, angle, planets)):
                            msk[x][y] = 0

                #if enough pixels have been found to calculate a noise value for this pixels radius,
                # the pixel value is divided by the standard deviation of pixels at that radius
                    try:
                    #if statement prevents a divide by zero warning message
                        if (NoiseMap[radius] == 0):
                            indiv[x][y] = np.nan
             
                        else:
                            indiv[x][y] = indiv[x][y]/NoiseMap[radius]

                        if makenoisemap == True:
                            noise[x][y] = NoiseMap[radius]

                #if no noise value has been calculated, pixel is given a nan value
                    except:
                        indiv[x][y] = np.nan

                # captures an array with just the planet pixels
                # use angular key equivalent to radial one (roughly mimic circular mask)
                    #planets_core = deepcopy(planets)
                    #planets_core[2][1] = np.arctan(planets_core[2][0]/planets[0][0])*180/np.pi
                #print("test radial mask is", planets_core[2][0], "azimuthal mask is now", planets_core[2][1], "instead of", planets[2][1])
                    if (isPlanet(radius, angle, planets)):
                        planet_pixels[x][y]=indiv[x][y]
                    ##includes only positive pixels under mask as planet pixels to avoid including self-subtraction regions in sums
                        if indiv[x][y] > 0:
                            planet_pixels_pos[x][y]=indiv[x][y]

                #store output for this method and # KL modes
                Output[methodctr,s,:,:] = indiv

            if checkmask==True:
                msks[s,:,:]=msk
            if makenoisemap==True:
                noises[methodctr, s,:,:]=noise
            snrs[methodctr,s]=np.nanmax(planet_pixels)
            snr_sums[methodctr,s] = np.nansum(planet_pixels_pos)
            print("max SNR under mask is", snrs[methodctr,s], "for slice", s)
            print("sum of SNRs under mask is", snr_sums[methodctr,s], "for slice", s)
        methodctr += 1

    print('method check', origmethod)
    #saves output to disk if saveOutput designated True
    if (saveOutput == True):
        newname = str(nameOutput(filename, outputName))
        fits.writeto(origmethod + "_"+newname, Output, overwrite=True)
        print("Wrote %s to "%newname + os.getcwd())

        if checkmask==True:
            maskedims = msks*inp
            fits.writeto(origmethod+'_'+newname[:-5]+'_masked.fits', maskedims, overwrite=True)

        if makenoisemap==True:
            fits.writeto(origmethod+'_noisemap.fits', noises, overwrite=True)

    #returns final SNR map
    if checkmask==True:
        return Output, snrs, snr_sums, maskedims
    else:
        return Output, snrs, snr_sums



def getPlanet(snrmap, sep, pa, _range):
    
    #try:
    modeDim, kldim, yDim, xDim = np.shape(snrmap)
    print("check - cube shape is", snrmap.shape)
    #except:
       # yDim, xDim = np.shape(snrmap)
       # zDim = 1
    global XCenter
    global YCenter
    XCenter = (xDim-1)/2
    YCenter = (yDim-1)/2  

    #establish x and y coordinates of planet at given pa and sep
    x = int(sep*math.cos(math.radians(pa+90))+XCenter)
    y = int(sep*math.sin(math.radians(pa+90))+YCenter)

    planet = -100000000

    for mode in range(modedim):
        for kl in range(kldim):
            print("check - mode = ", mode, "kl slice = ", kl)
            for i in range (x-_range, x+_range):
                for j in range (y-_range, y+_range):
                    #finds highest pixel under mask for this planet and returns it
                    if (snrmap[mode][kl][j][i] > planet):
                        planet = snrmap[mode][kl][j][i]
                
    if (planet == -100000000):
        return np.nan
    print("check - dimensions of planet object are", planet.shape)
    return planet


def getNoise(noise_, sli, rad, _range):

    ##pull noise pixels from stdev dictionary
    noisepix = noise_[rad]
    print(len(noisepix), 'noise pixels')

    try:
        zDim, yDim, xDim = np.shape(snrmap)
    except:
        yDim, xDim = np.shape(snrmap)
        zDim = 1
    global XCenter
    global YCenter
    XCenter = (xDim - 1) / 2
    YCenter = (yDim - 1) / 2

    x = int(rad * math.cos(math.radians(pa + 90)) + XCenter)
    y = int(rad * math.sin(math.radians(pa + 90)) + YCenter)

    planet = -100000000

    for i in range(x - _range, x + _range):
        for j in range(y - _range, y + _range):
            if (snrmap[sli][j][i] > planet):
                planet = snrmap[sli][j][i]

    if (planet == -100000000):
        return np.nan
    return planet


def getMedian(snrmap, sli, rad, pa, _range):
    try:
        zDim, yDim, xDim = np.shape(snrmap)
    except:
        yDim, xDim = np.shape(snrmap)
        zDim = 1
    global XCenter
    global YCenter
    XCenter = (xDim - 1) / 2
    YCenter = (yDim - 1) / 2

    x = int(rad * math.cos(math.radians(pa + 90)) + XCenter)
    y = int(rad * math.sin(math.radians(pa + 90)) + YCenter)

    planet = -100000000

    for i in range(x - _range, x + _range):
        for j in range(y - _range, y + _range):
            if (snrmap[sli][j][i] > planet):
                planet = snrmap[sli][j][i]

    if (planet == -100000000):
        return np.nan
    return planet


