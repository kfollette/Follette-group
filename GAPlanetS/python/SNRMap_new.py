import numpy as np
import math
import statistics as stat 
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.convolution as conv
import pyklip.klip as klip
import os
import sys
import time
from copy import deepcopy
    

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
    Data = fits.getdata(filename)
    Head = fits.getheader(filename)
    return Data, Head



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
    
    Output:
    isplanetpix: Boolean indicating whether there is a planet at the given r, theta location
    whichplanet: index corresponding to which planet in the list it belongs to

    Written by:
    Clare Leonard 2016
    
    Last Modified:
    2/25/20 by KBF - output now tuple returning boolean and number of planet. Simplified so isplanetpix is just false by default
                     instead of having multiple test conditions
    2/26/20 by KBF - fixed bug in planets spanning PA=0
    3/2/20 by KBF - adding core pix for averaging
    
    """

    #initialize vairable defaults
    isplanetpix = False
    iscorepix = False
    whichplanet = np.nan

    #stores lists found in 'planets' tuple as separate variables
    rads, PAs, wid = planets

    #stores both arguements of 'wid' parameter in separate variables
    r_wid, pa_wid = wid

    for x in range (len(rads)):
           
        #checks to see if point falls within masked radii
        if ((radius < rads[x] + r_wid) and (radius > rads[x] - r_wid)):

            #converts position angle and upper and lower angle limits t fall between 0 and 360 degrees
            PA = convertAngle(PAs[x])
            theta1 = convertAngle(PA - pa_wid)
            theta2 = convertAngle(PA + pa_wid)
            core_wid = np.arctan(r_wid/radius)*180/np.pi
            core1 = convertAngle(PA - core_wid)
            core2 = convertAngle(PA + core_wid)

            #returns true if the point falls within the bounds of the angle limits, as well as within specified radii
            if(inWedge(theta, theta1, theta2)):
                isplanetpix = True
                if(inWedge(theta, core1, core2)):
                    iscorepix = True
                #records which planet in the list this pixel corresponds to
                whichplanet = x
                  
    return(isplanetpix, iscorepix, whichplanet)


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
            
            #if planets are present, mask them
            if planets != False:
                isplanetpix, iscorepix, whichplanet = isPlanet(radius, angle, planets)
            else:
                isplanetpix = False

            #adds pixel values to radial profile dictionary with the radius as key. ignores masked pixels. 
            if(not isplanetpix and not np.isnan(indiv[x][y])):
                
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




def create_map(filename, fwhm, head = None, smooth = False, planets=False, saveOutput = True, outputName=False, ctrlrad=30, method = 'all', checkmask=False, makenoisemap=False, sigma = 5, verbose = False):
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
    Sept 2019 by KBF - misc. cleanup, added additional SNR methodology to maps (median = noise), which are now 4D, added sums of snrs under mask as return
    """
    #print('this is the REPAIRED SNRMap code')

    #checks data type of 'filename'
    # if 'filename' is a string, assumes it is a filepath and reads in file
    if filename[-5:] == '.fits':
        if verbose is True:
            print("found fits file", filename)
        inp, head = read_file(filename)
        outname = filename[:-5]+'_'+method+'snrmap'
        
    #if data type is not a string, reads in python object holding data
    else:
        inp = np.copy(filename)
        if head == None:
            head = fits.Header()
            head['NKLMODES'] = str(inp.shape[0])
            if outputName==False:
                print("Must specify outputName if passing in an existing array")
                return
            else:
                outname=outputName[:-5]

    #smooth input image by specified amount
    if smooth > 0:
        nkl = inp.shape[0]
        for kl in np.arange(nkl):
            inp[kl,:,:] = klip.nan_gaussian_filter(inp[kl,:,:], smooth)

  
    #gets size of pixel value array
    try:
        kldim, ydim, xdim = np.shape(inp)
    except:
        ydim, xdim = np.shape(inp)
        kldim = 1

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

    Output = np.zeros((nmethods, kldim, ydim, xdim))

    if checkmask == True:
        msks = np.ones((kldim,ydim,xdim))
        msk = np.ones((ydim, xdim))

    if makenoisemap == True:
        noises = np.ones((nmethods,kldim,ydim,xdim))
        noise = np.ones((ydim, xdim))

    #pull values from planets needed for spurious pixels at same radius counts
    if planets != False:
        rads, PAs, wid = planets
        nplanets = len(rads)

        #make some blank arrays for populating with planet statistics
        snrs = np.zeros((nmethods,kldim,nplanets))
        snr_sums = np.zeros((nmethods,kldim,nplanets))
        snr_spurious = np.ones((nmethods, kldim, 5)) * np.nan
        planet_pixels = np.ones((ydim, xdim, nmethods, kldim, nplanets)) * np.nan
        planet_pixels_pos = np.ones((ydim, xdim, nmethods, kldim, nplanets)) * np.nan
        npospix=np.zeros((nmethods,kldim,nplanets))

    #initialize method counter for two SNR computation methods
    methodctr = 0
    for method in methods:

        #KL mode loop
        for s in range (kldim):
            try:
                indiv = inp[s,:,:]
            except:
                indiv = inp

            #creates dictionary holding the noise value at each radius for this method and kl mode
            NoiseMap = noisemap(indiv, planets, fwhm, method=method)
            fivesig = 0
            fivesig_atmask=0
            fivesig_inmask = 0
            allplanetpix = 0
            notplanetpix = 0
            
            # loops through all pixels in array
            for x in range (xdim):
                for y in range (ydim):
                    #converts x/y indices to polar coordinates
                    radius, angle = toPolar(x,y)

                #if enough pixels have been found to calculate a noise value for this pixels radius,
                # the pixel value is divided by the standard deviation of pixels at that radius
                    try:
                    #if statement prevents a divide by zero warning message
                        if (NoiseMap[radius] == 0):
                            indiv[x][y] = np.nan
                 
                        else:
                            indiv[x][y]/=NoiseMap[radius]

                        if makenoisemap == True:
                            noise[x][y] = NoiseMap[radius]

                    #if no noise value can be calculated, pixel is given a nan value
                    except:
                        indiv[x][y] = np.nan

                    #BEGIN PLANET STUFF
                    #if known or injected planets are present, mask them and measure their SNRs
                    if planets != False:

                        #returns whether this location has a planet and, if so, which planet in the list
                        isplanetpix, iscorepix, p = isPlanet(radius, angle, planets)

                        if checkmask==True:
                            if (isplanetpix == True):
                                msk[x][y] = 0
                    # captures an array with just the planet pixels
                    # use angular key equivalent to radial one (roughly mimic circular mask)
                        #planets_core = deepcopy(planets)
                        #planets_core[2][1] = np.arctan(planets_core[2][0]/planets[0][0])*180/np.pi
                    #print("test radial mask is", planets_core[2][0], "azimuthal mask is now", planets_core[2][1], "instead of", planets[2][1])

                        if isplanetpix == True:
                            planet_pixels[x][y][methodctr][s][p]=indiv[x][y]
                        ##includes only positive pixels under mask as planet pixels to avoid including self-subtraction regions in sums
                            if iscorepix == True:
                                if indiv[x][y]>0:
                                    planet_pixels_pos[x][y][methodctr][s][p]=indiv[x][y]
                                    npospix[methodctr][s][p]+=1

                        #count up how many pixels OUTSIDE the mask have >5 sigma values
                        if not isplanetpix:
                            notplanetpix += 1
                            if indiv[x][y] > sigma:
                                fivesig+=1

                        #count up how many pixels OUTSIDE the mask AND inside the control radius have >5 sigma values
                        if (radius>fwhm) and (radius<ctrlrad) and (np.isnan(indiv[x][y])==False) :
                            
                            if not isplanetpix:
                                if indiv[x][y] > sigma:
                                    fivesig_atmask += 1
                        #count up how many pixels INSIDE the mask have >5 sigma values
                        if (radius>fwhm):
                            
                            if isplanetpix:
                                allplanetpix += 1
                                if indiv[x][y] > sigma:
                                    fivesig_inmask += 1                     

                #store output for this method and # KL modes
                Output[methodctr,s,:,:] = indiv

                if checkmask==True:
                    msks[s,:,:]=msk
                if makenoisemap==True:
                    noises[methodctr, s,:,:]=noise

                #calculate and store planet data
                if planets != False:
                    for p in np.arange(nplanets):
                        snrs[methodctr,s,p]=np.nanmax(planet_pixels[:,:,methodctr,s,p])
                        snr_sums[methodctr,s,p] = np.nansum(planet_pixels_pos[:,:,methodctr,s,p])/npospix[methodctr,s,p]

                    snr_spurious[methodctr,s,:]=[fivesig, fivesig_atmask, fivesig_inmask, allplanetpix, notplanetpix]
                #print("max SNR under mask is", snrs[methodctr,s], "for slice", s)
                #print("sum of SNRs under mask is", snr_sums[methodctr,s], "for slice", s)
                    max_toprint = [round(n,5) for n in snrs[methodctr,s]]
                    head["MX"+method[0:2]+'_'+'KL'+str(s)]= str(max_toprint)
                    sums_toprint = [round(n,1) for n in snr_sums[methodctr,s]]
                    head["SM"+method[0:2]+'_'+'KL'+str(s)]= str(sums_toprint)
                    head["EX" + method[0:2] + '_' + 'KL'+str(s)] = str(snr_spurious[methodctr,s,:])

        methodctr += 1

    #saves output to disk if saveOutput designated True
    if (saveOutput == True):

        fits.writeto(outname+'.fits', Output, head, overwrite=True)
        print("Wrote", outname)

        if checkmask==True:
            fits.writeto(outname+'_corepix.fits', planet_pixels_pos, overwrite=True)  
            fits.writeto(outname+'_maskpix.fits', planet_pixels, overwrite=True)   
            #print(planet_pixels_pos.shape, planet_pixels.shape) 
            maskedims = msks*inp
            fits.writeto(outname+'_masked.fits', maskedims, head, overwrite=True)

        if makenoisemap==True:
            fits.writeto(outname+'_noisemap.fits', noises, head, overwrite=True)

    #returns final SNR map
    if planets != False:
        if checkmask==True:
            return Output, snrs, snr_sums, snr_spurious, maskedims
        else:
            return Output, snrs, snr_sums, snr_spurious
    
    #return only snrmap if no known or injected planets present
    else:
        return Output


def getPlanet_peak(snrmap, sep, pa, _range):
    #takes a 2d map - note have to loop over kl modes and modes

    #try:
    yDim, xDim = np.shape(snrmap)
    #print("check - cube shape is", snrmap.shape)
    #except:
       # yDim, xDim = np.shape(snrmap)
       # zDim = 1
    global XCenter
    global YCenter
    XCenter = (xDim-1)/2
    YCenter = (yDim-1)/2  

    #translate pa and sep into (y,x) for MagAO image geometry
    x = int(sep*math.cos(math.radians(pa+90))+XCenter)
    y = int(sep*math.sin(math.radians(pa+90))+YCenter)

    planet = -100000000

    #for mode in range(modedim):
        #for kl in range(kldim):
            #print("check - mode = ", mode, "kl slice = ", kl)
    for i in range (x-_range, x+_range):
        for j in range (y-_range, y+_range):
            #finds highest pixel under mask for this planet and returns it
            if (snrmap[j][i] > planet):
                planet = snrmap[j][i]
                
    if (planet == -100000000):
        return np.nan
    #print("planet SNR is", planet)
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


