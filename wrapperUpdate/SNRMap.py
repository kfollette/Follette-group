import numpy as np
import math
import statistics as stat 
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import sys
import time 

XCenter = 224.5
YCenter = 224.5

def nameOutput(filename, output):
    
    if (output == None):
        
        #if no output name is specified at runtime and input data is a fits file
        if(isinstance(filename, str)):
            outputName = filename[:-5] + "_SNRMap.fits"
        
        #if data type is not a string (file name), names output file after date and time 
        else:
            outputName = "SNRMap_"  +(time.strftime("%d-%m-%Y")) +"_" +(time.strftime("%H-%M-%S")) + ".fits"
        
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
    This function takes in the polar coordinates of a point to be tested and a touple containing lists of parameters for planets in the data to be masked.
    
    Reuired Inputs:
    1. Integer radius of point to be tested
    2. Angle coordinate of point to be tested
    3. Tuple containing the following lists:
        a. List of radial coordinates of planets in data
        b. List of corresponding position angles of planets in data (must be same length of a)
        c. List containing radial thickness of deired mask on either side of the planet, followed by the disired angular thickness
            
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



def stdevMap(indiv, planets):
    
    """
    This function takes a filename and a list of parameters for objects to mask and outputs a dictionary object of integer value radii pointing to the standard deviation for pixel values in the image at that radius from the center.
    
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
    #creates empty dictionary objects to store unmaked pixel values and radial standard deviations
    stdevs_ = {}
    radialProfs = {}
    
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
                else:
                    radialProfs[radius] = [indiv[x][y],]
        
                     
    
    #loops through each key in radial profile dictionary, and takes standard deviation of list of pixel values
    #adds standard deviation to stdevs_ dictionary with radius as the key
    #ignores data points if there are too few at a certain radius to take a standard deviation. These pixels will eventually become nans
    for r in radialProfs.keys():
        try: 
            stdevs_[r]= np.nanstd(radialProfs[r])
        except: 
            pass
        
    #returns dictionary holding standard deviations
    return stdevs_
    



def create_map(filename, planets = None, saveOutput = True, outputName = None):
    """
    creates signal to noise ratio map of image.
    
    Required Input:
    1. String containing filename of original klipped image OR object containing data already taken from original klipped image

    Optional Inputs:
    1. Tuple containing the following lists:
        a. List of radial coordinates of planets in data
        b. List of corresponding position angles of planets in data (must be same length of a)
        c. List containing radial thickness of deired mask on either side of the planet, followed by the disired angular thickness
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
    6/28/2017
    """
    
    #checks data type of 'filename'
    # if 'filename' is a string, assumes it is a filepath and reads in file
    if(isinstance(filename, str)):
        indiv = read_file(filename)
        
    #if data type is not a string, reads in python object holding data
    else:
        indiv = filename
        
    #creates dictionary holding the standard deviation of pixlel values at each radius 
    stdMap = stdevMap(indiv, planets)
  
    #gets size of pixel value array
    xDim, yDim = np.shape(indiv)  
    global XCenter
    global YCenter 
    XCenter = (xDim-1)/2
    YCenter = (yDim-1)/2

    #loops through all pixels in array
    for x in range (xDim): 
        for y in range (yDim):
            
            #converts indeces to polar coordinates
            radius, angle = toPolar(x,y,)
            
            #use for debugging if you want to see where the mask is:
            #if (isPlanet(radius, angle, planets)):
                #indiv[x][y] = np.nan
           
            #if enough pixels have been found to calculate a standard deviation for this pixels radius, the pixel value is divided by the standard deviation of pixels at that radius
            try:
                #if statement prevents a divide by zero warning message
                if (stdMap[radius] == 0):
                    indiv[x][y] = np.nan
                else:
                    indiv[x][y] = indiv[x][y]/stdMap[radius]
                
                #debugging step to show noise map:
                #indiv[x][y] = stdMap[radius]
     
     
            #if no standard deviation has been calculated, pixel is given a nan value
            except:
                indiv[x][y] = np.nan
    
    #saves output to disk if saveOutput designated True
    if (saveOutput == True):
        hdu = fits.PrimaryHDU(indiv)
        hdulist = fits.HDUList([hdu])
        newname = str(nameOutput(filename, outputName))
        hdulist.writeto(newname, overwrite=True)
        print("Wrote %s to "%newname + os.getcwd())


    #returns final SNR map            
    return indiv
    
    



def getPlanet(snrmap, rad, pa, _range):
    
    xDim, yDim = np.shape(snrmap)
    global XCenter
    global YCenter
    XCenter = (xDim-1)/2
    YCenter = (yDim-1)/2  

    x = int(rad*math.cos(math.radians(pa+90))+XCenter)
    y = int(rad*math.sin(math.radians(pa+90))+YCenter)
    print(x)
    print(y)


    planet = -100000000
   
    for i in range (x-_range, x+_range):
        for j in range (y-_range, y+_range):
            if (snrmap[i][j] > planet):
                planet = snrmap[i][j]
            
    return planet



#print(getPlanet('med_HD142527_8Apr14short_SDI_a7m3-5KLmodes.fits', 220, 215, 10))





