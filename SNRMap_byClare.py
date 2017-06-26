
# coding: utf-8

# In[33]:

#import SNRMap_helper as snr
import numpy as np
import math
import statistics as stat 


# In[34]:

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


# In[35]:

def toPolar(x, y):
    r = int(((xCen-x)**2+(yCen-y)**2)**.5)
    theta = math.degrees(math.atan((y-yCen)/(x-xCen)))
    if(x < xCen): 
         theta = theta + 180
    elif(y< yCen):
          theta = theta + 360
    return (r,theta)


# In[36]:

def stdevMap(indiv, r1, r2, theta1, theta2):
    stdevs_ = {}
    radialProfs = {}
    xDim, yDim = np.indiv.shape()
    xCen = int(xDim/2)
    yCen = int(yDim/2)
    for x in range (xDim): 
        for y in range (yDim):         
            radius, angle = toPolar(x,y)   
            if((radius >= r1 or radius <= r2) and (theta1 >= angle or theta2 <= angle)):
                if (radialProfs.has_key(radius)):
                    radialProfs[radius].append(indiv[x][y])
                    
    for r in radialProfs.keys():
        stdevs_[r]= stat.stdev(radialProfs[r])
        
    return stdevs_
    


# In[41]:

def create_map(filename, r1=0, r2=0, theta1=0, theta2=0, saveOutput=false):
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
        
    
    stdMap = stdevMap(indiv, r1, r2, theta1, theta2)
    
    xDim, yDim = np.indiv.shape()
    
    for x in range (xDim): 
        for y in range (yDim):
            radius, angle = toPolar(x,y)
            indiv[x][y] = indiv[x][y]/stdMap[radius]
            
    return indiv
    
    


# In[42]:

def getPlanet(filename, x, y, range):
    stdMap = create_map(filename)
    
    planet = -100000000
    
    for i in range (x-range, x+range):
        for j in range (y-range, y+range):
            if (stdMap[i][j] > planet):
                highestPixel = planet[i][j]
                
    return planet


# In[ ]:



