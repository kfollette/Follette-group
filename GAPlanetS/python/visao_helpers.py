"""
Contains helpers for the first couple image processing steps in the visao pipeline
Includes:
    visao_getimtypes
    requad_angles
    make_gaussian
    fitgaussian
    fits_dateconv

all these functions should automatically be imported into the visao_process file

"""


import glob
from astropy.io import fits
import numpy as np
import math
import scipy
import time
import os
from scipy.ndimage.filters import gaussian_filter
import image_registration as imreg
import scipy.ndimage.interpolation as sci
import scipy.optimize as optimize



"""
##################################################################################################################
NAME: visao_getimtypes

PURPOSE: Extracts image type and other headers from all the VisAO images in a directory

DESCRIPTION: Generates a list of files in a directory, which can be specified with the subdir keyword, which match the
        regex "V47_*.fits".  The prefix can be changes with the prefix keyword (e.g. to 'dsub_').  At minimum returns
        the filenames and the image types of the files (science, dark, etc.).  To avoid re-generating the list of files,
        and instead use the file names passed in as fnames, set the keyword usefnames.

INPUT KEYWORDS:
    prefix      :  the filename prefix, default is 'V47_'
    subdir      :  if set, operates on images in the specified directory, otherwise uses the current directory
    region      :  a 4 element vector specifying an area of the images to return. [x0,x1, y0,y1]
    goodims     :  a vector if indices for the images you want to retrieve, if not set all images are retrieved
    usefnames   : load the images specified in fnames, rather than getting all images in subdir
    double      : load images as doubles.  default is float.

OPTIONAL WORDS:     ['EXPTIME', 'VFOCPOS', 'AOLOOPST', 'ROTOFF', 'VGIMXPOS', 'VGIMYPOS', 'AVGWFE', 'STDWFE',
                    'VW1POSN', 'VFW2POSN', 'VFW3POSN', 'DATEOBS', 'AOREC', 'UT', 'AM', 'HA', 'FRAMENO', 'FGDMATS',
                    'ORIGXCEN', 'ORIGYCEN', 'GOODCENT', 'ROTANG']

OPTIONAL MAP KEYS:   {'ims': [], 'prefix': [], 'subdir': [], 'region': [], 'goodims': [],
                        'usefnames': [], 'double': []}

OUTPUTS:
    fnames      :   The file names of the images
    imtypes     :   Image types, 0=sci, 1=acq, 2=dark, 3=sky, 4=flat
    return_dictionary : dictionary for all other key words (ex return_dictionary['AOLOOPST'] gives all values for that

EXAMPLE:
    fnames = []
    imtypes = []
    visao_getimtypes(fnames, imtypes, 'EXPTIME', subdir = 'raw')

HISTORY:
    CREATOR: 2012-11-14 by Jared Males, jrmales@email.arizona.edu
    PY TRANS: 2016-07-06 by Wyatt Mullen, wmullen1@stanford.edu
"""
def visao_getimtypes(fnames=None, imtypes=None, *keywords, **keysMap):

    #for keyword in keywords:
    #probably some better way to do this, but not clear yet
    if not ('subdir' in keysMap):
        keysMap['subdir'] = '.'

    if not 'usefnames' in keywords or len(fnames) == 0:
        keysMap['usefnames'] = 0
        if not 'prefix' in keysMap:
            keysMap['prefix'] = 'V47_'
        srchstr = (keysMap['subdir'] + '/' + keysMap['prefix'] + '*.fits')
        srchstr = srchstr.strip()
        #searches for the files
        files = glob.glob(srchstr)
        fnames += files
    else:
        keysMap['usefnames'] = 1

    if 'goodims' in keysMap and keysMap['goodims'] > 0:
        fnames = []
        for image in keysMap['goodims']:
            fnames.append(image)

    #if 'ims' in keys:
        #Do a lot of stuff that I'm not going to try
        #inability to uses class 'mrdfits' hindered full functionality

    #Creates a dictionary that we will be able to return with all the information
    returnDictionary = {}
    for word in keywords:
        returnDictionary[word] = []

    #specifies number of images we will work with
    nims = len(fnames)
    for item in range (0,nims):
        #inability to use classes 'statusline' and 'mrdfits' here hindered full functionality

        fitsData = fits.open(fnames[item])
        head = fitsData[0]

        #helps specify the image type
        imageType = head.header['VIMTYPE'].strip()
        if imageType == 'SCIENCE': imtypes.append(0)
        if imageType == 'ACQUISITION': imtypes.append(1)
        if imageType == 'DARK': imtypes.append(2)
        if imageType == 'SKY': imtypes.append(3)
        if imageType == 'FLAT': imtypes.append(4)

        # should take care of all the different arrays created.
        # Will return a map of the keywords specified attached to their arrays
        for word in keywords:

            if word == 'AOLOOPST':
                if head.header['AOLOOPST'] == 'CLOSED':
                    returnDictionary['AOLOOPST'].append(1)
                else:
                    returnDictionary['AOLOOPST'].append(0)

            #converts dates from 2014-04-09T03:42:12.759794 to 2014099.1543143494 and puts them in array
            elif word == 'DATE-OBS':
                dobs = head.header[word]
                returnDictionary[word].append(fitsdate_conv(dobs))

            elif word == 'FGDMATS':
                ts = head.header[word]
                returnDictionary[word].append(fitsdate_conv(ts))

            #Writes different values to the array based on the gain (0 is high, 3 is low)
            elif word == 'V47GAIN':
                g = head.header[word]
                if g == 'HIGH': returnDictionary[word].append(0)
                if g == 'MEDHIGH': returnDictionary[word].append(1)
                if g == 'MED': returnDictionary[word].append(2)
                if g == 'LOW': returnDictionary[word].append(3)

            else:
                if word in head.header:
                    returnDictionary[word].append(head.header[word])

    #all values we need are in this "return dictionary"
    #if you need to access them in another function, use whatever dictionary you have assigned to the function
    #call the word you needed and get the list or whatever
    return returnDictionary

"""
##################################################################################################################
NAME: requad_angles

PURPOSE/DESCRIPTION: Adjusts the vector of angles such that it is continuous if it wraps around 0 or 360 by moving
                    a portion of it to be lt 0

INPUTS: angles   :  vector of angles Q such that 0<= Q < 360 (2pi).  Must be monotonic, and have less than 360 (2pi)
                    degree total change

KEYWORDS: radians    :  if set, then angles are in radians with range 0<=0<2pi

OUTPUT: returns modified vector

HISTORY:
    CREATOR: 2013-01-22 by Jared Males, jrmales@email.arizona.edu
    PY TRANS: 2016-07-07 by Wyatt Mullen, wmullen1@stanford.edu
"""
def requad_angles(angles, *keyWords):

    subval = 360.
    if 'radians' in keyWords: subval = 2. * math.pi

    #don't modify angles
    rangles = angles

    #finds the first time the angle changes in the data
    ind = 1
    while (ind < len(angles) - 1) and angles[ind] == angles[0]:
        ind += 1

    if angles[ind] > angles[0]:
        min_angle = min(angles)
        if min_angle < angles[0]:
            for x in range (0,len(angles)):
                if angles[x] >= angles[0]:
                    rangles[x] = angles[x] - subval
    else:
        max_angle = max(angles)
        if max_angle > angles[0]:
            for x in range (0,len(angles)):
                if angles[x] >= angles [0]:
                    rangles[x] = angles[x] - subval

    return rangles


"""
##################################################################################################################
NAME: fitsdate_conv

PURPOSE: Converts a string date in fits format to a double

INPUT: datestr  : the input date in format 'YYYY-MM-DDTHH:MM:SS.SS...'

OUTPUT: returns the date in double precision days

HISTORY:
    CREATOR: 2013-01-21 by Jared Males, jrmales@email.arizona.edu
    PY TRANS: 2016-07-06 by Wyatt Mullen, wmullen1@stanford.edu
"""
def fitsdate_conv(date):

    regYear = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    leapYear = [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    year = int(date[0:4])
    month = int(date[5:7])
    day = int(date[8:10])
    hour = int(date[11:13])
    minute = int(date[14:16])
    second = float(date[17:])
    partDay = 0

    if year % 4 == 0:
        for x in range (1,month):
            partDay += leapYear[x]
    else:
        for x in range (1,month):
            partDay += regYear[x]

    day += partDay

    newDate = year * 1000 + day
    partialDay = ((((second / 60) + minute) / 60) + hour) / 24
    newDate += partialDay

    return newDate


"""
##################################################################################################################
NAME: make_Gaussian

Purpose: Creates a gaussian array of the specified dimensions, fwhm, and center

DESCRIPTION: Creates an (x,y) array with values ranging from 1 to x and 1 to y. Uses the definition of a gaussian
        to then populate that array with gaussian spread values

INPUTS:
    dimX:   x dimension (width) of gaussian
    dimY:   y dimension (length) of gaussian
    fwhm:   full width half max of gaussian, if not specified automatic value is 10
    center: (x,y) coordinates of center of gaussian

OUTPUTS:
    2-d array of gaussian values

HISTORY:
    Outline: http://stackoverflow.com/questions/7687679/how-to-generate-2d-gaussian-with-python
    Modified: 2016-07-13 by Wyatt Mullen, wmullen1@stanford.edu
"""
def make_Gaussian(dimX, dimY, fwhm = 10, center=None):

    x = np.arange(0, dimX, 1, float)
    y = np.arange(0, dimY, 1, float)
    y = y[:,np.newaxis]

    if center is None:
        x0 = dimX / 2
        y0 = dimY / 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)




"""
##################################################################################################################
NAME: 2D Gaussian Functions: gaussian, moments, fitgaussian (main function)

PURPOSE/DESCRIPTION: Fits a 2 gaussian to the given array, cannot be a fits image

INPUT: use fitgaussian
    data:   an image array, preferably centered with a 2d gaussian pattern

OUTPUT:
    x/y cen:    coordinates of the center
    x/y fwhm:   full width half maxes in each dimension
    height:      max value of the gaussian fit

HISTORY:
    SOURCE: http://scipy-cookbook.readthedocs.io/items/FittingData.html
"""
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p