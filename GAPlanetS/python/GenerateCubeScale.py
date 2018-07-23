#Written by Elijah Spiro
#6/29/16 - 11/29/16
#Takes in H-Alpha and continuous spectrum data cubes. Outputs a file of scale values that will best suit the data.              
#Works on data cube sets               
#Modified Aug 2017 by KBF to accomodate arbitrary image size, add overwrite=True to fits output, remove masking (made obsolete by specification of zones with arbitrary radii), rename some files

#TODO:
#first plot includes a bunch of things that are not showing up
#way too many print functions to terminal. should be removed or relaces with statuslines.
#debug and radii should be required input to function call rather than user-prompted input
#Add wavefront error to ratio plot (plot1())

print("NOW RUNNING 'GenerateCubeScale.py'")

##############################################################################
#                             IMPORT LIBRARIES                               #
##############################################################################
                  
import os
import sys
import math
from astropy.io import fits 
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import statistics

###############################################################################
#                    FUNCTION AND VARIABLE DECLARATIONS                       #
###############################################################################

#Arrays to hold the 1758 brightest pixels
maxStarValuesHA = []
maxStarValuesCONT = []
scaledMaxValues1 = []
medianPixelsHA = []
medianPixelsCONT = []

#Prompt for debug info
def debug():
    global debug
    global InnerRadius
    global OuterRadius
    debug = input("Run in debug mode? y/n ")
    #debug = "y"
    InnerRadius = int(input("Enter an inner radius (distance from center)"))
    OuterRadius = int(input("Enter an outer radius (distance from center)"))

#Establish the directory we want to work from
def switchDirectories():
    print("Starting in directory: " + os.getcwd())
    global argnum 
    argnum = 0
    pathToFiles = sys.argv[1]
    try:
        while (not pathToFiles[-5:] == '%done'):
            argnum += 1
            pathToFiles = pathToFiles + " " + sys.argv[1+argnum]
        pathToFiles = pathToFiles[:-5]
        
    except:
        pathToFiles = sys.argv[1]
        argnum = 0
    os.chdir(pathToFiles)
    print("Switched into working directory: " + os.getcwd())
    printDivider()

#Visually organize program flow
def printDivider():
    print("--------------------------------------------------------------------------------")

"""
#Made obsolete by specification of arbitrary radii
#Requires the existence of file "mask.fits" in working directory
def applyMask(preMaskDataCube):
    cubeSize = len(preMaskDataCube)
    dim = preMaskDataCube.shape[1]
    maskedArray = np.zeros((cubeSize,dim,dim))
    hdulist = fits.open("CUBE_MASK.fits")
    maskData = hdulist[0].data
    for z in range(0, len(preMaskDataCube)):
    #for z in range(0,10):
        for y in range(0, dim):
            for x in range(0, dim):
                #print("preMaskDataCube[" + str(x) + "][" +str(y) + "][" + str(z) + "] = " + str(preMaskDataCube[x][y][z])  +  "     |     maskData[" + str(x) + "][" + str(y) + "] = " + str(maskData[x][y]))
                maskedArray[z][y][x] = preMaskDataCube[z][y][x] * maskData[z][y][x]
        if debug == "y":
            print("Applying mask to image " + str(z) + " of " + str(len(preMaskDataCube)))
    return maskedArray
"""

#Create 3D array from fits file cube
def readFitsCubeToArray(filename, filetype):
    global dataHA
    global dataCONT
    print("Now reading in data cube \'"+str(filename)+"\'")
    hdulist = fits.open(filename)
    hdulist.info()
    if filetype == 1:
        dataHA = hdulist[0].data
    elif filetype == 2:
        dataCONT = hdulist[0].data
    printDivider()
    #NOTE: Array is in format data[Z][Y][X]

#Identifies the maximum value pixel in every file within a data cube
def findMaxPixel(data, filetype):
    numFiles = 10
    if filetype == 1:
        print("Scanning " + str(numFiles)  +  " images in Hydrogen-Alpha image cube to evaluate brightest pixels")
    elif filetype == 2:
        print("Scanning " + str(numFiles) + " images in continuous spectrum image cube to evaluate brightest pixels") 
    count = 0
    for z in range(numFiles):
        maxVal = 0
        tempX = 0
        tempY = 0
        dim=data.shape[1]
        if debug == "y":
            print("Scanning image " + str(z))
        for y in range(dim):
            for x in range(dim):
                if data[z][y][x] > maxVal:
                    tempX = x
                    tempY = y
                    maxVal = data[z][y][x]
        if debug == "y":
            print("Max value found is " + str(maxVal) + " at x=" + str(tempX) + " , y=" + str(tempY))
        if filetype == 1:
            maxStarValuesHA.append(maxVal)
        elif filetype == 2:
            maxStarValuesCONT.append(maxVal)
        elif filetype == 3:
            scaledMaxValues1.append(maxVal)
        count = count + 1
    print("Done. " + str(count) + " brightest pixels added to array")
    printDivider()

#Generates a plot of H-Alpha vs Continuous Spectrum maximum pixel values
def plotPeakValues():
    print("Generating peak pixel value chart")
    plt.figure(1)
    plt.subplot(211)
    plt.plot(maxStarValuesHA, 'ro')    
    
    plt.subplot(212)
    plt.plot(maxStarValuesCONT, 'bo')
    
    ratioValues = []
    for i in range(0,len(maxStarValuesHA)):
        ratio = maxStarValuesHA[i] / maxStarValuesCONT[i]
        ratioValues.append(ratio)
        if debug == "y":
            print("Ratio of image " + str(i+1)  + " measured: "  + str(ratio))
    plt.xlabel("Image Number")
    plt.ylabel("Brightness")
    plt.show()
    scaleBySinglePixelRatio(ratioValues)

def scaleBySinglePixelRatio(scaleArray):
    #scaledArray is the maxStarValuesCONT array scaled UP to "match" HA values
    dim=dataCONT.shape[1]
    scaledArray = np.zeros((dim,dim,1758))
    try:
        for z in range(0, len(dataCONT)):
            for y in range(0, len(dataCONT[z])):
                for x in range(0, len(dataCONT[z][y])):
                    scaledValue = dataCONT[z][y][x] * scaleArray[z]
                    scaledArray[z][y][x] = scaledValue      
                    if debug == "y":
                        print(str(x) + " " + str(y) + " "  +  str(z) +" | " + str(dataCONT[x][y][z]) + " --> "  + str(scaledValue))                           
    except:
        print("done")

#Takes in a 3D array and spits out a 2D array of median values
def median(dataCube, filetype):
    if filetype == 1:
        print("Taking median of H-Alpha data")
    elif filetype == 2:
        print("Taking median of continuous spectrum data")
    medianImage = np.median(dataCube, axis=0)
    print("Median completed; shape of new data is " + str(medianImage.shape))
    if filetype == 1:
        medianPixelsHA = medianImage
        return medianImage
        print(medianPixelsHA[100][100])
    elif filetype == 2:
        medianPixelsCONT = medianImage
        return medianImage


def radarr(xdim,ydim):
    '''
    make an image of size xdim x ydim where every pixel has a value equal to its distance from the center of the array
    '''
    im = np.zeros((xdim, ydim))
    xcen = (xdim - 1) / 2.
    ycen = (ydim - 1) / 2.
    for index1 in range(xdim):
        for index2 in range(ydim):
            im[index1, index2]=np.sqrt((xcen-index1)**2. + (ycen-index2)**2.)
    #hdu=fits.PrimaryHDU(im)
    #hduList = fits.HDUList([hdu])
    #hduList.writeto('radarr_test.fits', overwrite=True)        
    return(im)

 
def annularmask(xdim, ydim, rin, rout):
    '''
    makes a xdim x ydim mask where pixels with distances from center between rin and
    rout have value 1, and everything else has been replaced by NaNs
    '''
    arr = radarr(xdim, ydim)
    new = np.ones((xdim, ydim))
    for y in range(0,ydim):
        for x in range(0,xdim):
            if arr[y][x] <= rin:
                new[y][x] = np.nan
            elif arr[y][x] >= rout:
                new[y][x] = np.nan
    #hdu=fits.PrimaryHDU(new)
    #hduList = fits.HDUList([hdu])
    #hduList.writeto('ctrlmask_test.fits', overwrite=True)
    return(new)
    

#Assumes halo radius of 100 pixels, returns the total value of center 31,415 pixels
def calculateTotalStarlight(image, i):
    dim=image.shape[1]
    mask = annularmask(dim, dim, InnerRadius, OuterRadius)
    maskedim = image[i]*mask
    sum = np.nansum(maskedim)
    print("Total starlight identified: ", sum)
    return sum

def calculateTotalStarlightSquared(image, i, cube):
    runningSum = 0
    dim=image.shape[1]
    cen = int((dim-1)/2)
    for y in range(cen-100,cen+100):
        for x in range(cen-100,cen+100):
            value = 0
            if cube == True:
                if image[i][y][x] == image[i][y][x]:
                    value = image[i][y][x]**2
                    value = math.sqrt(value)
            elif cube == False:
                if image[y][x] == image[y][x]:
                    value = image[y][x]**2
                    value = math.sqrt(value)
            runningSum = runningSum + value
    #print("Total starlight identified: " + str(runningSum))
    return runningSum
             
#Takes 2D array and ratio, scales every pixel
def multiplyByRatio(array, ratio, i):
    dim=array.shape[1]
    scaledArray = np.zeros((dim,dim))
    ratio = ratio
    for y in range(0,dim):
        for x in range(0,dim):
            scaledArray[y][x] = (array[i][y][x] * ratio)
    return scaledArray
                     
#Takes H-Alpha and Continuous spectrum 2D arrays, subtracts each pixel between the two     
def subtractImage(HAlpha, Cont, i):
    dim=HAlpha.shape[1]
    print('test', HAlpha.shape[1])
    subtractedArray = []
    for y in range(0,dim):
        for x in range(0,dim):
            subtractedArray.append(HAlpha[i][y][x] - Cont[y][x])
    return subtractedArray

#Creates histogram of pixel count vs. number of pixels
def plotHistogram(subtractedValues, ratio):
    print("Plotting histogram of subtracted values")
    plt.grid(True)
    xbins3 = np.linspace(-10,10,200)
    plt.hist(subtractedValues, bins=xbins3, color='purple', alpha=.75)
    plt.xlabel("Pixel Value")
    plt.ylabel("Number of Pixels")
    plt.title("H-Alpha : Continuous Spectrum = " + str(ratio))
    plt.show() 
  
#Creates graph of subtracted data-- used for analysis and testing of scaling factor
def plotFlatline(subtractedValues, ratio):
    print("Plotting flatline of subtracted values")
    plt.grid(True)
    plt.title("H-Alpha : Continuous Spectrum = " + str(ratio))
    plt.ylabel("Pixel Value")
    plt.xlabel("Pixels in Sequential Order (Center implies starlight)")
    plt.plot(subtractedValues, 'g.')
    plt.show()

#Calculates total value of 1D array (accounts for negative numbers by making everything positive)
def calculateTotalRemainingStarlight(array):
    runningSum = 0
    for i in range(len(array)):
        tempVal = array[i]**2
        tempVal = math.sqrt(tempVal)
        runningSum = runningSum + tempVal
    return runningSum

#Generate plot of scaling factor vs. starlight remaining after SDI; used to visualize results
def plotResults(dataX, dataY, idealRatio, idealStarlight):
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.scatter(dataX, dataY, color="red")
    plt.ylabel("Total Remaining Starlight (Counts)")
    plt.xlabel("Scale Factor")
    ax.set_xlim([1.15,.95])
    plt.plot(dataX, np.poly1d(np.polyfit(dataX, dataY, 100))(dataX))
    #ax.annotate('(%s)' % str(idealRatio), xy=(idealRatio, idealStarlight), xytext=(3, 1.5), arrowprops=dict(facecolor='black', shrink=0.05))
    plt.grid(True)
    plt.title("Scale Factor vs. Remaining Pixel Count (" + str(idealRatio) + ")")
    plt.show()

#Recursive function to generate best guess of scaling ratio. Parameters described below:
#guess: previous ratio used                                                    
#result: previous result obtained
#amountToStep: amount to increment by (used in form: newRatio = oldRatio * (1 - amountToStep))
#medianImHA + medianImCONT: initial data files passed along to be reused each step in ratio testing                    
#wrongTries: keeps track of how many times the program has overstepped to shut off after 10 consecutive bad guesses
#dataX + dataY: arrays of data points to produce final plot                                                                                                                       
def step(guess, result, amountToStep,  medianImHA, medianImCONT, wrongTries, dataX, dataY, firstCall, flipped, i):
    if debug == "y":
        print("[STEP 10 OF 14] Flipped: " + str(flipped) + " | FirstCall: " + str(firstCall))
    if flipped == False:
        ratio = guess * (1 - amountToStep)
    elif flipped == True:
        ratio = guess * (1 + amountToStep)
    print("[STEP 10 OF 14] Testing with ratio: " + str(ratio))
    scaledCONT = multiplyByRatio(medianImCONT, ratio, i) 
    subtractedValues = subtractImage(dataHA, scaledCONT, i)
    subtractedValues = np.array(subtractedValues)
    subtractedValues = subtractedValues[~np.isnan(subtractedValues)]
    starlight = calculateTotalRemainingStarlight(subtractedValues)
    if starlight < result:
        wrongTries = 0
        print("[STEP 10 OF 14] Total starlight from new ratio: " + str(starlight) + "-- improvement, stepping again in this direction")
        dataX.append(ratio)
        dataY.append(starlight)
        printDivider()
        if flipped == False:
            step(ratio, starlight, .01, medianImHA, medianImCONT, wrongTries, dataX, dataY, False, False, i)
        elif flipped == True:
            step(ratio, starlight, .01, medianImHA, medianImCONT, wrongTries, dataX, dataY, False, True, i)
    elif starlight > result:
        if firstCall == False:
            if wrongTries < 10:
                print("[STEP 10 OF 14] Total starlight from new ratio: " + str(starlight) + "-- this is worse, trying again with smaller increment")
                dataX.append(ratio)
                dataY.append(starlight)
                printDivider()
                if flipped == False:
                    step(guess, result, (amountToStep/2), medianImHA, medianImCONT, (wrongTries+1), dataX, dataY, False, False, i)
                elif flipped == True:
                    step(guess, result, (amountToStep/2), medianImHA, medianImCONT, (wrongTries+1), dataX, dataY, False, True, i)
            elif wrongTries >= 5:
                print("[STEP 10 OF 14] Total starlight from new ratio: " + str(starlight) + "-- this is worse, final ratio is " + str(guess))
                idealRatio = guess
                idealStarlight = result
                finalRatios.append(idealRatio)
                #plotResults(dataX, dataY, idealRatio, idealStarlight)    
        elif firstCall == True:
            print("[STEP 10 OF 14] Stepped in the wrong direction initially. Trying again, in o\
ppositte direction")
            step(guess, result, amountToStep, medianImHA, medianImCONT, 0, dataX, dataY, False, True, i)


def obtainRatios(dataHA, dataCONT, i):
    ratio = dataHA / dataCONT
    return ratio


def makeFitsOfRatios(data):                                                                                                                                                                     
    print("[STEP 11 OF 14] Generating ratio list")
    hdu = fits.PrimaryHDU(data)
    hduList = fits.HDUList([hdu])
    hduList.writeto("scale_factors.fits", overwrite=True)

def SDICubeMR(scales, dataHA, dataCONT):
    dataSize = len(dataHA)
    dim=dataHA.shape[1]
    SDIImage = np.zeros((dataSize,dim,dim)) #Z,Y,X                                                                                                                                                         
    for z in range(len(SDIImage)):
        for y in range(dim):
            for x in range(dim):
                scaledCONT = dataCONT[z][y][x] * scales[z]
                subtractedValue = dataHA[z][y][x] - scaledCONT
                SDIImage[z][y][x] = subtractedValue
        print("[STEP 12 OF 14] Building image " + str(z) + " of " + str(dataSize) +  " in data cube")
    hdu = fits.PrimaryHDU(SDIImage)
    hduList = fits.HDUList([hdu])
    hduList.writeto("SDI_CUBE_INDIVRATIOS.fits", overwrite=True)

def plot1(errorData, medianScale):
    print("[STEP 13 OF 14] Producing ratio plot-- please close graph window to continue")
    hdulist = fits.open("scale_factors.fits")
    ratioData = hdulist[0].data
    medianData = []
    dataSize = len(errorData)
    for i in range(dataSize):
        medianData.append(medianScale)
    fig = plt.figure(1)
    fig = plt.subplot(111)
    plt.plot(ratioData, 'b-')
    plt.plot(medianData, 'r-', linewidth=2.0)
    plt.plot(errorData, 'y-')
    plt.legend(['Individual Ratios', 'Median Data Ratio', 'Wavefront Error'], loc='upper left')
    fig.set_facecolor('black')
    plt.ylabel("Scale Factor")
    plt.xlabel("Image Number")
    plt.title("Ratio Plot")
    plt.show()

def plot2():
    print("[STEP 14 OF 14] Producing starlight plot-- please close graph window to finish program execution")
    hdulist = fits.open("SDI_CUBE_INDIVRATIOS.fits")
    ratioCubeMR = hdulist[0].data
    hdulist = fits.open("SDI_CUBE_FIXEDSCL.fits")
    ratioCubeSR = hdulist[0].data
    hdulist = fits.open("SDI_MEDIAN_FIXEDSCL.fits")
    ratioMedianSR = hdulist[0].data
    array1 = []
    array2 = []
    array3 = []
    dataSize = len(ratioCubeSR)
    for i in range(dataSize):
        array1.append(calculateTotalStarlightSquared(ratioCubeMR, i, True))
        array2.append(calculateTotalStarlightSquared(ratioCubeSR, i, True))
        array3.append(calculateTotalStarlightSquared(ratioMedianSR, 0, False))
        print("Done appending " + str(i))
    fig = plt.figure(1)
    fig = plt.subplot(111)
    plt.plot(array1, 'b-')
    plt.plot(array2, 'y-')
    plt.plot(array3, 'r-', linewidth=2.0)
    plt.legend(['SDI_CUBE_MR', 'Cube median scale', 'Median image scale'], loc='upper left')
    fig.set_facecolor('black')
    plt.ylabel("Starlight Remaining (Counts)")
    plt.xlabel("Image Number")
    plt.title("Starlight Plot")

    #median = np.nanmedian()

    plt.show()

def getWaveFrontErrors():
    errors = []
    os.chdir("raw/")
    for filename in os.listdir(os.getcwd()):
        if debug == "y":
            print("Opening " + str(filename))
        hdulist = fits.open(filename)
        error = hdulist[0].header['AVGWFE']
        if debug == "y":
            print("Average wavefront error: " + str(error))
        errors.append(error/90)
    os.chdir("../")
    return errors

def getMedianScale():
    hdulist = fits.open("SDI_MEDIAN_FIXEDSCL.fits")
    scale = hdulist[0].header['MEDSCL']
    return scale

###############################################################################
#                            BEGIN PROGRAM FLOW                               #
###############################################################################

#Prepare for calculations
debug() #Prompts user whether or not to display progress as program runs
argnum = 0
switchDirectories() #Navigates to the proper working directory (NOTE you'll have to change this function to fit your machine)
readFitsCubeToArray(str(sys.argv[2+argnum]), 1)
readFitsCubeToArray(str(sys.argv[3+argnum]), 2)


#try:
#    errors = getWaveFrontErrors()
#except:
errors = []

"""
#Applies masks to each image in both data cubes 
print("Applying mask to H-Alpha data")
dataHA = applyMask(dataHA)
print("Finished masking H-Alpha data")
printDivider()
print("Applying mask to Continuum data")
dataCONT = applyMask(dataCONT)
print("Finished masking Continuum data")
printDivider()

#print("[STEP 9 OF 14] Reading in masked data files")
#hdulist = fits.open("MASKED_DATA_HA.fits")
#dataHA = hdulist[0].data
##hdulist.writeto('test.fits', overwrite=True)
#hdulist = fits.open("MASKED_DATA_CONT.fits")
#dataCONT = hdulist[0].data
"""

finalRatios = []
dataSize = len(dataHA)

#for i in range(len(dataHA[0][0])):
for i in range(dataSize):
    for time in range(5):
        printDivider()
    print("NOW RUNNING IMAGE " + str(i))
    for time2 in range(5):
        printDivider()
    #Obtain initial ratios
    totalStarlightHA = calculateTotalStarlight(dataHA, i)
    totalStarlightCONT = calculateTotalStarlight(dataCONT, i)
    ratio = obtainRatios(totalStarlightHA, totalStarlightCONT, i)
    if debug == "y":
        print("Initial best guess of image " + str(i)  + "'s ratio of H-Alpha to Continuous Spectrum is " + str(ratio))
    #Test ratio
    scaledCONT = multiplyByRatio(dataCONT, ratio, i)
    subtractedValues = subtractImage(dataHA, scaledCONT, i)
    subtractedValues = np.array(subtractedValues)
    subtractedValues = subtractedValues[~np.isnan(subtractedValues)]
    starlight = calculateTotalRemainingStarlight(subtractedValues)
    if debug == "y":
        print("Total starlight from image " + str(i)  +  "'s initial ratio: " + str(starlight))
    printDivider()
    #Kick off recursion
    step(ratio, starlight, .01, dataHA, dataCONT, 0, [], [], True, False, i)

for i in range(len(finalRatios)):
    print(str(i) + ": " + str(finalRatios[i]))

makeFitsOfRatios(finalRatios)

SDICubeMR(finalRatios, dataHA, dataCONT)

#for i in range(10):
 #   print()
#print(finalRatios)

plot1(errors, getMedianScale()) #Ratio plot                                     
plot2() #Starlight plot 


