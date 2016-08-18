#Written by Elijah Spiro
#6/29/16 - 7/12/16
#Takes in H-Alpha and continuous spectrum data cubes. Outputs a single value of scale that will best suit the data.              
#Works on median data sets

print("NOW RUNNING 'GenerateMedianScale.py'")
               
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
    #debug = input("Run in debug mode? y/n ")
    debug = "y"

#Establish the directory we want to work from
def switchDirectories():
    print("Starting in directory: " + os.getcwd())
    os.chdir(str(sys.argv[1]))
    #os.chdir("../HD100546/12Apr14A/")
    #os.chdir("../HD142527/HD142527/8Apr14/MERGED_long_sets/")
    #os.chdir("../HD142527/HD142527/15May15/")
    print("Switched into working directory: " + os.getcwd())
    printDivider()

#Visually organize program flow
def printDivider():
    print("--------------------------------------------------------------------------------")

#Requires the existence of file "mask.fits" in working directory
def applyMask(preMaskData):
    print("[STEP 3 OF 14] Masking median images")
    maskedArray = np.zeros((450,450))
    hdulist = fits.open("MEDIAN_MASK.fits")
    maskData = hdulist[0].data
    for y in range(450):
        for x in range(450):
            maskedArray[x][y] = preMaskData[x][y] * maskData[x][y]
            #print (str(x) + " " + str(y) + " " + str(maskedArray[x][y]))
    return maskedArray

def applyMaskCube(preMaskDataCube, filetype):
    dataSize = len(preMaskDataCube)
    maskedArray = np.zeros((dataSize,450,450))
    hdulist = fits.open("CUBE_MASK.fits")
    maskData = hdulist[0].data
    for z in range(0, len(preMaskDataCube)):
    #for z in range(0,10):                                                                                                                                                                                  
        for y in range(0, len(preMaskDataCube[0])):
            for x in range(0, len(preMaskDataCube[0][0])):                                                                                                                                                                                                         
                maskedArray[z][y][x] = preMaskDataCube[z][y][x] * maskData[z][y][x]
        if debug == "y":
            if filetype == 1:
                print("[STEP 5 OF 14] Applying mask to image " + str(z) + " of " + str(len(preMaskDataCube)) + " (H-Alpha)")
            elif filetype == 2:
                print("[STEP 5 OF 14] Applying mask to image " + str(z) + " of " + str(len(preMaskDataCube)) + " (Continuous)")
    hdu = fits.PrimaryHDU(maskedArray)
    hduList = fits.HDUList([hdu])
    if filetype == 1:
        hduList.writeto("MASKED_DATA_HA.fits")
    elif filetype == 2:
        hduList.writeto("MASKED_DATA_CONT.fits")
    return maskedArray


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
        if debug == "y":
            print("Scanning image " + str(z))
        for y in range(450):
            for x in range(450):
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
    scaledArray = np.zeros((450,450,1758))
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
        print("[STEP 1 OF 14] Taking median of H-Alpha data")
    elif filetype == 2:
        print("[STEP 1 OF 14] Taking median of continuous spectrum data")
    medianImage = np.nanmedian(dataCube, axis=0)
    print("[STEP 1 OF 14] Median completed; shape of new data is " + str(medianImage.shape))
    if filetype == 1:
        medianPixelsHA = medianImage
        return medianImage
        print(medianPixelsHA[100][100])
    elif filetype == 2:
        medianPixelsCONT = medianImage
        return medianImage

#Assumes halo radius of 100 pixels, returns the total value of center 31,415 pixels
def calculateTotalStarlight(image):
    runningSum = 0
    for y in range(125,325):
        for x in range(125,325):
            if image[x][y] == image[x][y]:
                runningSum = runningSum + image[x][y]
    print("Total starlight identified: " + str(runningSum))
    return runningSum
             
#Takes 2D array and ratio, scales every pixel
def multiplyByRatio(array, ratio):
    scaledArray = np.zeros((450,450))
    ratio = ratio
    for y in range(0,450):
        for x in range(0,450):
            scaledArray[x][y] = (array[x][y] * ratio)
    return scaledArray
                     
#Takes H-Alpha and Continuous spectrum 2D arrays, subtracts each pixel between the two     
def subtractImage(HAlpha, Cont):
    subtractedArray = []
    for y in range(0,450):
        for x in range(0,450):
            subtractedArray.append(HAlpha[x][y] - Cont[x][y])
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
def step(guess, result, amountToStep,  medianImHA, medianImCONT, wrongTries, dataX, dataY, firstCall, flipped):
    if debug == "y":
        print("[STEP 6 OF 14] Flipped: " + str(flipped) + " | FirstCall: " + str(firstCall))
    if flipped == False:
        ratio = guess * (1 - amountToStep)
    elif flipped == True:
        ratio = guess * (1 + amountToStep)
    print("[STEP 6 OF 14] Testing with ratio: " + str(ratio))
    scaledCONT = multiplyByRatio(medianImCONT, ratio)
    subtractedValues = subtractImage(medianImHA, scaledCONT)
    subtractedValues = np.array(subtractedValues)
    subtractedValues = subtractedValues[~np.isnan(subtractedValues)]
    starlight = calculateTotalRemainingStarlight(subtractedValues)
    if starlight < result:
        wrongTries = 0
        print("[STEP 6 OF 14] Total starlight from new ratio: " + str(starlight) + "-- improvement, stepping again in this direction")
        dataX.append(ratio)
        dataY.append(starlight)
        printDivider()
        if flipped == False:
            step(ratio, starlight, .01, medianImHA, medianImCONT, wrongTries, dataX, dataY, False, False)
        elif flipped == True:
            step(ratio, starlight, .01, medianImHA, medianImCONT, wrongTries, dataX, dataY, False, True)
    elif starlight > result:
        if firstCall == False:
            if wrongTries < 10:
                print("[STEP 6 OF 14] Total starlight from new ratio: " + str(starlight) + "-- this is worse, trying again with smaller increment")
                dataX.append(ratio)
                dataY.append(starlight)
                printDivider()
                if flipped == False:
                    step(guess, result, (amountToStep/2), medianImHA, medianImCONT, (wrongTries+1), dataX, dataY, False, False)
                elif flipped == True:
                    step(guess, result, (amountToStep/2), medianImHA, medianImCONT, (wrongTries+1), dataX, dataY, False, True)
            elif wrongTries >= 10:
                print("[STEP 6 OF 14] Total starlight from new ratio: " + str(starlight) + "-- this is worse, final ratio is " + str(guess))
                idealRatio = guess
                idealStarlight = result
                #plotResults(dataX, dataY, idealRatio, idealStarlight)
                finalScale.append(idealRatio)
        elif firstCall == True:
            print("[STEP 6 OF 14] Stepped in the wrong direction initially. Trying again, in o\
ppositte direction")
            step(guess, result, amountToStep, medianImHA, medianImCONT, 0, dataX, dataY, False, True)

def createMaskFile(medianImHA):
    #Scan median image
    if debug == "y":
        print("[STEP 2 OF 14] Building median image mask")
    centerX = 225
    radius = 0
    continueScanning = True
    while continueScanning == True:
        if medianImHA[225][centerX + radius] > 15000:
            radius = radius + 1
        else:
            continueScanning = False
            finalRadius = radius
    #Write file
    data = np.zeros((450,450))
    for y in range(450):
        for x in range(450):
            if ((x - 225)**2 + (y - 225)**2) < 35**2:
                if ((x - 225)**2 + (y - 225)**2) > finalRadius**2:
                    data[y][x] = 1
                else:
                    data[y][x] = float('NaN')
            else:
                data[y][x] = float('NaN')
    hdu = fits.PrimaryHDU(data)
    hduList = fits.HDUList([hdu])
    hduList.writeto("MEDIAN_MASK.fits")

def createMaskFiles(cube):
    #Scan image cube                                                                                                                                                                                     
    centerX = 225                                                                                                                                                                                            
    data = np.zeros((len(cube),450,450))
    for z in range(len(cube)):
        if debug == "y":
            print("[STEP 4 OF 14] Scanning cube image " + str(z) + " of " + str(len(cube)) + " for mask creation")
        radius = 0
        continueScanning = True
        while continueScanning == True:
            if cube[z][225][centerX + radius] > 15000:
                radius = radius + 1
            else:
                continueScanning = False
                finalRadius = radius
        if debug == "y":
            print("[STEP 4 OF 14] Writing data to cube mask file (RADIUS = " + str(finalRadius) + " pixels)")
        for y in range(450):
            for x in range(450):
                if ((x - 225)**2 + (y - 225)**2) < 35**2:
                    if ((x - 225)**2 + (y - 225)**2) > finalRadius**2:
                        data[z][y][x] = 1
                    else:
                        data[z][y][x] = float('NaN')
                else:
                    data[z][y][x] = float('NaN')
    hdu = fits.PrimaryHDU(data)
    hduList = fits.HDUList([hdu])
    hduList.writeto("CUBE_MASK.fits")



def SDIMedianSR(scale, medianImHA, medianImCONT):
    print("[STEP 7 OF 14] Writing scaled median image SDI")
    SDIImage = np.zeros((450,450))
    for y in range(450):
        for x in range(450):
            scaledCONT = medianImCONT[x][y] * scale
            subtractedValue = medianImHA[x][y] - scaledCONT
            SDIImage[x][y] = subtractedValue
    hdu = fits.PrimaryHDU(SDIImage)
    hduList = fits.HDUList([hdu])
    prihdr = hduList[0].header
    prihdr['MEDSCL'] = scale
    hduList.writeto("SDI_MEDIAN_SR.fits")
                               

def SDICubeSR(scale, dataHA, dataCONT):
    dataSize = len(dataHA)
    SDIImage = np.zeros((dataSize,450,450)) #Z,Y,X 
    for z in range(len(SDIImage)):
        for y in range(450):
            for x in range(450):
                scaledCONT = dataCONT[z][y][x] * scale
                subtractedValue = dataHA[z][y][x] - scaledCONT
                SDIImage[z][y][x] = subtractedValue
        print("[STEP 8 OF 14] Writing scaled image " + str(z) + " of " + str(dataSize) + " to cube")
    hdu = fits.PrimaryHDU(SDIImage)
    hduList = fits.HDUList([hdu])
    prihdr = hduList[0].header
    prihdr['MEDSCL'] = scale
    hduList.writeto("SDI_CUBE_SR.fits")


def SDICubeMR(scales, dataHA, dataCONT):
    dataSize = len(dataHA)
    SDIImage = np.zeros((dataSize,450,450)) #Z,Y,X                                  
    for z in range(len(SDIImage)):
        for y in range(450):
            for x in range(450):
                scaledCONT = dataCONT[z][y][x] * scales[z]
                subtractedValue = dataHA[z][y][x] - scaledCONT
                SDIImage[z][y][x] = subtractedValue
        print("Building image " + str(z))
    hdu = fits.PrimaryHDU(SDIImage)
    hduList = fits.HDUList([hdu])
    hduList.writeto("SDI_CUBE_MR.fits")

###############################################################################
#                            BEGIN PROGRAM FLOW                               #
###############################################################################

#Prepare for calculations
debug() #Prompts user whether or not to display progress as program runs
switchDirectories() #Navigates to the proper working directory (NOTE you'll have to change this function to fit your machine)
#readFitsCubeToArray("Line_clip450_reg.fits", 1) #Data is stored in dataHA
#readFitsCubeToArray("Cont_clip450_reg.fits", 2) #Data is stored in dataCONT
readFitsCubeToArray(str(sys.argv[2]), 1)
readFitsCubeToArray(str(sys.argv[3]), 2)


#Obtain an initial ratio
medianImHA = median(dataHA, 1)
medianImCONT = median(dataCONT, 2)
createMaskFile(medianImHA)
medianImHA = applyMask(medianImHA)
medianImCONT = applyMask(medianImCONT)
createMaskFiles(dataHA)
dataHA = applyMaskCube(dataHA, 1)
dataCONT = applyMaskCube(dataCONT, 2)

totalStarlightHA = calculateTotalStarlight(medianImHA)
totalStarlightCONT = calculateTotalStarlight(medianImCONT)
ratio = totalStarlightHA / totalStarlightCONT
print("Initial best guess of ratio of H-Alpha to Continuous Spectrum is " + str(ratio))

ratio2 = ratio/2

#Test ratio
scaledCONT = multiplyByRatio(medianImCONT, ratio)
subtractedValues = subtractImage(medianImHA, scaledCONT)
subtractedValues = np.array(subtractedValues)
subtractedValues = subtractedValues[~np.isnan(subtractedValues)]
starlight = calculateTotalRemainingStarlight(subtractedValues)
print("Total starlight from initial ratio: " + str(starlight))
printDivider()

scaledCONT2 = multiplyByRatio(medianImCONT, ratio2)
subtractedValues2 = subtractImage(medianImHA, scaledCONT2)
subtractedValues2 = np.array(subtractedValues2)
subtractedValues2 = subtractedValues2[~np.isnan(subtractedValues2)]
starlight2 = calculateTotalRemainingStarlight(subtractedValues2)

finalScale = []

#Kick off recursion
step(ratio, starlight, .01, medianImHA, medianImCONT, 0, [], [], True, False)
step(ratio2, starlight2, .01, medianImHA, medianImCONT, 0, [], [], True, False)

output = (finalScale[0] + finalScale[1]) / 2

SDIMedianSR(output, medianImHA, medianImCONT)
SDICubeSR(output, dataHA, dataCONT)

