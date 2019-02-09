#import os
#import sys
#import math
from astropy.io import fits
from astropy.modeling import models, fitting
import numpy as np
import matplotlib.pyplot as plt
#import glob
from scipy.optimize import curve_fit
from scipy.integrate import quad
#import pandas as pd
from numpy import sqrt
import image_registration as ir

def readFitsCubeToArray(filename):
    original = fits.getdata(filename)
    print("read fits file into data")
    medianImage = np.nanmedian(original, axis=0)
    data = medianImage
    #fits.writeto("test_median.fits",data,overwrite=True)
    print("created median")
    return original, data

def findGhostSingleIm(gX,gY, data):
    """
    takes gX and gY, guesses for the location of the ghost image
    returns the actual coordinates of the approximate ghost center 
    and the pixel value at that point
    """
    #xGhost = 383
    #yGhost = 217
    maxVal = 0
    tempX = 0
    tempY = 0
    
    #finds center based on brightest pixel around location of ghost
    for y in range(gY-15, gY+16):
        for x in range(gX-15,gX+16):
            if data[y,x] > maxVal:
                tempX = x
                tempY = y
                maxVal = data[y,x]
    print("Calculated ghost center is " + str(maxVal) + " at x=" + str(tempX) + " , y=" + str(tempY))
    print("returning")
    return (tempX, tempY, maxVal)

def makeSquare(xc,yc,data,filename="test_square.fits"):
    foundX, foundY, maxVal = findGhostSingleIm(xc,yc,data)
    maskedArray = np.zeros((31,31))
    for y in range(-15,16):
        for x in range(-15,16):
            maskedArray[y+15][x+15] = data[foundY+y][foundX+x]
    fits.writeto(filename, maskedArray, overwrite=True)
    return maskedArray, foundX, foundY

def moffat(xcen,ycen,amp,wid,power):
    return models.Moffat2D(amplitude=amp, x_0=xcen, y_0=ycen, gamma=wid, alpha=power)

def makeMoffatGhost(array,xcen,ycen,amp,wid,power):
    global best_vals2
    y,x = np.mgrid[:31,:31]
    z = array
    #xcen,ycen,mv = findGhostSingleIm(gx,gy)
    model = moffat(xcen,ycen,amp,wid,power)
    print(model.x_0)
    lfit = fitting.LevMarLSQFitter()
    fit = lfit(model, x, y, z)
    '''
        plt.figure(figsize=(8, 2.5))
        plt.subplot(1, 3, 1)
        plt.imshow(z, origin='lower', interpolation='nearest', vmin=-1e4, vmax=5e4)
        plt.title("Data")
        plt.subplot(1, 3, 2)
        plt.imshow(fit(x, y), origin='lower', interpolation='nearest', vmin=-1e4, vmax=5e4)
        plt.title("Model")
        plt.subplot(1, 3, 3)
        plt.imshow(z - fit(x, y), origin='lower', interpolation='nearest', vmin=-1e4, vmax=5e4)
        plt.title("Residual")
    '''
    print("x center is " + str(fit.x_0))
    print("y center is " + str(fit.y_0))
    return fit.x_0, fit.y_0

def normalizeSquare(array):
    maxVal = 0
    for y in range(31):
        for x in range(31):
            if array[y][x] > maxVal:
                maxVal = array[y][x]
    for y in range(31):
        for x in range(31):
            array[y][x] = array[y][x] / maxVal
            
def shift(x,y,original):
    size = original.shape[0]
    shiftX, shiftY = x-15, y-15
    for z in range(size):
        shifted_image = ir.fft_tools.shift2d(original[z],-shiftX,-shiftY)
        original[z] = shifted_image
        if(z%100 == 0):
            print("shift number " + str(z))  

################

def ghostIsolation(filename, gx, gy, amp, wid, power):

    original, data = readFitsCubeToArray(filename)
    maskedArray, foundX, foundY = makeSquare(gx,gy,data)

    #convert to maskedArray's coordinates
    foundXSmall = foundX+15-gx
    foundYSmall = foundY+15-gy

    xCenter, yCenter = makeMoffatGhost(maskedArray,foundXSmall,foundYSmall,amp,wid,power)
    print(xCenter, yCenter)

    shift(xCenter,yCenter, original)

    medianImage = np.nanmedian(original, axis=0)
    data = medianImage
    print("created median after shifting")
    maskedArray2, a, b = makeSquare(foundX,foundY,data,"ghost.fits")
    #normalizeSquare(maskedArray2)
    #fits.writeto("test_final.fits", maskedArray2, overwrite=True)
    
    '''
    with open("moffat.txt","w") as text_file:
        text_file.write(str(amp) + " " + str(ampErr) + " " + str(wid) + " " + str(widErr) + " " + str(area) + " " + str(areaErr))
        text_file.close()
    '''
