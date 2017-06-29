import os
import sys
import math
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit
from scipy.integrate import quad
import pandas as pd
from numpy import sqrt

def readFitsCubeToArray(filename):
    global data
    hdulist = fits.open(filename)
    data = hdulist[0].data
    hdulist.close()
    print("read fits file into data")
    medianImage = np.nanmedian(data, axis=0)
    data = medianImage
    print("created median")

def findGhostSingleIm():
    global ghostX
    global ghostY
    #count = 0
    xGhost = 383
    yGhost = 217
    maxVal = 0
    tempX = 0
    tempY = 0
    for y in range(yGhost-15, yGhost+16):
        for x in range(xGhost-15,xGhost+16):
            if data[y,x] > maxVal:
                tempX = x
                tempY = y
                maxVal = data[y,x]
    print("Calculated ghost center is " + str(maxVal) + " at x=" + str(tempX) + " , y=" + str(tempY))
    ghostX = tempX
    ghostY = tempY
    return (tempX, tempY, maxVal)

def makeSquare():
    global maskedArray
    maskedArray = np.zeros((31,31))
    for y in range(-15,16):
        for x in range(-15,16):
            maskedArray[y+15][x+15] = data[ghostY+y][ghostX+x]

def isolateGhostArray():
     center = []
     for x in range (0, 30):
         value = maskedArray[15,x]
         center.append(value)
     print(center)
     return center

def modifyDataSet():
    y = data
    length = len(y)
    x = list(range(len(y)))
    i = 0 
    while i < (len(x)):
        if (y[i]>15500 or x[i]== np.nan):
            x.remove(x[i])
            y.remove(y[i])
        else:
            i = i+1
    if (len(y) == length):
        return x,y, 'unsat'
    else:
        return x, y, 'sat'

def moffat(x, amp, cen, wid, pow):
    return (amp*(((x-cen)/wid)**(2)+1)**(-pow))

def makeMoffatGhost(array):
    try:
        #an array of x values corresponding to the size of the array
        new_x = list(range(len(array)))

        #uses the ghost array as y values
        new_y = array

        # initial (dumb) guess for [amp, cen, wid] of moffat
        init_vals = [1, 15, 10, 0.5]
        #curvefit spits out two arrays - one with the fit parameters (best_vals) and one with the covariance matrix (covar1)
        best_vals1, covar1 = curve_fit(moffat, new_x, new_y, p0=init_vals)
        print(best_vals1)
        
        #add noise
        #read noise
        rd_noise = np.repeat(0.02, len(new_x))

        #photon noise
        #just so that there are no negative values, add a constant. this is an oversimplification, but will work for now.
        floor=np.repeat(abs(min(new_y)),len(new_x)) 
        pt_noise = sqrt(new_y+floor) #add me - /sqrt(nimages)
        #noise adds in quadrature
        e = sqrt(rd_noise**2 + pt_noise**2)
        #plt.errorbar(new_x, new_y, yerr=e, fmt='none')
        best_vals2, covar2 = curve_fit(moffat, new_x, new_y, p0=init_vals, sigma=e)
        #best_vals2


        print("amp =", best_vals2[0], "+/-", covar2[0,0]**0.5)
        print ("cen =", best_vals2[1], "+/-", covar2[1,1]**0.5)
        print ("wid =", best_vals2[2], "+/-", covar2[2,2]**0.5)
        plt.figure(1)
        plt.errorbar(new_x, new_y, yerr=e, fmt='none', ecolor = 'r')
        #np.linspace?
        xfine = np.linspace(0., 30., 30000)  # define values to plot the fit for
        #plt.plot(xfine, gaussian(xfine, best_vals2[0], best_vals2[1], best_vals2[2]), 'r-')
        plt.plot(xfine, moffat(xfine, *best_vals2), 'c-')
        plt.xlabel('Pixel Number')
        plt.ylabel('Pixel Value')
        plt.legend(('moffat best fit', 'data'), loc='upper left')
        #plt.title("%s %s ghost" %(name, date))
        print(str(best_vals2))
        #plt.savefig('%s_%s_ghost.png' %(name,date))
        #plt.clf()

        def moff(x):
            return moffat(x, *best_vals2)
        
        area, areaErr = quad(moff, 0, 119)
        amp = best_vals2[0]
        ampErr = covar2[0,0]
        print("the area under the curve is", str(area))
        wid = best_vals2[2]
        widErr = covar2[2,2]
        #ghostPeaks.append(amp)
        return amp, ampErr, wid, widErr, area, areaErr

    except (RuntimeError):
            print ('signal/noise too low; no optimal fit found')
            plt.figure(2)
            xfine = np.linspace(0., 120., 30000)  # define values to plot the fit for
            #plt.plot(xfine, gaussian(xfine, best_vals2[0], best_vals2[1], best_vals2[2]), 'r-')
            plt.xlabel('Pixel Number')
            plt.ylabel('Pixel Value')
            plt.legend(('moffat best fit', 'data'), loc='upper left')
            #plt.title("No Star Fit for %s %s" %(name, date))
            #plt.savefig('%s_%s_ghost.png' %(name,date))

            plt.clf()

            return 'no fit', 'no fit','no fit','no fit','no fit','no fit'

def normalizeSquare():
    maxVal = 0
    for y in range(31):
        for x in range(31):
            if maskedArray[y][x] > maxVal:
                maxVal = maskedArray[y][x]
    for y in range(31):
        for x in range(31):
            maskedArray[y][x] = maskedArray[y][x] / maxVal


readFitsCubeToArray("Line_clip450_flat_reg_circsym_nocosmics.fits")
findGhostSingleIm()
makeSquare()
normalizeSquare()
#modifyDataSet()
amp, ampErr, wid, widErr, area, areaErr = makeMoffatGhost(isolateGhostArray())
hdu = fits.PrimaryHDU(maskedArray)
hduList = fits.HDUList([hdu])
hduList.writeto("test.fits", overwrite=True)
with open("moffat.txt","w") as text_file:
    text_file.write(str(amp) + " " + str(ampErr) + " " + str(wid) + " " + str(widErr) + " " + str(area) + " " + str(areaErr))
    text_file.close()
