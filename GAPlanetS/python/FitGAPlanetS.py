#program searches through the GAPlanetS_working directory on GPIcruncher to find image cubes that have been run through visao_circlesym.pro and fits star and ghost radial brightness profiles from these images to Moffat profiles. Then runs statistics on these fits. 


##############################################################################
#                             IMPORT LIBRARIES                               #
##############################################################################
                
import os
import math
import matplotlib
matplotlib.use('nbagg')  # required for interactive plottin
from astropy.io import fits 
import matplotlib.pyplot as plt
import numpy as np
import lmfit
from scipy.optimize import curve_fit
from numpy import sqrt, pi, exp, linspace
from lmfit import Model
import os.path
get_ipython().magic('pylab')
get_ipython().magic('matplotlib inline')
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
from __future__ import print_function, division
from astropy.table import Table
from matplotlib import style
style.use('ggplot')  
from os.path import expanduser
import glob
import pandas as pd


#switches directory, parameter is directory path 
def switchDirectories(directory):
    os.chdir(directory)
    
#move to home directory
def goHome():
    home = expanduser('~')
    while (not os.getcwd() == home):
        switchDirectories('..')
        
#Visually organize program flow
def printDivider():
    print("--------------------------------------------------------------------------------")
    
#reads in fits file
def readFits(filename):
    hdulist = fits.open(filename)
    data = hdulist[0].data
    return data

#Create 3D array from fits file cube
def readFitsCubeToArray(filename, filetype):
    global dataHA
    global dataCONT
    print("Now reading in data cube \'"+str(filename)+"\'")
    hdulist = fits.open(filename)
    return hdulist[0].data
    #NOTE: Array is in format data[Z][Y][X]
    
#creates median image from data cube
def median(dataCube, filetype):
    if filetype == 1:
        print("Taking median of H-alpha data")
    elif filetype == 2:
        print("taking median of continuous spectrum data")
    medianImage = np.nanmedian(dataCube, axis=0)
    print("Median completed; shape of new data is " + str(medianImage.shape))
    return medianImage

#creates array from center row of pixels of a 450x450 image
def isolateCentralArray(data, filetype):
     center = []
     for x in range (170,280):
         value1 = data[x][224]
         value2 = data[x][225]
         value = (value1 + value2)/2
         center.append(value)
     return center

#locates the center of a ghost image by looking around an approximated center for the brightest pixel in a 10 pixel radius
def findGhostSingleIm(data, filetype):
    if filetype == 1:
        print("Scanning Hydrogen-Alpha image to evaluate center of stellar ghost")
    elif filetype == 2:
        print("Scanning images in continuous spectrum image to evaluate center of stellar ghost")
    count = 0
    #xGhost = int(input("Enter the x coordinate of an approximated center of the ghost: "))
    #yGhost = int(input("Enter the y coordinate of an approximated center of the ghost: "))
    xGhost = 383
    yGhost = 217
    print("Locating the brightest pixel of the ghost")
    maxVal = 0
    tempX = 0
    tempY = 0
    for y in range(yGhost-10, yGhost+10):
        for x in range(xGhost-10,xGhost+10):
            if data[y,x] > maxVal:
                tempX = x
                tempY = y
                maxVal = data[y,x]
    print("Calculated ghost center is " + str(maxVal) + " at x=" + str(tempX) + " , y=" + str(tempY))
    print(str(filetype))
    return (tempX, tempY, maxVal)

#creates 50 pix array from center of ghost
def isolateGhostArray(data, filetype):
     center = []
     if (filetype == 1):
         for x in range (GhostxHA-25, GhostxHA+25):
             value = data[GhostyHA,x]
             center.append(value)
     if(filetype == 2):
         for x in range (GhostxCONT-25, GhostxCONT+25):
             value = data[GhostyCONT,x]
             center.append(value)
     return center

#modifies data sets by removing non-linear values. An independent variable is then added, an array containing all numbers from 0 to the length of the y variable. Values are removed at the same indeces as the removed y values.
#determines whether the data set is saturated
def modifyDataSet(y):
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
    
#mfoffat function definition, takes a data set as well as an approximated center, width, fwhm, and exponent as parameters
def moffat(x, amp, cen, wid, pow):
    return (amp*(((x-cen)/wid)**(2)+1)**(-pow))

#fits a moffat profile to ghost array found by isolateGhostArray funciton, takes this array as a parameter
def makeMoffatGhost(aghostCenterHA):
    try:
        #an array of x values from 0 to 50
        new_x = list(range(len(aghostCenterHA)))

        #uses the ghost array as y values
        new_y = new_y = aghostCenterHA
                
        #removes linear star halo
        haloSlope = -1*(new_y[0]-new_y[len(new_x)-1])/(len(new_x)-1)
        haloAdj = -1*(new_y[0]-new_y[len(new_x)-1])
        print('the approximate slope of the linear star halo across the ghost is: %s' %haloAdj)
        
        if(haloAdj < -1):
            for x in range (len(new_x)):
                new_y[x]= new_y[x] - new_x[x]*haloSlope    
        
        # initial (dumb) guess for [amp, cen, wid] of gaussian
        init_vals = [100, 25, 3, 2]
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

        yres = new_y[0] 
        for x in range (len(new_y)):
            new_y[x]= new_y[x]-yres 
        
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
        xfine = np.linspace(0., 50., 30000)  # define values to plot the fit for
        #plt.plot(xfine, gaussian(xfine, best_vals2[0], best_vals2[1], best_vals2[2]), 'r-')
        plt.plot(xfine, moffat(xfine, *best_vals2), 'c-')
        plt.xlabel('Pixel Number')
        plt.ylabel('Pixel Value')
        plt.legend(('moffat best fit', 'data'), loc='upper left')
        plt.title("%s %s ghost" %(name, date))
        print(str(best_vals2))
        
        if (len(extra) == 0):
            plt.savefig('%s_%s_ghost' %(name,date))
        else:
            plt.savefig('%s_%s_%s_ghost' %(name,date, extra))
        
        goHome()
        if (len(extra) == 0):
            plt.savefig('../GAPlanetS/GAPlanetS_Working/curveFits/%s_%s_ghost' %(name,date))
        else:
            plt.savefig('../GAPlanetS/GAPlanetS_Working/curveFits/%s_%s_%s_ghost' %(name,date, extra))
        plt.clf()

        def moff(x):
            return moffat(x, *best_vals2)
        area, areaErr = quad(moff, 0, 119)
        amp = best_vals2[0]
        ampErr = covar2[0,0]**.5
        print("the area under the curve is", str(area))
        wid = best_vals2[2]
        widErr = covar2[2,2]**.5
        ghostPeaks.append(amp)
        if (amp>0):
            return amp, ampErr, wid, widErr, area, areaErr
        else:
            return 'no fit', 'no fit','no fit','no fit','no fit','no fit'

    except (RuntimeError):
            print ('signal/noise too low; no optimal fit found')
            
                #an array of x values from 0 to 50
            new_x = list(range(len(aghostCenterHA)))

            #uses the ghost array as y values
            new_y = new_y = aghostCenterHA

            #removes linear star halo
            haloSlope = -1*(new_y[0]-new_y[len(new_x)-1])/(len(new_x)-1)
            print('the approximate slope of the linear star halo across the ghost is: %s' %haloSlope)

            if(haloSlope < -1):
                for x in range (len(new_x)):
                    new_y[x]= new_y[x] - new_x[x]*haloSlope

            yres = new_y[0] 
            for x in range (len(new_y)):
                new_y[x]= new_y[x]-yres     

            #add noise
            #read noise
            rd_noise = np.repeat(0.02, len(new_x))

            #photon noise
            #just so that there are no negative values, add a constant. this is an oversimplification, but will work for now.
            floor=np.repeat(abs(min(new_y)),len(new_x)) 
            pt_noise = sqrt(new_y+floor) #add me - /sqrt(nimages)

            #noise adds in quadrature
            e = sqrt(rd_noise**2 + pt_noise**2)
            
            plt.figure(2)
            xfine = np.linspace(0., 120., 30000)  # define values to plot the fit for
            #plt.plot(xfine, gaussian(xfine, best_vals2[0], best_vals2[1], best_vals2[2]), 'r-')
            plt.errorbar(new_x, new_y, yerr=e, fmt='none', ecolor = 'r')
            plt.xlabel('Pixel Number')
            plt.ylabel('Pixel Value')
            plt.legend(('moffat best fit', 'data'), loc='upper left')
            plt.title("No ghost Fit for %s %s" %(name, date))
            
            if (len(extra) == 0):
                plt.savefig('%s_%s_ghost' %(name,date))
            else:
                plt.savefig('%s_%s_%s_ghost' %(name,date, extra))
        
            goHome()
            if (len(extra) == 0):
                plt.savefig('../GAPlanetS/GAPlanetS_Working/curveFits/%s_%s_ghost' %(name,date))
            else:
                plt.savefig('../GAPlanetS/GAPlanetS_Working/curveFits/%s_%s_%s_ghost' %(name,date, extra))

            plt.clf()

            return 'no fit', 'no fit','no fit','no fit','no fit','no fit'
        
#fits moffat profile to array of pixels accross the center of the star
def makeMoffatStar(acenterMedHA, aXValsHA):
    try:
        new_x = aXValsHA

        new_y = new_y = acenterMedHA

        # initial (dumb) guess for [amp, cen, wid] of gaussian
        init_vals = [100000, 25, 6, .5]
        #curvefit spits out two arrays - one with the fit parameters (best_vals) and one with the covariance matrix (covar1)
        best_vals1, covar1 = curve_fit(moffat, new_x, new_y, p0=init_vals)
        print(best_vals1)
        #print(covar)
        #add noise

        #read noise
        rd_noise = np.repeat(0.02, len(new_x))

        #photon noise
        #just so that there are no negative values, add a constant. this is an oversimplification, but will work for now.
        floor=np.repeat(abs(min(new_y)),len(new_x)) 
        pt_noise = sqrt(new_y+floor) #add me - /sqrt(nimages)
        #this noise, as you will see below, seems way too big given how clean your median radial profile will be. This is because 
        #when you combine images, photon noise is reduced by the square root of the number of exposures that went into your image
        
        yres = new_y[0] 
        for x in range (len(new_y)):
            new_y[x]= new_y[x]-yres
        
        #noise adds in quadrature
        e = sqrt(rd_noise**2 + pt_noise**2)
        #plt.errorbar(new_x, new_y, yerr=e, fmt='none')
        best_vals2, covar2 = curve_fit(moffat, new_x, new_y, p0=init_vals, sigma=e, maxfev=10000)
        #best_vals2
        #z = covar1 - covar2
        #z
        ###print("a =", popt[0], "+/-", pcov[0,0]**0.5)
        ###print ("b =", popt[1], "+/-", pcov[1,1]**0.5)


        print("amp =", best_vals2[0], "+/-", covar2[0,0]**0.5)
        print ("cen =", best_vals2[1], "+/-", covar2[1,1]**0.5)
        print ("wid =", best_vals2[2], "+/-", covar2[2,2]**0.5)
        plt.figure(2)
        plt.errorbar(new_x, new_y, yerr=e, fmt='none', ecolor = 'r')

        #np.linspace?
        xfine = np.linspace(0., 120., 30000)  # define values to plot the fit for
        #plt.plot(xfine, gaussian(xfine, best_vals2[0], best_vals2[1], best_vals2[2]), 'r-')
        plt.plot(xfine, moffat(xfine, *best_vals2), 'c-')
        plt.xlabel('Pixel Number')
        plt.ylabel('Pixel Value')
        plt.legend(('moffat best fit', 'data'), loc='upper left')
        plt.title("%s %s star" %(name, date))
        print(str(best_vals2))
        
        if (len(extra) == 0):
            plt.savefig('%s_%s_star' %(name,date))
        else:
            plt.savefig('%s_%s_%s_star' %(name,date, extra))
        
        goHome()
        if (len(extra) == 0):
            plt.savefig('../GAPlanetS/GAPlanetS_Working/curveFits/%s_%s_star' %(name,date))
        else:
            plt.savefig('../GAPlanetS/GAPlanetS_Working/curveFits/%s_%s_%s_star' %(name,date, extra))

        plt.clf()

        def moff(x):
            return moffat(x, *best_vals2)
        area, areaErr = quad(moff, 0, 119)
        amp = best_vals2[0]
        ampErr = covar2[0,0]**.5
        print("the area under the curve is", str(area))
        wid = best_vals2[2]
        widErr = covar2[2,2]**.5

        return amp, ampErr, wid, widErr, area, areaErr
   
    except (RuntimeError):
        print ('signal/noise too low; no optimal fit found')
        plt.figure(2)
        xfine = np.linspace(0., 120., 30000)  # define values to plot the fit for
        plt.xlabel('Pixel Number')
        plt.ylabel('Pixel Value')
        plt.legend(('moffat best fit', 'data'), loc='upper left')
        plt.title("No Star Fit for %s %s" %(name, date))
        
        if (len(extra) == 0):
            plt.savefig('%s_%s_star' %(name,date))
        else:
            plt.savefig('%s_%s_%s_star' %(name,date, extra))
        
        goHome()
        if (len(extra) == 0):
            plt.savefig('../GAPlanetS/GAPlanetS_Working/curveFits/%s_%s_star' %(name,date))
        else:
            plt.savefig('../GAPlanetS/GAPlanetS_Working/curveFits/%s_%s_%s_star' %(name,date, extra))

        plt.clf()
        
        return 'no fit', 'no fit','no fit','no fit','no fit','no fit'
    
#looks at directory path to find star name, and date a given set of images was taken
#function assumes directory path format: '~/../GAPlanetS/GAPlanetS_Working/star_name/observation_date/setinfo'
def getName():
    #creates string from directory path
    directory = str(os.getcwd())
    
    #finds index of "Desktop" in directory string
    start = directory.find('GAPlanetS_Working') +10
    
    #array to hold star name
    name = []
    #array to hold observation date
    date = []
    
    stop = 0
    
    #array to hold set information 
    extra = []
    
    #loops through the directory path starting after 'GAPlanetS_Working' adding each following character to 'name' array; loop breaks once '/' is found
    for x in range (start, len(directory)):
        if (directory[x] == '/'):
            #stop value shows where to start looking for observation date
            stop = x
            break
        else:
            name.append(directory[x])
    stop2=0

    #looks for date if the name directory contains subdirectories
    if (stop >0):
        #loops through the directory path starting after name directory adding each following character to 'date' array; loop breaks once '/' is found
        for x in range (stop+1, len(directory)):
            if (directory[x] == '/'):
                stop2 = x
                break
            else:
                date.append(directory[x])
                
    #looks for set info if the date directory contains subdirectories
    if (stop2 >0):
        #loops through the directory path starting after name directory adding each following character to 'date' array; loop breaks once '/' is found
        for x in range (stop2+1, len(directory)):
            if (directory[x] == '/'):
                stop = x
                break
            else:
                extra.append(directory[x])
                            
    #changes 'name' and 'date' arrays to strings
    name = ''.join(name)
    date = ''.join(date)
    extra = ''.join(extra)
    #prints name of star and observation date
    print(name)
    print(date)
    if (len(extra)>0):
        print(extra)
    #returns name and date values
    return name, date, extra

#calcuulated mean and standard deviation, while also propagating errors for each
def calcStats(data, errors):
    acc = 0
    meanErrSq = 0
    for x in range (len(data)-1):
        acc = acc + data[x]
        meanErrSq = meanErrSq + errors[x]**2
    mean = acc / (len(data)-1) 
    mean = mean[0]
    meanErr = meanErrSq**.5/(len(data)-1)
    meanErr = meanErr[0]
    print('the mean of the data is :')
    print('%s +\- %s' %('{:0.4e}'.format(mean), '{:0.4e}'.format(meanErr)))
    DevSum = 0
    DevSumErr= 0
    for x in range (len(data)-1):
        dev = (mean-data[x])
        devSq=(mean-data[x])**2
        DevSum = DevSum + devSq
        devErr = (meanErr**2 + errors[x]**2)**.5
        devSqErr = 2*(devErr/dev)*devSq
        DevSumErr = DevSumErr + devSqErr**2
    DevSumErr = DevSumErr**.5
    varianceErr = DevSumErr/(len(data)-2)
    variance = DevSum/(len(data)-2)
    stDev= variance**.5
    errSTD = .5*(varianceErr/variance)*stDev
    stDev = stDev[0]
    errSTD = errSTD[0]
    print('the standard deviation of the data is :')
    print('%s +\- %s' %('{:0.4e}'.format(stDev), '{:0.4e}'.format(errSTD)))
    return mean, meanErr, stDev, errSTD
        
#removes elements from an array consisting of string 'no fit'
def removeBadFits(data):
    #creates array sshowing which indeces do not contain "no fit" strings
    idxs = numpy.any(data != 'no fit', axis=1)
    #creates array conaining only elememts from data at given indeces
    data = data[idxs, :]
    #returns modified array
    return data

###############################################
##                  main                     ##
###############################################

#change directories to start in home directory
goHome()

#creates an array containing a list of all possible directory paths extending from the GAPlantS directory
directoryList = [x[0] for x in os.walk('../GAPlanetS/GAPlanetS_Working')]


#creates a directory called curveFits to hold graphs if it does not already exist 
if not os.path.exists('../GAPlanetS/GAPlanetS_Working/curveFits'):
    os.makedirs('../GAPlanetS/GAPlanetS_Working/curveFits')

#creates a directory called tableData to hold data frames if it does not already exist 
if not os.path.exists('../GAPlanetS/GAPlanetS_Working/tableData'):
    os.makedirs('../GAPlanetS/GAPlanetS_Working/tableData')

#creates empty pandas data frame with appropriate columns for unsaturated data sets 
tdUnsat = pd. DataFrame(columns = ('name; date; set', 'saturated?', 'star peak', 'star peak error', 'ghost peak', 'ghost peak error', 'star fwhm', 'star fwhm error', 'ghost fwhm', 'ghost fwhm error', 'total starlight', 'total starlight error', 'total light from ghost', 'light from ghost error', 'ghost/star peak ratio', 'peak ratio error', 'ghost/star total light ratio', 'total light ratio error'))
#sets index to the name of the star 
tdUnsat.set_index(['name; date; set'], inplace= True)

#creates empty pandas data frame with appropriate columns for saturated data sets
tdSat = pd. DataFrame(columns = ('name; date; set', 'saturated?', 'star peak', 'star peak error', 'ghost peak', 'ghost peak error', 'star fwhm', 'star fwhm error', 'ghost fwhm', 'ghost fwhm error', 'total starlight', 'total starlight error', 'total light from ghost', 'light from ghost error', 'ghost/star peak ratio', 'peak ratio error', 'ghost/star total light ratio', 'total light ratio error'))
tdSat.set_index('name; date; set', inplace= True)

print('CALCULATING BEST FITS:')
print('')
#for loop loops through all possible directory paths 
for x in range (len(directoryList)-1):
    #changes to next directory
    switchDirectories(directoryList[x])
    ##creates string from directory path
    #directory = str(directoryList[x])
    #directory = directory.replace('/','_')
   
    #creates array of all files ending in "circsym.fits"
    circsymFiles = glob.glob('*circsym.fits')
    
    #continues if directory contains any files ending in "circsym.fits"
    if (len(circsymFiles) > 0):
        #prints directory path
        print(os.getcwd())
        
        #saves name of star and observation date to variables
        name, date, extra = getName()
        
        #arrays to hold peaks of various stars
        starPeaks = []
        ghostPeaks = []

        #stores data cubes that have been run through 'circlesym.pro' as arrays
        dataHA = readFitsCubeToArray("Line_clip450_reg_circsym.fits", 1) #Data is stored in dataHA                                                                                            
        dataCONT = readFitsCubeToArray("Cont_clip450_reg_circsym.fits", 2) #Data is stored in dataCONT 

        #creates median image from data cubes if files do not already exist                                                                                                                                        
        if (not os.path.isfile('medianHA.fits')):
            dataMedHA = median(dataHA, 1) #takes median of HA cube                                                                                                                   
            hdu = fits.PrimaryHDU(dataMedHA)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto('medianHA.fits')
        else:
            dataMedHA = readFits('medianHA.fits')

        if (not os.path.isfile('medianCONT.fits')):
            dataMedCONT = median(dataCONT, 2) #takes median of CONT cube                                                                                                             
            hdu = fits.PrimaryHDU(dataMedCONT)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto('medianCONT.fits')
        else:
            dataMedCONT = readFits('medianCONT.fits')

        #arrays to hold arrays showing centers of stars                                
        centerMedHA = isolateCentralArray(dataMedHA, 1)
        centerMedCONT = isolateCentralArray(dataMedCONT,2)

        #converts non-linear values in saturated data to NaNs
        #determines whether image is saturated
        XValsHA, centerMedHA, saturationS = modifyDataSet(centerMedHA)
        XValsCONT, centerMedCONT, saturationG = modifyDataSet(centerMedCONT)

        #finds the center of the 'ghost' in HA image
        GhostxHA, GhostyHA, HAghostPeak = findGhostSingleIm(dataMedHA, 1)
        #makes array accross the center of ghost in HA image
        ghostCenterHA = isolateGhostArray(dataMedHA, 1)
        printDivider()
        print('ghost fits:')
        #fits the ghost to a moffat profile, extracts the aplitude, fwhm, and area under the curve as well as the errors in these measurements, and saves graphs to current directory as png files
        gAmp, gAmpErr, gWid, gWidErr, gArea, gAreaErr = makeMoffatGhost(ghostCenterHA)
        printDivider()
        print('star fits:')
        #fits the ghost to a moffat profile, extracts the aplitude, fwhm, and area under the curve as well as the errors in these measurements, and saves graphs to current directory as png file
        sAmp, sAmpErr, sWid, sWidErr, sArea, sAreaErr = makeMoffatStar(centerMedHA, XValsHA)
        printDivider()
        printDivider()
        #moves back to home directory
        goHome()
        
        try:
            ampRatio = gAmp/sAmp
            areaRatio = gArea/sArea
            ampRatErr = ((gAmpErr/gAmp)**2 + (sAmpErr/sAmp)**2)**.5
            areaRatErr = ((gAreaErr/gArea)**2 + (sAreaErr/sArea)**2)**.5
        except(TypeError):
            ampRatio = 'no fit'
            areaRatio = 'no fit'
            ampRatErr = 'no fit'
            areaRatErr = 'no fit'
            
        #adds fit statistics to data frame
        if (saturationS ==  'unsat'):
            tdUnsat.loc['%s; %s; %s' %(name, date, extra)] = (saturationS, sAmp,sAmpErr,gAmp,gAmpErr,sWid, sWidErr, gWid,gWidErr,sArea,sAreaErr,gArea,gAreaErr, ampRatio, ampRatErr, areaRatio, areaRatErr)
       
        if (saturationS ==  'sat'):
            tdSat.loc['%s; %s; %s' %(name, date, extra)] = (saturationS, sAmp,sAmpErr,gAmp,gAmpErr,sWid, sWidErr, gWid,gWidErr,sArea,sAreaErr,gArea,gAreaErr, ampRatio, ampRatErr, areaRatio, areaRatErr)
       
    #moves back to home directory
    goHome()    

#organizes output
printDivider()
printDivider()
tdUnsat = tdUnsat.round(4)

#saves sat data frame to html 
tdSat.to_html('../GAPlanetS/GAPlanetS_Working/tableData/tableSatData.html')

#saves unsat data frame to html 
tdUnsat.to_html('D../GAPlanetS/GAPlanetS_Working/tableData/tableUnsatData.html')

#shows data frame of unsaturated data
tdUnsat.round(4)

#shows data frame with saturated data
tdSat

#organizes output
printDivider()
printDivider()
print('RUNNING STATS:')
print('')

#creates array from unsat data frame column containing ratios of ghost to star total light
lightRatio = tdUnsat.as_matrix(columns=tdUnsat.columns[15:16])

#creates array from unsat data frame column containing errors for above calculations
lightRatioErr = tdUnsat.as_matrix(columns=tdUnsat.columns[16:17])

#removes data containing bad fits
lightRatio = (removeBadFits(lightRatio))
lightRatioErr = (removeBadFits(lightRatioErr))

#prints total light ratios
print('The ratios between total light from each ghost to that of its respective star are as follows: ')
print(", ".join(map(str, ('%.4f' % elem for elem in lightRatio))))

#finds standard deviation and mean for total light ratio data, also finds errors for each calculation
lightRatioAv, lightRatioAvErr, lightRatioSTD, lightRatioSTDErr = calcStats(lightRatio, lightRatioErr)

#organizes output
printDivider()
print('')

#creates array from unsat data frame column containing ratios of ghost to star peak
ampRatio = tdUnsat.as_matrix(columns=tdUnsat.columns[13:14])

#creates array from unsat data frame column containing errors for above calculations
ampRatioErr = tdUnsat.as_matrix(columns=tdUnsat.columns[14:15])

#removes data containing bad fits
ampRatio = (removeBadFits(ampRatio))
ampRatioErr = (removeBadFits(ampRatioErr))

#prints peak ratios
print('The ratios between peak of each ghost to that of its respective star are as follows: ')
print(", ".join(map(str, ('%.4f' % elem for elem in ampRatio))))

#finds standard deviation and mean for amp ratio data, also finds errors for each calculation
ampRatioAv, ampRatioAvErr, ampRatioSTD, ampRatioSTDErr = calcStats(ampRatio, ampRatioErr)

#organizes output
printDivider()
print('')

#creates empty pandas data frame for statistics
tdStats= pd. DataFrame(columns = ('ratio', 'mean ratio', 'mean error', 'standard deviation', 'stdev error'))
#sets index to the ratio title 
tdStats.set_index('ratio', inplace= True)

#adds light ratio stats to data frame
tdStats.loc['total light'] = (lightRatioAv, lightRatioAvErr, lightRatioSTD, lightRatioSTDErr)

#adds peak ratio stats to data frame
tdStats.loc['peak value'] = (ampRatioAv, ampRatioAvErr, ampRatioSTD, ampRatioSTDErr)

#saves stats data frame to html 
tdStats.to_html('../GAPlanetS/GAPlanetS_Working/tableData/tableUnsatStats.html')

#shows stats data frame
tdStats
