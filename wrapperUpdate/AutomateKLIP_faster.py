#Clare Leonard                                                                
#Version 1.0 - 6/21/17                                                          




##################################################################
#############                                        #############
#############                 IMPORTS                #############
#############                                        #############
##################################################################

import glob
import inspect
import os                                      
import instruments.MagAO as MagAO                                   
import parallelized as parallelized
import numpy as np
import sys             
import pyklip.fm as fm
import klip as klip
from astropy.io import fits
import SNRMap as snr   
import time
import warnings
from astropy.utils.exceptions import AstropyWarning

warnings.filterwarnings('ignore', category=AstropyWarning, append=True)


##################################################################
#############                                        #############
#############               SAVE FILES               #############
#############                                        #############
################################################################## 
 

def writeData(indiv, allParams = False, snrmap = False, pre = ''): 
    #function writes out fits files and writes important information to fits headers
    
    hdu = fits.PrimaryHDU(indiv)
    hdulist = fits.HDUList([hdu])
 
    
    if (allParams):
    #creates new strings to add parameter information to file names
        annuli_fname = annuli
        annuli_head = annuli
        movement_fname = movement
        movement_head = movement
        subsections_fname = subsections
        subsections_head = subsections
       
        #if program iterates over several parameter values, formats these for fits headers and file names
        if (isinstance(annuli, tuple)):
            annuli_fname = str(annuli[0]) + '-' + str(annuli[1]) + 'x' + str(annuli[2])
            annuli_head = str(annuli[0]) + 'to' + str(annuli[1]) + 'by' + str(annuli[2])  
        if (isinstance(movement, tuple)):
            movement_fname = str(movement[0]) + '-' + str(movement[1]) + 'x' + str(movement[2])
            movement_head = str(movement[0]) + 'to' + str(movement[1]) + 'by' + str(movement[2])
        if (isinstance(subsections, tuple)):
            subsections_head = str(subsections[0]) + 'to' + str(subsections[1]) + 'by' + str(subsections[2])
            subsections_fname = str(subsections[0]) + '-' + str(subsections[1]) + '-' + str(subsections[2])
    else:
        annuli_head = a
        movement_head = m
        subsections_head = s
        annuli_fname = a
        movement_fname = m
        subsections_fname = s

    #shortens file path to bottom 4 directories so it will fit in fits header
    try:
        pathToFiles_short = '/'.join(pathToFiles.split(os.path.sep)[-4:])
    except:
        pathToFiles_short = pathToFiles
            
    #adds info to fits headers
    prihdr.set('annuli', str(annuli_head))
    prihdr.set('movement', str(movement_head))
    prihdr.set('subsctns', str(subsections_head))
    prihdr.set('klmodes', str(klmodes))
    prihdr.set('filepath', str(pathToFiles_short))
 
    if(snrmap):
        rad, pa, wid = mask 
        prihdr.set('mask_rad', str(rad))
        prihdr.set('mask_pa', str(pa))
        prihdr.set('mask_wid', str(wid))
  
        prihdr.set('smooth_val', str(_smooth))

        prihdr.set('FWHM', str(FWHM))
   
    hdulist[0].header = prihdr
    
    #writes out files
    hdulist.writeto(str(pathToFiles) + "_klip/" + str(pre)  + outputFileName + "_a" + str(annuli_fname) + "m" + str(movement_fname) + "s" + str(subsections_fname) + "iwa" + str(iwa) + '_klmodes-all.fits', clobber=True)




##################################################################
#############                                        #############
#############               GET INPUTS               #############
#############                                        #############
################################################################## 

#value adjusts argument numbering in case of white space in file path 
argnum = 0


pathToFiles = sys.argv[1]

#if the file path has white space in it, recognizes the end of the filepath by the phrase '%finish'
#If the phrase '%finish' does not occur, leaves pathToFiles as the first argument
try:
    while (not pathToFiles[-7:] == '%finish'):
        argnum += 1
        pathToFiles = pathToFiles + " " + sys.argv[1+argnum]
    pathToFiles = pathToFiles[:-7]
        
except:
    pathToFiles = sys.argv[1]
    argnum = 0
    
print("File Path = " + pathToFiles)   
print()

#create directory to save ouput to
if not os.path.exists(pathToFiles + "_klip"):
    os.makedirs(pathToFiles + "_klip")
    os.chmod(pathToFiles + "_klip", 0o777)

print("Parameters to explore:")

annuli_start = int(sys.argv[5+argnum])
annuli_stop = int(sys.argv[6+argnum])
annuli_inc = int(sys.argv[7+argnum])
#creates touple for easier eventual string formatting when saving files
annuli = (annuli_start, annuli_stop, annuli_inc)
#if only one parameter is iterated over, makes sure increment is 1 and changes touple to single int
if(annuli_start == annuli_stop):
    annuli_inc = 1;
    annuli = annuli_start
    print("Annuli = " +str(annuli_start))    
else:
    print("Annuli: start = %s; end = %s; increment = %s " %(str(annuli_start), str(annuli_stop), str(annuli_inc)))



movement_start = float(sys.argv[8+argnum])
movement_stop = float(sys.argv[9+argnum])
movement_inc = float(sys.argv[10+argnum])
#creates touple for easier eventual string formatting when saving files
movement = (movement_start, movement_stop, movement_inc)
#if parameter is not set to change, makes sure increment is 1 and changes touple to single int
if(movement_start == movement_stop):
    movement_inc = 1;
    movement = movement_start
    print("Movement = " +str(movement_start))
    
else:
    print("Movement: start = %s; end = %s; increment = %s " %(str(movement_start), str(movement_stop), str(movement_inc)))


subsections_start = int(sys.argv[11+argnum])
subsections_stop = int(sys.argv[12+argnum])
subsections_inc = int(sys.argv[13+argnum])
#creates touple for easier eventual string formatting when saving files
subsections = (subsections_start, subsections_stop, subsections_inc)
#if parameter is not set to change, makes sure increment is 1 and changes touple to single int
if(subsections_start == subsections_stop):
    subsections_inc = 1;
    subsections = subsections_start
    print("Subsections = " +str(subsections_start))    
else:
    print("Subsections: start = %s; end = %s; increment = %s " %(str(subsections_start), str(subsections_stop), str(subsections_inc)))

print()
    
iwa = int(sys.argv[2+argnum])
print("IWA = " + str(iwa))

print("KL Modes = " + str(list(map(int, sys.argv[3+argnum].split(",")))))
klmodes = list(map(int, sys.argv[3+argnum].split(",")))

FWHM = float(sys.argv[14+argnum])
print("Star FWHM = " + str(FWHM))

_smooth = float(sys.argv[15+argnum])
print("Smoothing Value = " + str(_smooth))

print()
print('Planet mask parameters:')

ra = list(map(int, sys.argv[16+argnum].split(",")))
print("Radius = " + str(ra))

pa = list(map(int, sys.argv[17+argnum].split(",")))
print("Position Angle = " + str(pa))

wid = list(map(int, sys.argv[18+argnum].split(",")))
print("Mask width (radial, angular): = " + str(wid))

#object to hold mask parameters for snr map 
mask = (ra, pa, wid)
                
print()

outputFileName = sys.argv[4+argnum]
print("Output FileName = " + outputFileName)


saveSNR = False
if (sys.argv[19+argnum] == 'true' or sys.argv[19+argnum] == 'True'):
    saveSNR = True    

    
singleAnn = False
if (sys.argv[20+argnum] == 'true' or sys.argv[20+argnum] == 'True'):
    singleAnn = True    
    
print()





##################################################################
#############                                        #############
#############            PERFORM AUTOMATION          #############
#############                                        #############
################################################################## 

print("reading: " + pathToFiles + "/*.fits")

print()

#grab header
hdulist = fits.open(pathToFiles + '/sliced_1.fits')
prihdr = hdulist[0].header
hdulist.close()
prihdr['rotoff'] = None 

#reads in files
filelist = glob.glob(pathToFiles + '/*.fits')
dataset = MagAO.MagAOData(filelist)

#set iwa
dataset.IWA = iwa

xDim = dataset._input.shape[2]
yDim = dataset._input.shape[1]
owa = min(xDim,yDim)/2
#creates cube to eventually hold average SNR data
snrCube = np.zeros((len(klmodes),int((subsections_stop-subsections_start)/subsections_inc+1), int((annuli_stop-annuli_start)/annuli_inc+1),int((movement_stop-movement_start)/movement_inc+1)))


#loop over annuli, movement, and subsection parameters

#keeps track of number of annuli values that have been tested, used for indexing
acount = 0

for a in range(annuli_start, annuli_stop+1, annuli_inc):
    
    dr = float(owa-iwa)/a
    all_bounds = [dr*rad+iwa for rad in range(a+1)]
    numAnn = a
    
    if(singleAnn):
        lowBound = max([b for b in all_bounds if (min(ra)>b)])
        upBound = min([b for b in all_bounds if (max(ra)<b)])
        all_bounds = [b for b in all_bounds if (b>=lowBound and b<=upBound)]
        numAnn = int((upBound-lowBound)/dr)
        dataset.IWA = lowBound
        dataset.OWA = upBound
        

        print("planet at: " + str(ra))
        print("lower bound: " + str(lowBound))
        print("upper bound: " + str(upBound))
        print("number of annuli: " + str(numAnn))
        
    #check to see if any planets fall very close to a zone boundary    
    for pl in ra:
        for b in all_bounds:
            if (b <= pl+FWHM/2 and b >= pl-FWHM/2):
                #dont run
    
    
    if( len([b for b in all bounds if ((b <= pl+FWHM/2 and b >= pl-FWHM/2) for pl in ra)]) == 0):
    
    
    #if ( (min(ra)-FWHM/2 >= lowBound)  and (max(ra)+FWHM/2 <= upBound) ):
    
        #keeps track of number of movement values that have been tested, used for indexing
        mcount = 0

        for m in np.arange(movement_start, movement_stop+1, movement_inc):

            scount = 0

            for s in range(subsections_start, subsections_stop+1, subsections_inc):

                print("Parameters: annuli = %d; movement = %s; subections = %d" %(a, m,s))

                #cube to hold median combinations of klipped images
                cube = np.zeros((len(klmodes),yDim,xDim))
                #creates cube to hold snr maps 
                snrMapCube = np.zeros((len(klmodes),yDim,xDim))


                runKLIP = True

                if (os.path.isfile(str(pathToFiles) + "_klip/med_" + outputFileName + "_a" + str(a) + "m" + str(m) + "s" + str(s) + "iwa" + str(iwa) + '_klmodes-all.fits')):
                    hdulist = fits.open(str(pathToFiles) + "_klip/med_" + outputFileName + "_a" + str(a) + "m" + str(m) + "s" + str(s) + "iwa" + str(iwa) + '_klmodes-all.fits')
                    klmodes2 = hdulist[0].header['klmodes'][1:-1]
                    klmodes2 = list(map(int, klmodes2.split(",")))

                    if (len([k for k in klmodes if not k in klmodes2]) == 0):
                        print("Found KLIP processed images for same parameters saved to disk. Reading in data.")
                        runKLIP = False 
                        for i in range(len(klmodes)):
                            cube[i,:,:] = hdulist[0].data[klmodes2.index(klmodes[i]),:,:]

                

                if (runKLIP):
                    print("Starting KLIP")
                    #run klip for given parameters
                    parallelized.klip_dataset(dataset, outputdir=(pathToFiles + "_klip/"), fileprefix=outputFileName, annuli=numAnn, subsections=s, movement=m, numbasis=klmodes, calibrate_flux=True, mode="ADI") 
                    #flips images
                    output = dataset.output[:,:,:,::-1]

                #keeps track of number of KL mode values that have been tested, used for indexing
                kcount = 0

                #iterates over kl modes
                for k in klmodes:

                    if (runKLIP):
                        #takes median combination of cube made with given number of KL modes
                        isolatedKL = np.nanmedian(output[kcount,:,:,:], axis=0)
                        #adds median image to cube 
                        cube[kcount,:,:] = isolatedKL

                    else:
                        isolatedKL = cube[kcount,:,:]

                    #makes SNR map 
                    snrmap = snr.create_map(isolatedKL,FWHM, smooth = _smooth, planets = mask, saveOutput = False)

                    #adds SNR map to 5d cube 
                    if (saveSNR):
                        snrMapCube[kcount,:,:] = snrmap 


                    planetSNRs = [snr.getPlanet(snrmap, 0, ra[x], pa[x], int(wid[0]/2)+1) for x in range (len(ra))]
                    planetSNR = np.mean(planetSNRs)

                    #add planet snr value to snrCube
                    snrCube[kcount,scount,acount,mcount] = planetSNR
                    kcount+=1

                if(runKLIP):
                    #write median combination cube to disk 
                    print("Writing median image combinations to " + pathToFiles + "_klip/")
                    writeData(snrCube, pre = 'med_')

                if (saveSNR):
                    print("Writing SNR maps to " + pathToFiles + "_klip/")
                    writeData(snrCube, snrmap = True, pre = 'snrmap_')
                print()

                scount+=1
            mcount+=1
    acount+=1


         
print("Writing average SNR values to " + pathToFiles + "_klip/")    
#write snr cube to disk 
writeData(snrCube, allParams = True, snrmap = True, pre = 'paramexplore_')

print()
print("KLIP automation complete")  
            