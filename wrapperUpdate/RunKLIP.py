##################################################################
#############                                        #############
#############                 IMPORTS                #############
#############                                        #############
##################################################################

import glob
import inspect
import os
#import pyklip.instruments.MAGAO as MAGAO
import instruments.MagAO as MagAO
#import pyklip.parallelized as parallelized
import parallelized as parallelized
import numpy as np
import sys
#import pyklip.klip as klip
import klip as klip
from astropy.io import fits
import warnings
from astropy.utils.exceptions import AstropyWarning
import SNRMap as snr

warnings.filterwarnings('ignore', category=AstropyWarning, append=True)


##################################################################
#############                                        #############
#############               SAVE FILES               #############
#############                                        #############
################################################################## 
 
def writeData(indiv, filepath, filename, annuli, movement, subsections, iwa, klmodes, mask = None, pre = '', suff = ""):    
    #function writes fits files and adds important info to headers
    hdu = fits.PrimaryHDU(indiv)
    hdulist = fits.HDUList([hdu])
    prihdr = hdulist[0].header
    
    #shorten file path to bottom 4 directories so it will fit in header
    pathToFiles_short  = ''
    numdir = 0
    for n in range (len(filepath)):
        pathToFiles_short = filepath[-1-n] + pathToFiles_short
        if (filepath[-1-n] == '/'):
            numdir += 1
        if (numdir >=4 ):
            break
    
    #add info to fits header
    prihdr.set('annuli', str(annuli))
    prihdr.set('movement', str(movement))
    prihdr.set('subsctns', str(subsections))
    prihdr.set('klmodes', str(klmodes))
    prihdr.set('filepath', str(pathToFiles_short))
    if (not mask == None):
        rad, pa, wid = mask 
        prihdr.set('mask_rad', str(rad))
        prihdr.set('mask_pa', str(pa))
        prihdr.set('mask_wid', str(wid))
    
    #write fits files
    hdulist.writeto(str(filepath) + "_KLIP/" + str(pre) + filename + "_a" + str(annuli) + "m" + str(movement) + "s" + str(subsections) + "iwa" + str(iwa) + str(suff) + '_KLmodes-all.fits' , clobber=True)

    

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

print("KLIP Parameters:")
try:
    while (not pathToFiles[-7:] == '%finish'):
        argnum += 1
        pathToFiles = pathToFiles + " " + sys.argv[1+argnum]
    pathToFiles = pathToFiles[:-7]
        
except:
    pathToFiles = sys.argv[1]
    argnum = 0

print("  File Path = " + pathToFiles) 

if not os.path.exists(pathToFiles + "_KLIP"):
    os.makedirs(pathToFiles + "_KLIP")
    os.chmod(pathToFiles + "_KLIP", 0o777)

iwa = int(sys.argv[2+argnum])
print("  IWA = " + str(iwa))

klmodes = list(map(int, sys.argv[3+argnum].split(",")))
print("  KL Modes = " + str(list(map(int, sys.argv[3+argnum].split(",")))))

annuli2 = int(sys.argv[4+argnum])
print("  Annuli = " + str(annuli2))

movement2 = float(sys.argv[5+argnum])
print("  Movement = " + str(movement2))

subsections2 = int(sys.argv[6+argnum])
print("  Subsections = " + str(subsections2))

outputFileName = sys.argv[7+argnum]     

SNR = False
if (sys.argv[8+argnum] == 'true' or sys.argv[8+argnum] == 'True'):
    SNR = True
    
saveData = False
if (sys.argv[9+argnum] == 'true' or sys.argv[9+argnum] == 'True'):
    saveData = True    
    
maskParams = None

if (SNR):
    try: 
        print()
        print("SNR Map Parameters:")
        FWHM = float(sys.argv[10+argnum])
        print('  Star FWHM = ' + str(FWHM))
        print('  Planet mask parameters:')
        print("    Radius = " + str(list(map(int, sys.argv[11+argnum].split(",")))))
        ra = list(map(int, sys.argv[11+argnum].split(",")))
        print("    Position Angle = " + str(list(map(int, sys.argv[12+argnum].split(",")))))
        pa = list(map(int, sys.argv[12+argnum].split(",")))
        print("    Mask width (radial, angular): = " + str(list(map(int, sys.argv[13+argnum].split(",")))))
        wid = list(map(int, sys.argv[13+argnum].split(",")))
        maskParams = (ra, pa, wid)
        
    except:
        pass
        


    
##################################################################
#############                                        #############
#############                 RUN KLIP               #############
#############                                        #############
##################################################################

print()
print("Reading: " + pathToFiles + "/*.fits")
filelist = glob.glob(pathToFiles + '/*.fits')

dataset = MagAO.MagAOData(filelist)
#set iwa
dataset.IWA = iwa

print()
print("Starting KLIP")


#run klip for given parameters
parallelized.klip_dataset(dataset, outputdir=(pathToFiles + "/.."), fileprefix=outputFileName, annuli=annuli2, subsections=subsections2, movement=movement2, numbasis=klmodes, calibrate_flux=True, mode="ADI")
           
#cube to hold median combinations of klipped images
dim = dataset.output.shape[2]

cube = np.zeros((len(klmodes),dim,dim))

#cube to hold SNR maps
SNRcube = np.zeros((len(klmodes),dim,dim))
            
#flips images
print("Now flipping KLIPed images")
dataset.output = dataset.output[:,:,:,::-1]
      
if (saveData):
    print("Writing KLIPed time series 4D cube to " + pathToFiles + "_KLIP")
    writeData(dataset.output, pathToFiles, outputFileName, annuli2, movement2, subsections2, iwa, klmodes, suff = "_uncombined")
    


#keeps track of number of KL mode values that have been tested, used for indexing
kcount = 0                       
#iterates over kl modes
for k in klmodes:

    #takes median combination of cube made with given number of KL modes
    isolatedKL = np.nanmedian(dataset.output[kcount,:,:,:], axis=0)
    #adds median image to cube 
    cube[kcount,:,:] = isolatedKL
    
    #creates SNR map if designated
    if (SNR):  
        print()
        print("Runing SNRMap on KLIPed data")
        SNRcube[kcount,:,:] = snr.create_map(isolatedKL, FWHM, planets = maskParams, saveOutput = False)
                
    kcount += 1
       
        
#write median combination cube to disk 
print()
print("Writing median KLIPed images to " + pathToFiles + "_KLIP")
writeData(cube, pathToFiles, outputFileName, annuli2, movement2, subsections2, iwa, klmodes, pre = "med_")

  
if (SNR):
    print("Writing SNR maps to " + pathToFiles + "_KLIP")
    writeData(SNRcube, pathToFiles, outputFileName, annuli2, movement2, subsections2, iwa, klmodes, mask = maskParams, pre = "SNRMap_")

        
print("KLIP completed")        

  
     
            





