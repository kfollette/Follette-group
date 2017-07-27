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
#############               GET INPUTS               #############
#############                                        #############
################################################################## 

#value adjusts argument numbering in case of white space in file path 
argnum = 0

pathToFiles = sys.argv[1]
print("File Path = " + pathToFiles)
#if filepath doesnt end in sliced, sontinues to add next arguements, helpful iin case of whitespace in file path
while (not pathToFiles[-6:] == 'sliced'):
    argnum += 1
    pathToFiles = pathToFiles + " " + sys.argv[1+argnum]

klmodes = list(map(int, sys.argv[3+argnum].split(",")))
print("KL Modes = " + str(list(map(int, sys.argv[3+argnum].split(",")))))

iwa = int(sys.argv[2+argnum])
print("IWA = " + str(iwa))

annuli2 = int(sys.argv[4+argnum])
print("Annuli = " + str(annuli2))

movement2 = float(sys.argv[5+argnum])
print("Movement = " + str(movement2))

subsections2 = int(sys.argv[6+argnum])
print("Subsections = " + str(subsections2))

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
        print('Planet mask parameters:')
        print("Radius = " + str(list(map(int, sys.argv[10+argnum].split(",")))))
        ra = list(map(int, sys.argv[10+argnum].split(",")))
        print("Position Angle = " + str(list(map(int, sys.argv[11+argnum].split(",")))))
        pa = list(map(int, sys.argv[11+argnum].split(",")))
        print("Mask width (radial, angular): = " + str(list(map(int, sys.argv[12+argnum].split(",")))))
        wid = list(map(int, sys.argv[12+argnum].split(",")))
        maskParams = (ra, pa, wid)
        
    except:
        pass
        


    
##################################################################
#############                                        #############
#############                 RUN KLIP               #############
#############                                        #############
##################################################################


print("Reading: " + pathToFiles + "/*.fits")
filelist = glob.glob(pathToFiles + '/*.fits')


dataset = MagAO.MagAOData(filelist)

dataset.IWA = iwa

print()
print("Starting KLIP")

#run klip for given parameters
parallelized.klip_dataset(dataset, outputdir=(pathToFiles + "/.."), fileprefix=outputFileName, annuli=annuli2, subsections=subsections2, movement=movement2, numbasis=klmodes, calibrate_flux=True, mode="ADI")
           
#cube to hold median combinations of klipped images
cube = np.zeros((len(klmodes),450,450))

#cube to hold SNR maps
SNRcube = np.zeros((len(klmodes),450,450))
            
#flips images
print("Now flipping KLIPed images")
dataset.output = dataset.output[:,:,:,::-1]
      
if (saveData):
    fits.writeto(pathToFiles + '/../' + outputFileName + "_a" + str(annuli2) + "m" + str(int(movement2)) + "s" + str(subsections2) + "iwa" + str(iwa) +'-KLmodes-all_uncombined.fits', cube, clobber=True)


#keeps track of number of KL mode values that have been tested, used for indexing
kcount = 0                       
#iterates over kl modes
for k in klmodes:

    #takes median combination of cube made with given number of KL modes
    isolatedKL = np.nanmedian(dataset.output[kcount,:,:,:], axis=0)
    
    if (SNR):
        SNRcube[kcount,:,:] = snr.create_map(isolatedKL, planets = maskParams, saveOutput = False)
            
    #adds median image to cube 
    cube[kcount,:,:] = isolatedKL
                
    kcount += 1
       
        
#write median combination cube to disk 
fits.writeto(pathToFiles + '/../med_' + outputFileName + "_a" + str(annuli2) + "m" + str(int(movement2)) + "s" + str(subsections2) + "iwa" + str(iwa) + '_KLmodes-all.fits', cube, clobber=True)
  
if (SNR):
    #write SNR cube to disk 
        fits.writeto(pathToFiles + '/../' + outputFileName + "_a" + str(annuli2) + "m" + str(int(movement2)) + "s" + str(subsections2) + "iwa" + str(iwa) + '_KLmodes-all_SNRMap.fits', SNRcube, clobber=True)
  
     
            





