##################################################################
#############                                        #############
#############                 IMPORTS                #############
#############                                        #############
##################################################################

import glob
import inspect
import os
#import pyklip.instruments.MAGAO as MAGAO
import MagAO as MagAO
#import pyklip.parallelized as parallelized
import parallelized as parallelized
import numpy as np
import sys
#import pyklip.klip as klip
import klip as klip
from astropy.io import fits


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
while (not pathToFiles[-1] == 'd' and not pathToFiles[-1] == '"'):
    argnum += 1
    pathToFiles = pathToFiles + " " + sys.argv[1+argnum]

klmodes = list(map(int, sys.argv[6+argnum].split(",")))
print("KL Modes = " + str(list(map(int, sys.argv[6+argnum].split(",")))))

iwa = int(sys.argv[3+argnum])
print("IWA = " + str(iwa))

annuli2 = int(sys.argv[2+argnum])
print("Annuli = " + str(annuli2))

movement2 = float(sys.argv[4+argnum])
print("Movement = " + str(movement2))

subsections2 = int(sys.argv[7+argnum])
print("Subsections = " + str(subsections2))

outputFileName = sys.argv[5+argnum]                                                                                                                                                                                                                                                                                                                                                                     
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

#run klip for given parameters
parallelized.klip_dataset(dataset, outputdir=(pathToFiles + "/.."), fileprefix=outputFileName, annuli=annuli2, subsections=subsections, movement=movement2, numbasis=klmodes, calibrate_flux=True, mode="ADI")
           
#cube to hold median combinations of klipped images
cube = np.zeros((5,450,450))
            
#flips images
print("Now flipping KLIPed images"
dataset.output = dataset.output[:,:,:,::-1]

#keeps track of number of KL mode values that have been tested, used for indexing
kcount = 0                       
#iterates over kl modes
for k in klmodes:

    #takes median combination of cube made with given number of KL modes
    isolatedKL = np.nanmedian(dataset.output[kcount,:,:,:], axis=0)

    #adds median image to cube 
    cube[kcount,:,:] = isolatedKL
                
    kcount+=1
       
        
#write median combination cube to disk 
fits.writeto('med_'+ outputFileName + "_a" + str(a) + "m" + str(m) + "s" str(subsections2) + iwa + str(iwa) + '_KLmodes-all.fits', cube, clobber=True)
  

     
            





