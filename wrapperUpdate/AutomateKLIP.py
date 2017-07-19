
# coding: utf-8

# In[ ]:

#Clare Leonard                                                                
#Version 1.0 - 6/21/17                                                          

import glob
import inspect
import os                                      
import MAGAO as MAGAO                                   
import parallelized as parallelized
import numpy as np
import sys                                                   
import klip as klip
from astropy.io import fits
import SNRMap as snr   
import time


#get inputs 

pathToFiles = sys.argv[1]
print("File Path = " + pathToFiles)

print()

print("Parameters to explore:")

annuli2_start = int(sys.argv[5])
annuli2_stop = int(sys.argv[6])
annuli2_inc = int(sys.argv[7])

if(annuli2_start == annuli2_stop):
    annuli2_inc = 1;
    print("Annuli = " +str(annuli2_start))
    
else:
    print("Annuli: start = %s; end = %s; increment = %s " %(str(annuli2_start), str(annuli2_stop), str(annuli2_inc)))



movement2_start = int(sys.argv[8])
movement2_stop = int(sys.argv[9])
movement2_inc = int(sys.argv[10])


if(movement2_start == movement2_stop):
    movement2_inc = 1;
    print("Movement = " +str(movement2_start))
    
else:
    print("Movement: start = %s; end = %s; increment = %s " %(str(movement2_start), str(movement2_stop), str(movement2_inc)))


subsections2_start = int(sys.argv[11])
subsections2_stop = int(sys.argv[12])
subsections2_inc = int(sys.argv[13])


if(subsections2_start == subsections2_stop):
    subsections2_inc = 1;
    print("Subsections = " +str(subsections2_start))
    
else:
    print("Subsections: start = %s; end = %s; increment = %s " %(str(subsections2_start), str(subsections2_stop), str(subsections2_inc)))

print()
    
iwa = int(sys.argv[2])
print("IWA = " + str(iwa))

print()

print("KL Modes = " + str(list(map(int, sys.argv[3].split(",")))))
klmodes = list(map(int, sys.argv[3].split(",")))

outputFileName = sys.argv[4]
#print("Output FileName = " + outputFileName)


print()

print("reading: " + pathToFiles + "/*.fits")

print()

filelist = glob.glob(pathToFiles + "/*.fits")
dataset = MAGAO.MAGAOData(filelist)

print()

snrCube = np.zeros((len(klmodes),int((annuli2_stop-annuli2_start)/annuli2_inc+1),int((movement2_stop-movement2_start)/movement2_inc+1)))











#loop over annuli, movement, and subsection parameters

print("running klip for parameters:")


#keeps track of number of annuli values that have been tested, used for indexing
acount = 0

for a in range(annuli2_start, annuli2_stop+1, annuli2_inc):
    
    #keeps track of number of movement values that have been tested, used for indexing
    mcount = 0
    
    for m in range(movement2_start, movement2_stop+1, movement2_inc):
        
        for s in range(subsections2_start, subsections2_stop+1, subsections2_inc):
            print("annuli = %d; movement = %d; subections = %d" %(a, m,s))
            #sys.stdout.flush()

            #run klip for given parameters
            parallelized.klip_dataset(dataset, outputdir="", fileprefix=outputFileName, annuli=a, subsections=s, movement=m, numbasis=klmodes, calibrate_flux=True, mode="ADI")

            
            #avgframe = np.nanmedian(dataset.output[1], axis=(0,1))
            #print("Shape of avgframe is " + str(avgframe.shape))
            #calib_frame = dataset.calibrate_output(avgframe)

            
            #cube to hold median combinations of klipped images
            cube = np.zeros((5,450,450))
            
            #keeps track of number of KL mode values that have been tested, used for indexing
            kcount = 0
            
            #flips images
            dataset.output = dataset.output[:,:,:,::-1]
            
            #iterates over kl modes
            for k in klmodes:
                
                #takes median combination of cube made with given number of KL modes
                isolatedKL = dataset.output[kcount,:,:,:]
                isolatedKL = np.nanmedian(isolatedKL, axis=0)
                
                #put together output name
                outputNameSNR = outputFileName + "_a" + str(a) + "m" + str(m) + "KL" + str(k) + "_SNRMap.fits"
                
                #object to hold mask parameters for snr map 
                mask = ([13,], [120,], [10, 15])
                
                #gets highest pixel value in snr map in the location of the planet 
                planetSNR = snr.getPlanet(snr.create_map(isolatedKL, planets = mask, outputName = outputNameSNR), 220, 215, 10)
                
                #adds median image to cube 
                cube[kcount,:,:] = isolatedKL
                
                #add planet snr value to snrCube
                snrCube[kcount,acount,mcount] = planetSNR
                kcount+=1
                                          
            #write median combination cube to disk 
            fits.writeto('med_'+ outputFileName + "_a" + str(a) + "m" + str(m) + '-KLmodes-all.fits', cube, clobber=True)
        mcount+=1
    acount+=1

#write snr cube to disk 
fits.writeto('med_'+ outputFileName + "_a" + str(annuli2_start) + "-" + str(annuli2_stop) + "x" + str(annuli2_inc) + "m" + str(movement2_start) + "-" + str(movement2_stop) + "x" + str(movement2_inc) + '-KLmodes-all_snrCube.fits', snrCube)     
            

