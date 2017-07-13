
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


# In[ ]:

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

snrCube = np.zeros((5,1,1))


#loop over annuli, movement, and subsection parameters

print("running klip for parameters:")
acount = 0
mcount = 0

for a in range(annuli2_start, annuli2_stop+1, annuli2_inc):
    for m in range(movement2_start, movement2_stop+1, movement2_inc):
        for s in range(subsections2_start, subsections2_stop+1, subsections2_inc):
            sys.stdout.write("annuli = %d; movement = %d; subections = %d                                      \r" %(a, m,s))
            #sys.stdout.flush()

            #run klip
            parallelized.klip_dataset(dataset, outputdir="", fileprefix=outputFileName, annuli=a, subsections=s, movement=m, numbasis=klmodes, calibrate_flux=True, mode="ADI")

            avgframe = np.nanmedian(dataset.output[1], axis=(0,1))
            print("Shape of avgframe is " + str(avgframe.shape))
            calib_frame = dataset.calibrate_output(avgframe)

            cube = np.zeros((5,450,450))
            kcount = 0
            for k in klmodes:
                isolatedKL = dataset.output[0,:,:,:]
                isolatedKL = nanmedian(isolatedKL, axis=0)
                outputNameSNR = outputFileName + "a" + str(a) + "m" + str(m) + "KL" + str(k) + "_SNRMap.fits"
                mask = ([13,], [120,], [10, 15])
                planetSNR = snr.getplanet(snr.create_map(isolatedKL, planets = mask, outputName = outputNameSNR))
                cube[k,:,:] = isolatedKL
                snrCube[kcount,acount,mcount] = planetSNR
                kcount+=1
                                          
            cube = cube[:,:,::-1]
            
            fits.writeto('med_'+ outputFileName + "a" + str(a) + "m" + str(m) + '-KLmodes-all.fits', cube, clobber=True)
        mcount+=1
    acount+=1
            
fits.writeto('med_'+ outputFileName + "a" + str(annuli2_start) + "-" + str(annul12_stop) + "x" + str(annuli2_inc) + "m" + str(annuli2_start) + "-" + str(annul12_stop) + "x" + str(annuli2_inc) + '-KLmodes-all_snrCube.fits', snrCube)     
            

