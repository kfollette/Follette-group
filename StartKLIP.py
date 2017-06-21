#Elijah Spiro
#Version 2.1 - 3/26/17

import glob
import inspect
import os
#import pyklip.instruments.MAGAO as MAGAO
import MAGAO as MAGAO
#import pyklip.parallelized as parallelized
import parallelized as parallelized
import numpy as np
import sys
#import pyklip.klip as klip
import klip as klip
from astropy.io import fits
#import SNRMap as SNR


pathToFiles = sys.argv[1]
print("File Path = " + pathToFiles)
annuli2 = int(sys.argv[2])
print("Annuli = " + str(annuli2))
iwa = int(sys.argv[3])
print("IWA = " + str(iwa))
movement2 = float(sys.argv[4])
print("Movement = " + str(movement2))
outputFileName = sys.argv[5] 
print("Output FileName = " + outputFileName)
print("KL Modes = " + str(list(map(int, sys.argv[6].split(",")))))
klmodes = list(map(int, sys.argv[6].split(",")))
subsections2 = int(sys.argv[7])
print("Subsections = " + str(subsections2))

print(pathToFiles + "/*.fits")

filelist = glob.glob(pathToFiles + "/*.fits")
dataset = MAGAO.MAGAOData(filelist)

parallelized.klip_dataset(dataset, outputdir="", fileprefix=outputFileName, annuli=annuli2, subsections=subsections2, movement=movement2, numbasis=klmodes, calibrate_flux=True, mode="ADI")

avgframe = np.nanmedian(dataset.output[1], axis=(0,1))
print("Shape of avgframe is " + str(avgframe.shape))
calib_frame = dataset.calibrate_output(avgframe)

print("Completed klipping. Rotating images")
hdulist = fits.open(outputFileName+"-KLmodes-all.fits")
cube = hdulist[1].data
hdulist.close()
cube = cube[:,:,::-1]
#med=np.nanmedian(cube, axis=0)
hdulist = fits.PrimaryHDU(cube)
hdulist.writeto("_"+outputFileName+"-KLmodes-all.fits")


print("Complete")
