#Elijah Spiro
#Version 1.0 - 9/17/16

import glob
import pyklip.instruments.GPI as GPI
import pyklip.instruments.MAGAO as MAGAO
import pyklip.parallelized as parallelized
import numpy as np
import sys
import pyklip.klip as klip
from astropy.io import fits

pathToFiles = sys.argv[1]
annuli2 = int(sys.argv[2])
IWA2 = int(sys.argv[3])
movement2 = float(sys.argv[4])
outputName = sys.argv[5] 
klmodes = sys.argv[6].split(',')
klmodes = list(map(int, klmodes))
subsections2 = int(sys.argv[7])

print(pathToFiles + "/*.fits")

filelist = glob.glob(pathToFiles + "/*.fits")
dataset = MAGAO.MAGAOData(filelist)

print("KLMODES ARE " + str(klmodes))

parallelized.klip_dataset(dataset, outputdir="", fileprefix=outputName, annuli=annuli2, subsections=subsections2, movement=movement2, numbasis=klmodes, calibrate_flux=True, mode="ADI")

avgframe = np.nanmean(dataset.output[1], axis=(0,1))
print("Shape of avgframe is " + str(avgframe.shape))
calib_frame = dataset.calibrate_output(avgframe)

print("Completed klipping. Rotating images")
hdulist = fits.open(outputName+"-KLmodes-all.fits")
cube = hdulist[1].data
hdulist.close()
cube = cube[:,:,::-1]
hdu = fits.PrimaryHDU(cube)
hduList = fits.HDUList([hdu])
prihdr = hduList[0].header
prihdr['Annuli'] = annuli2
prihdr['IWA'] = IWA2
prihdr['Movement'] = movement2
prihdr['KLModes'] = str(klmodes)
prihdr['Subsects'] = subsections2
hduList.writeto("_"+outputName+"-KLmodes-all.fits")

print("Complete")
