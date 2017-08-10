#slice up h-a data cube, feed that in and get it to run (incorporate rotoff into corresponding headers)

import glob
import pyklip.instruments.GPI as GPI
import pyklip.instruments.MAGAO as MAGAO
import pyklip.parallelized as parallelized
import numpy as np
import pyklip.klip as klip
from astropy.io import fits

#filelist = glob.glob("20141218_H_Spec/*.fits")
print("Enter the path to the working sliced directory (ending with 'sliced/'):")
dir = input()
print("Enter a name for the ouput file:")
name = input()
print("Enter number of annuli:")
annuli2 = input()
print("Enter movement parameter:")
movement2 = input()
filelist = glob.glob(str(dir) + "*.fits")
#filelist = glob.glob("spiral/sliced/*.fits")
dataset = MAGAO.MAGAOData(filelist)
#dataset = GPI.GPIData(filelist)

outputFileName = str(name)

parallelized.klip_dataset(dataset, outputdir="", fileprefix=outputFileName, annuli=int(annuli2), subsections=1, movement=int(movement2), numbasis=[1,5,10,50], calibrate_flux=True, mode="ADI")

print("Shape of dataset.output is " + str(dataset.output.shape))
print("Shape of dataset.output[1] is " + str(dataset.output[1].shape))
avgframe = np.nanmean(dataset.output[1], axis=(0,1))
print("Shape of avgframe is " + str(avgframe.shape))
calib_frame = dataset.calibrate_output(avgframe)

print("Shape of calib_frame: " + str(calib_frame.shape))
#seps, contrast = klip.meas_contrast(calib_frame, dataset.IWA, 1.1/GPI.GPIData.lenslet_scale, 3.5)

print("Completed klipping. Rotating images")
hdulist = fits.open(outputFileName+"-KLmodes-all.fits")
cube = hdulist[1].data
hdulist.close()
cube = cube[:,:,::-1]
"""
newCube = []
for i in range(len(cube)):
    newCube.append(cube[len(cube)-i-1])
"""
hdulist = fits.PrimaryHDU(cube)
hdulist.writeto(outputFileName+"-KLmodes-all.fits", clobber=True)


cubetwo = dataset.output[0,:,:,:]
cubetwo = cubetwo[:,:,::-1]
med=np.nanmedian(cubetwo, axis=0)
fits.writeto('med_'+ outputFileName + '-1KLmodes.fits', med, clobber=True)

#cubethree = dataset.output[1,:,:,:]
#cubethree = cubethree[:,:,::-1]
#med=np.nanmedian(cubethree, axis=0)
#fits.writeto('med_'+ outputFileName + '-2KLmodes.fits', med, clobber=True)

#cubefour = dataset.output[2,:,:,:]
#cubefour = cubefour[:,:,::-1]
#med=np.nanmedian(cubefour, axis=0)
#fits.writeto('med_'+ outputFileName + '-3KLmodes.fits', med, clobber=True)

#cubefive = dataset.output[3,:,:,:]
#cubefive = cubefive[:,:,::-1]
#med=np.nanmedian(cubefive, axis=0)
#fits.writeto('med_'+ outputFileName + '-4KLmodes.fits', med, clobber=True)

cubesix = dataset.output[1,:,:,:]
cubesix = cubesix[:,:,::-1]
med=np.nanmedian(cubesix, axis=0)
fits.writeto('med_'+ outputFileName + '-5KLmodes.fits', med, clobber=True)

cubeseven = dataset.output[2,:,:,:]
cubeseven = cubeseven[:,:,::-1]
med=np.nanmedian(cubeseven, axis=0)
fits.writeto('med_'+ outputFileName + '-10KLmodes.fits', med, clobber=True)

#cubeeight = dataset.output[6,:,:,:]
#cubeeight = cubeeight[:,:,::-1]
#med=np.nanmedian(cubeeight, axis=0)
#fits.writeto('med_'+ outputFileName + '-20KLmodes.fits', med, clobber=True)

cubenine = dataset.output[3,:,:,:]
cubenine = cubenine[:,:,::-1]
med=np.nanmedian(cubenine, axis=0)
fits.writeto('med_'+ outputFileName + '-50KLmodes.fits', med, clobber=True)

#cubeten = dataset.output[3,:,:,:]
#cubeten = cubeten[:,:,::-1]
#med=np.nanmedian(cubeten, axis=0)
#fits.writeto('med_'+ outputFileName + '-100KLmodes.fits', med, clobber=True)


print("Complete")
