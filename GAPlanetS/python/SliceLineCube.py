#Take in data cube
#Slice it in to X individual files, print corresponding rotoff parameters to headers from separate file

from astropy.io import fits
import os
import numpy as np

print("Enter directory filepath from " + str(os.getcwd()))
dir = input()
os.chdir(str(dir))
hdulist = fits.open("Line_clip450_flat_reg_circsym_nocosmics.fits")
cube = hdulist[0].data
hdulist.close()

hdulist = fits.open("rotoff_nocosmics.fits")
rotoffs = hdulist[0].data
hdulist.close()

if not os.path.exists("sliced_Line/"):
    os.makedirs("sliced_Line/")

os.chdir("sliced_Line/")

for z in range(len(cube)):
    newFITS = np.zeros((450,450))
    for y in range(450):
        for x in range(450):
            newFITS[y][x] = cube[z][y][x]
    hdu = fits.PrimaryHDU(newFITS)
    hdulist = fits.HDUList([hdu])
    prihdr = hdulist[0].header
    prihdr.set('rotoff', str(rotoffs[z]))
    hdulist.writeto("sliced_"+str(z+1)+".fits")
    print("Sliced line image " + str(z+1))

print("Done slicing")
