#Take in data cube
#Slice it in to X individual files, print corresponding rotoff parameters to headers from separate file

from astropy.io import fits
import os
import numpy as np

print("Enter absolute directory filepath")
dir = input()
print("Enter name for output directory ")
outdir = input()
print("Enter filename to be sliced ")
fname = input()
print("Enter name of rotoff cube ")
rotname = input()

os.chdir(str(dir))
hdulist = fits.open(str(fname))
cube = hdulist[0].data
hdulist.close()
dim = cube.shape[1]

hdulist = fits.open(str(rotname))
rotoffs = hdulist[0].data
hdulist.close()

if not os.path.exists(str(outdir)):
    os.makedirs(str(outdir))

os.chdir(str(outdir))
if len(cube) != len(rotoffs):
    print("the specified rotoff cube is not the same length as the z dimension of the image cube")

else:
    for z in range(len(rotoffs)):
        newFITS = np.zeros(dim, dim)
        for y in range(dim):
            for x in range(dim):
                newFITS[y][x] = cube[z][y][x]
        hdu = fits.PrimaryHDU(newFITS)
        hdulist = fits.HDUList([hdu])
        prihdr = hdulist[0].header
        prihdr.set('rotoff', str(rotoffs[z]))
        hdulist.writeto("sliced_"+str(z+1)+".fits", clobber=True)
        print("Sliced image " + str(z+1))

print("Done slicing")
