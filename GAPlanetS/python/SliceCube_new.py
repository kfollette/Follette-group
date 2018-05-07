#Take in data cube
#Slice it in to X individual files, print corresponding rotoff parameters to headers from separate file

from astropy.io import fits
import os
import numpy as np

def slice(filename, output='sliced', rotoff='rotoff_nocosmics.fits'):
    #os.chdir(str(filepath))
    hdulist = fits.open(str(filename))
    cube = hdulist[0].data
    header = hdulist[0].header
    wv = header['WLENGTH']
    hdulist.close()
    dim = cube.shape[1]

    hdulist = fits.open(str(rotoff))
    rotoffs = hdulist[0].data
    hdulist.close()

    if not os.path.exists(str(output)):
        os.makedirs(str(output))

    current = os.getcwd
    os.chdir(str(output))
    if len(cube) != len(rotoffs):
        print(len(cube),len(rotoffs))
        print("the specified rotoff cube is not the same length as the z dimension of the image cube")

    else:
        for z in range(len(rotoffs)):
            newFITS = np.zeros((dim, dim))
            for y in range(dim):
                for x in range(dim):
                    newFITS[y][x] = cube[z][y][x]
            hdu = fits.PrimaryHDU(newFITS)
            hdulist = fits.HDUList([hdu])
            prihdr = hdulist[0].header
            prihdr.set('rotoff', str(rotoffs[z]))
            prihdr.set('WLENGTH', wv)
            hdulist.writeto("sliced_"+str(z+1)+".fits", clobber=True)
            print("Sliced image " + str(z+1))

    print("Done slicing")
    os.chdir(current)


