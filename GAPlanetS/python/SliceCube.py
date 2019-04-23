#Take in data cube
#Slice it in to X individual files, print corresponding rotoff parameters to headers from separate file

from astropy.io import fits
import os
import numpy as np
import sys 

def slice(dirPath, filename, output='sliced', rotoff='rotoff_nocosmics.fits'):
    filename = dirPath + '/' + filename
    output = dirPath + '/' + output
    rotoff = dirPath + '/' + rotoff

    hdulist = fits.open(str(filename))
    cube = hdulist[0].data
    header = hdulist[0].header
    hdulist.close()
    dim = cube.shape[1]

    hdulist = fits.open(str(rotoff))
    rotoffs = hdulist[0].data
    hdulist.close()

    if not os.path.exists(str(output)):
        os.makedirs(str(output))

    current = str(os.getcwd())
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
            header.set('rotoff', str(rotoffs[z]))
            hdulist[0].header = header
            hdulist.writeto("sliced_"+str(z+1)+".fits", overwrite=True)
            #print("Sliced image " + str(z+1))

    print("Done slicing")
    os.chdir(current)
    
    
    
def main():
    try:
        dir = sys.argv[1]
        fname = sys.argv[2]
        outdir = sys.argv[3]
        rotname = sys.argv[4]
        
    except:    
        print("Enter absolute directory filepath")
        dir = input()
        print("Enter name for output directory ")
        outdir = input()
        print("Enter filename to be sliced ")
        fname = input()
        print("Enter name of rotoff cube ")
        rotname = input()


    slice(dir, fname, output=outdir, rotoff=rotname)



if __name__ == "__main__": main()
