"""
redo of python code using image_registration.
will eventually be added to to do the circular symmetry thing as well.

Registers images

Last modified: 7/1/17
"""
from astropy.io import fits
import numpy as np
import image_registration as ir

def center_image(file, ref):
    global original 
    #can add in special filename stuff like in the original later; keeping things simple to test
    hdulist = fits.open(file)
    original = hdulist[0].data
    original = original[0:300]
    header_data = fits.getheader(file,0)
    hdulist.close()
    print("read in cube")
    
    global reg_cube
    reg_cube = register(original, ref)
    print("registered cube")
    center = find_max(reg_cube[ref])
    print(center)
    clip_cube = clip(reg_cube, center[0], center[1])
    print("writing to fits file")
    hdu = fits.PrimaryHDU(clip_cube, header_data)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("shifted_cube.fits", overwrite=True)
    #med = np.nanmedian(original, axis=0)
    #print("constructed median image")
    
    #center = sym_center(med)

def register(cube, ref):
    size = cube.shape[0]
    for z in range(size):
        shift = ir.chi2_shift(cube[ref],cube[z])
        shiftX, shiftY = shift[0], shift[1]
        
        #print(shift)
        
        shifted_image = ir.fft_tools.shift2d(cube[z],-shiftX,-shiftY)
        cube[z] = shifted_image
        if(z%100 == 0):
            print("shift number " + str(z))
    return cube

def clip(cube, centerX, centerY):
    clipped_cube = np.zeros((cube.shape[0],451,451))
    for z in range(cube.shape[0]):
        for y in range(-225,226):
            for x in range(-225,226):
                clipped_cube[z][y+225][x+225] = cube[z][centerY+y][centerX+x]
    return clipped_cube

def find_max(image):
    m = 0
    centerx = 0
    centery = 0
    for y in range(int(image.shape[0]/2-100),int(image.shape[0]/2+100)):
        for x in range(int(image.shape[1]/2-100),int(image.shape[1]/2+100)):
            if(image[y][x] > m):
                m = image[y][x]
                centerx = x
                centery = y
    return centerx, centery
            
"""    
def sym_center(image):
"""  

###TESTING###
center_image("HD142527_8Apr14_Line_flat_preproc.fits",2)