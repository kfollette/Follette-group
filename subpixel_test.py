"""
subpixel shifting
"""

from astropy.io import fits
import sys
import image_registration as ir

###INPUTS###
#will not be command line arguments in the future, this is temporary
#specific frames selected just for test use
hdulist = fits.open(sys.argv[1])
original = hdulist[0].data[500]
later_image = hdulist[0].data[600]
header_data = fits.getheader(sys.argv[1],0)
hdulist.close()

###TEST WITH FAKE OFFSET###
'''
offset = ir.tests.make_offset_extended(original,3.6,3.6,noise=10)
hdu = fits.PrimaryHDU(offset, header_data)
hdulist = fits.HDUList([hdu])
hdulist.writeto("test_offset.fits", overwrite=True)
hdulist.close()


shift = ir.chi2_shift(original,offset)
shiftX, shiftY = shift[0], shift[1]

print(shift)

shifted_image = ir.fft_tools.shift2d(offset,-shiftX,-shiftY)
hdu = fits.PrimaryHDU(shifted_image, header_data)
hdulist = fits.HDUList([hdu])
hdulist.writeto("reshifted.fits", overwrite=True)
'''

###TESTING WITH ACTUAL CUBE IMAGES###

shift = ir.chi2_shift(original,later_image)
shiftX, shiftY = shift[0], shift[1]

print(shift)

shifted_image = ir.fft_tools.shift2d(later_image,-shiftX,-shiftY)
hdu = fits.PrimaryHDU(shifted_image, header_data)
hdulist = fits.HDUList([hdu])
hdulist.writeto("reshifted2.fits", overwrite=True)