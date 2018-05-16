## Make post-KLIP SDI images for parameter exploration
## Written by Kate Follette 5/14/18

import numpy as np
from astropy.io import fits

annuli = [1,2]
movement = [1.,2.]
#annuli = list(np.arange(0,36))
#movement = list(np.arange(0,23))
subsections = [1]
iwa = [0]
scale=[1]
klmodes = [1,2,5,10,20,50,100]
filename = 'LkCa15_18Nov16'
filepath='/Users/kfollette/Shared_Dropbox/Follette-Lab/GAPlanetS_Data/GAPlanetS_2017Reprocessing/LkCa15/18Nov16/merged/fwhm38cut/' 
#filename='LkCa15_16Nov14_1e1fake_nobpl'
#filepath='/Users/kfollette/Shared_Dropbox/Follette-Lab/GAPlanetS_Data/GAPlanetS_2017Reprocessing/LkCa15/16Nov14/wfe171pt9/'

for a in np.arange(len(annuli)):
    for m in np.arange(len(movement)):
        for s in np.arange(len(subsections)):
            print(a,m,s)
            line = fits.getdata(str(filepath)+'sliced_Line_KLIP/' +'med_'  + filename +'_Line'+ '_a' + str(annuli[a]) + "m" + str("{:.1f}".format(movement[m])) + "s" + str(subsections[s]) + "iwa" + str(iwa[0])  +  '_KLmodes-all.fits', clobber=True)
            hdr = fits.getheader(str(filepath)+'sliced_Line_KLIP/' +'med_'  + filename +'_Line'+ '_a' + str(annuli[a]) + "m" + str("{:.1f}".format(movement[m])) + "s" + str(subsections[s]) + "iwa" + str(iwa[0])  +  '_KLmodes-all.fits', clobber=True)
            cont=fits.getdata(str(filepath)+'sliced_Cont_KLIP/' +'med_'  + filename +'_Cont'+ "_a" + str(annuli[a]) + "m" + str("{:.1f}".format(movement[m])) + "s" + str(subsections[s]) + "iwa" + str(iwa[0])  +  '_KLmodes-all.fits', clobber=True)
            SDI = line - scale*cont
    #        hdu = fits.PrimaryHDU(SDI)
     #       hdulist=fits.HDUList([hdu])
            fits.writeto(str(filepath)+'sliced_SDI_KLIP/med_SDI'+"_a" + str(annuli[a]) + "m" + str("{:.1f}".format(movement[m])) + "s" + str(subsections[s]) + "iwa" + str(iwa[0]) +'_KLmodes-all.fits', SDI, hdr, clobber=True) 

