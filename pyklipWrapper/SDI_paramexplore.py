## Make post-KLIP SDI images for parameter exploration
## Written by Kate Follette 5/14/18

import numpy as np
from astropy.io import fits

#annuli = [1,2,3,4,5,6]
#movement = [1,2,3,4,5,6,7,8,9,10,11,12]
annuli = list(np.arange(1,31))
movement = list(np.arange(0,23))
subsections = [1]
iwa = [0]
scale=[1.6]
klmodes = [1,2,5,10,20,50,100]
linefilename = 'LkCa15_16Nov14_Line_1e2fake_nobpl'
contfilename = 'LkCa15_16Nov14_Cont'
#filepath='/Users/kfollette/Shared_Dropbox/Follette-Lab/GAPlanetS_Data/GAPlanetS_2017Reprocessing/LkCa15/18Nov16/merged/fwhm38cut/' 
#filename='LkCa15_16Nov14_1e1fake_nobpl'
#filepath='/Users/kfollette/Shared_Dropbox/Follette-Lab/GAPlanetS_Data/GAPlanetS_2017Reprocessing/LkCa15/16Nov14/wfe171pt9/'
filepath='/Users/kfollette/Shared_Dropbox/Follette-Lab/GAPlanetS_Data/GAPlanetS_2017Reprocessing/LkCa15/KLIP_experiments/'
linedir='lkca15_16nov14_sliced_line_1e2fake_pa100_nobpl_KLIP'
contdir='lkca15_16nov14_sliced_cont_KLIP'
outdir='SDI_1e2pl_pa100_KLIP/'

for a in np.arange(len(annuli)):
    for m in np.arange(len(movement)):
        for s in np.arange(len(subsections)):
            print(a,m,s)
            line = fits.getdata(str(filepath)+linedir +'/med_'  + linefilename + '_a' + str(annuli[a]) + "m" + str("{:.1f}".format(movement[m])) + "s" + str(subsections[s]) + "iwa" + str(iwa[0])  +  '_KLmodes-all.fits', clobber=True)
            hdr = fits.getheader(str(filepath)+linedir +'/med_'  + linefilename + '_a' + str(annuli[a]) + "m" + str("{:.1f}".format(movement[m])) + "s" + str(subsections[s]) + "iwa" + str(iwa[0])  +  '_KLmodes-all.fits', clobber=True)
            cont=fits.getdata(str(filepath)+contdir +'/med_'  + contfilename + "_a" + str(annuli[a]) + "m" + str("{:.1f}".format(movement[m])) + "s" + str(subsections[s]) + "iwa" + str(iwa[0])  +  '_KLmodes-all.fits', clobber=True)
            SDI = line - scale*cont
    #        hdu = fits.PrimaryHDU(SDI)
     #       hdulist=fits.HDUList([hdu])
            fits.writeto(str(filepath)+str(outdir)+'med_SDI'+"_a" + str(annuli[a]) + "m" + str("{:.1f}".format(movement[m])) + "s" + str(subsections[s]) + "iwa" + str(iwa[0]) +'_KLmodes-all.fits', SDI, hdr, clobber=True) 

