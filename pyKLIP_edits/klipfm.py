import glob
import numpy as np
import pyklip.instruments.MAGAO as MAGAO
from astropy.io import fits
import pyklip.fmlib.fmpsf as fmpsf
import pyklip.fm as fm


# read in the data into a dataset
filelist = glob.glob("sliced2/*.fits")
dataset = MAGAO.MAGAOData(filelist)
print("Dataset.wvs  = " + str(dataset.wvs))


# read in instrumental PSF
#hdulist = fits.open("InstrumentalPSF.fits")
hdulist = fits.open("IPSF2.fits")
dataset.psf = np.array(hdulist[0].data)
hdulist.close()

# setup FM guesses
numbasis = np.array([1, 5, 10, 20, 50, 100]) # KL basis cutoffs you want to try
guesssep = 20 # estimate of separation in pixels
guesspa = 120 # estimate of position angle, in degrees
guessflux = 1e-2 # estimated contrast
dn_per_contrast = 1 # DN/contrast ratio. For GPI, this is dataset.dn_per_contrast
guessspec = np.array([1]) # should be 1-D array with number of elements = np.size(np.unique(dataset.wvs))

# initialize the FM Planet PSF class
fm_class = fmpsf.FMPlanetPSF(dataset.input.shape, numbasis, guesssep, guesspa, guessflux, dataset.psf,
                             np.unique(dataset.wvs), dn_per_contrast, star_spt='A6', spectrallib=[guessspec])


# PSF subtraction parameters
# You should change these to be suited to your data!
outputdir = "." # where to write the output files
prefix = "klip-fm_test" # fileprefix for the output files
annulus_bounds = [[guesssep-15, guesssep+15]] # one annulus centered on the planet
subsections = 1 # we are not breaking up the annulus
padding = 0 # we are not padding our zones
movement = 0 # we are using an conservative exclusion criteria of 4 pixels

# run KLIP-FM
fm.klip_dataset(dataset, fm_class, outputdir=outputdir, fileprefix=prefix, numbasis=numbasis,
                annuli=annulus_bounds, subsections=subsections, padding=padding, movement=movement)






print("DONE RUNNING KLIPFM.PY")




# Mysterious "fm_class.save_fmout" is where it's failing. I have no idea where this file/library/package is #


# Display an image of the FMPSF from fm_class after it's made








#On evernote, write up a note about SNR and about pyklip
#Describe MAGAO.py and .ini, what's in both of them, how they're different from gpi etc. (1 wavelength at a time isntead of 39), what changes ive made to the pyklip code, 
