
import glob
import numpy as np
import pyklip.instruments.MagAO as MagAO
from astropy.io import fits

filelist = glob.glob("/Users/awatson18/Desktop/Line_sliced/*.fits")
dataset = MagAO.MagAOData(filelist)
psf = fits.getdata("/Users/awatson18/Desktop/test_square2.fits")

psf2 = np.zeros((1,psf.shape[0],psf.shape[1]))
psf2[0] = psf

# setup FM guesses
# You should change these to be suited to your data!
numbasis = np.array([100]) # KL basis cutoffs you want to try
guesssep = 13 # estimate of separation in pixels
guesspa = 120 # estimate of position angle, in degrees
guessflux = 1e-2 # estimated contrast
dn_per_contrast = 2.38 # factor to scale PSF to star PSF. For GPI, this is dataset.dn_per_contrast
#guessspec # should be 1-D array with number of elements = np.size(np.unique(dataset.wvs))

# initialize the FM Planet PSF class
import pyklip.fmlib.fmpsf as fmpsf
fm_class = fmpsf.FMPlanetPSF(dataset.input.shape, numbasis, guesssep, guesspa, guessflux, psf2, np.unique(dataset.wvs), dn_per_contrast, star_spt='A6')
#spectrallib=[guessspec])

# PSF subtraction parameters
# You should change these to be suited to your data!
outputdir = "/Users/awatson18/Desktop/BKA_output" # where to write the output files
prefix = "HD142527_8Apr14_Line_test1" # fileprefix for the output files
annulus_bounds = [[guesssep-15, guesssep+15]] # one annulus centered on the planet
subsections = 1 # we are not breaking up the annulus
padding = 0 # we are not padding our zones
movement = 9 # we are using an conservative exclusion criteria of 4 pixels

# run KLIP-FM
import pyklip.fm as fm
fm.klip_dataset(dataset, fm_class, mode="ADI", outputdir=outputdir, fileprefix=prefix, numbasis=numbasis, annuli=annulus_bounds, subsections=subsections, padding=padding, movement=movement)
