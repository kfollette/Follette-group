
import glob
import numpy as np
import pyklip.instruments.MagAO as MagAO
from astropy.io import fits
import sys

line = sys.argv[1]
ghost = sys.argv[2]
output = sys.argv[3]
pre = sys.argv[4]
KL = int(sys.argv[5])
sep = float(sys.argv[6])
pa = float(sys.argv[7])
flux = float(sys.argv[8])
annulus1 = sys.argv[9]
annulus2 = sys.argv[10]
move = float(sys.argv[11])
scale = float(sys.argv[12])

annulus = [int(annulus1),int(annulus2)]

filelist = glob.glob(line)
dataset = MagAO.MagAOData(filelist)
psf = fits.getdata(ghost)

psf2 = np.zeros((1,psf.shape[0],psf.shape[1]))
psf2[0] = psf

# setup FM guesses
# You should change these to be suited to your data!
numbasis = np.array([KL]) # KL basis cutoffs you want to try
guesssep = sep # estimate of separation in pixels
guesspa = pa # estimate of position angle, in degrees
guessflux = flux # estimated contrast
dn_per_contrast = np.zeros((dataset.input.shape[0]))
for i in range (dn_per_contrast.shape[0]):
    dn_per_contrast[i] = scale # factor to scale PSF to star PSF. For GPI, this is dataset.dn_per_contrast
guessspec = np.array([1]) # should be 1-D array with number of elements = np.size(np.unique(dataset.wvs))

# initialize the FM Planet PSF class
import pyklip.fmlib.fmpsf as fmpsf
fm_class = fmpsf.FMPlanetPSF(dataset.input.shape, numbasis, guesssep, guesspa, guessflux, psf2, np.unique(dataset.wvs), dn_per_contrast, star_spt='A6', spectrallib=[guessspec])

# PSF subtraction parameters
# You should change these to be suited to your data!
outputdir = output # where to write the output files
prefix = pre # fileprefix for the output files
annulus_bounds = [annulus] # one annulus centered on the planet
subsections = 1 # we are not breaking up the annulus
padding = 0 # we are not padding our zones
movement = move 

# run KLIP-FM
import pyklip.fm as fm
fm.klip_dataset(dataset, fm_class, mode="ADI", outputdir=outputdir, fileprefix=prefix, numbasis=numbasis, annuli=annulus_bounds, subsections=subsections, padding=padding, movement=movement)
