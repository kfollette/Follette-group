# original implementation written by Alex Watson circa 2017
# expanded by William Balmer
# last edited 9/18/2020
# based on pyklip documentation found here:
# https://pyklip.readthedocs.io/en/latest/bka.html

# import statements
import sys
import glob
import warnings
import numpy as np
from astropy.io import fits

import pyklip.instruments.MagAO as MagAO

import GhostIsolation as Ghost

warnings.filterwarnings('ignore')

if __name__ == '__main__':  # This is a very important precaution for Windows

    if len(sys.argv) > 1:
        # if user wants to run Alex's script from command line they'll have
        # included command line arguments which are accessed by sys.argv

        # Alex's Code

        line = sys.argv[1]  # path to MagAO data
        ghost = sys.argv[2]  # path to instrumental psf
        output = sys.argv[3]  # path to output dir (must exist, doesn't create)
        pre = sys.argv[4]  # prefix for output files
        KL = int(sys.argv[5])  # KL mode
        sep = float(sys.argv[6])  # guess separation
        pa = float(sys.argv[7])  # guess position angle
        flux = float(sys.argv[8])  # guess contrast
        annulus1 = sys.argv[9]  # inner annuli
        annulus2 = sys.argv[10]  # outer annuli
        move = float(sys.argv[11])  # KLIP movement param
        scale = float(sys.argv[12])  # inst psf scale factor

        annulus = [int(annulus1), int(annulus2)]

        filelist = glob.glob(str(line))
        dataset = MagAO.MagAOData(filelist)
        psf = fits.getdata(ghost)

        psf2 = np.zeros((1,psf.shape[0],psf.shape[1]))
        psf2[0] = psf

        # setup FM guesses
        # You should change these to be suited to your data!
        numbasis = np.array([KL])  # KL basis cutoffs you want to try
        guesssep = sep  # estimate of separation in pixels
        guesspa = pa  # estimate of position angle, in degrees
        guessflux = flux  # estimated contrast
        dn_per_contrast = np.zeros((dataset.input.shape[0]))
        for i in range (dn_per_contrast.shape[0]):
            dn_per_contrast[i] = scale  # factor to scale PSF to star PSF.
        guessspec = np.array([1])  # should be 1-D array with number of elements = np.size(np.unique(dataset.wvs))

        # initialize the FM Planet PSF class
        import pyklip.fmlib.fmpsf as fmpsf
        fm_class = fmpsf.FMPlanetPSF(dataset.input.shape, numbasis, guesssep,
                                     guesspa, guessflux, psf2,
                                     np.unique(dataset.wvs), dn_per_contrast,
                                     star_spt='A6', spectrallib=[guessspec])

        # PSF subtraction parameters
        # You should change these to be suited to your data!
        outputdir = output  # where to write the output files
        prefix = pre  # fileprefix for the output files
        annulus_bounds = [annulus]  # one annulus centered on the planet
        subsections = 1  # we are not breaking up the annulus
        padding = 0  # we are not padding our zones
        movement = move

        # run KLIP-FM
        import pyklip.fm as fm
        fm.klip_dataset(dataset, fm_class, mode="ADI", outputdir=outputdir, fileprefix=prefix, numbasis=numbasis, annuli=annulus_bounds, subsections=subsections, padding=padding, movement=movement)

    else:
        # run William's script which asks for inputs

        # check if they have a ghost/psf nearby

        ghostpath = input('Enter the path to your instrumental psf (enter \'none\' to generate one): ')
        if ghostpath == 'none':
            cubepath = input('enter path to your MagAO image cube: ')
            Ghost.ghostIsolation(cubepath, 380, 220, 10, 10, 10)
            ghostpath = 'ghost.fits'
        # ask for path to sliced data
        filepaths = input('Input the filepaths to your sliced data: ')
        filelist = glob.glob(filepaths)
        dataset = MagAO.MagAOData(filelist)
        output = input('Input your output path: ')
        pre = input('Input your output files prefixes: ')
        KL = int(input('Input your desired KL mode: '))
        sep = float(input('Input your separation guess as float: '))
        pa = float(input('Input your position angle guess as float: '))
        flux = float(input('Input your contrast guess as float: '))
        annulus1 = input('Input your inner annulus: ')
        annulus2 = input('Input your outer annulus: ')
        move = float(input('Input KLIP movement parameter: '))
        scale = float(input('Input factor to scale inst. PSF by (usually 1): '))

        annulus = [int(annulus1), int(annulus2)]

        # filelist = glob.glob(str(line))
        # dataset = MagAO.MagAOData(filelist)
        # psf = fits.getdata(ghost)

        # shape instrumental psf how pyklipFM wants
        psf = fits.getdata(ghostpath)
        psf2 = np.zeros((1,psf.shape[0],psf.shape[1]))
        psf2[0] = psf

        # setup FM guesses
        # You should change these to be suited to your data!
        numbasis = np.array([KL])  # KL basis cutoffs you want to try
        guesssep = sep  # estimate of separation in pixels
        guesspa = pa  # estimate of position angle, in degrees
        guessflux = flux  # estimated contrast
        dn_per_contrast = np.zeros((dataset.input.shape[0]))
        for i in range (dn_per_contrast.shape[0]):
            dn_per_contrast[i] = scale  # factor to scale PSF to star PSF.
        guessspec = np.array([1])  # should be 1-D array with number of elements = np.size(np.unique(dataset.wvs))

        # initialize the FM Planet PSF class
        import pyklip.fmlib.fmpsf as fmpsf
        fm_class = fmpsf.FMPlanetPSF(dataset.input.shape, numbasis, guesssep,
                                     guesspa, guessflux, psf2,
                                     np.unique(dataset.wvs), dn_per_contrast,
                                     star_spt='A6', spectrallib=[guessspec])

        # PSF subtraction parameters
        # You should change these to be suited to your data!
        outputdir = output  # where to write the output files
        prefix = pre  # fileprefix for the output files
        annulus_bounds = [annulus]  # one annulus centered on the planet
        subsections = 1  # we are not breaking up the annulus
        padding = 0  # we are not padding our zones
        movement = move

        # run KLIP-FM
        import pyklip.fm as fm
        fm.klip_dataset(dataset, fm_class, mode="ADI", outputdir=outputdir, fileprefix=prefix, numbasis=numbasis, annuli=annulus_bounds, subsections=subsections, padding=padding, movement=movement)
        print('Done constructing forward model! You are ready to MCMC.')
