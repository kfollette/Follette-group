import os
import re
import subprocess
import glob
import astropy.io.fits as fits

from astropy import wcs
from astropy.modeling import models, fitting
import numpy as np
import scipy.ndimage as ndimage
import scipy.stats

import sys
from copy import copy
import configparser as ConfigParser

from scipy.interpolate import interp1d

class MAGAOData(object):
    
    """
    A sequence of P1640 Data. Each P1640Data object has the following fields and functions 
    Args:
        filepaths: list of filepaths to occulted files
        skipslices: a list of datacube slices to skip (supply index numbers e.g. [0,1,2,3])
        corefilepaths: a list of filepaths to core (i.e. unocculted) files, for contrast calc
        spot_directory: (None) path to the directory where the spot positions are stored. Defaults to P1640.ini val
    Attributes:
        input: Array of shape (N,y,x) for N images of shape (y,x)
        centers: Array of shape (N,2) for N centers in the format [x_cent, y_cent]
        filenums: Array of size N for the numerical index to map data to file that was passed in
        filenames: Array of size N for the actual filepath of the file that corresponds to the data
        PAs: Array of N for the parallactic angle rotation of the target (used for ADI) [in degrees]
        wvs: Array of N wavelengths of the images (used for SDI) [in microns]. For polarization data, defaults to "None"
        wcs: Array of N wcs astormetry headers for each image.
        IWA: a floating point scalar (not array). Specifies to inner working angle in pixels
        output: Array of shape (b, len(files), len(uniq_wvs), y, x) where b is the number of different KL basis cutoffs
        spot_flux: Array of N of average satellite spot flux for each frame
        contrast_scaling: Flux calibration factors (multiply by image to "calibrate" flux)
        flux_units: units of output data [DN, contrast]
        prihdrs: not used for P1640, set to None
        exthdrs: Array of N P1640 headers (these are written by the P1640 cube extraction pipeline)
    Methods:
        readdata(): reread in the data
        savedata(): save a specified data in the P1640 datacube format (in the 1st extension header)
        calibrate_output(): calibrates flux of self.output
    """

    #I'm marking things that I'm not sure if we need with a "#!"

    ##########################
   ### Class Initialization ###
    ##########################
    #Some static variables to define the MAGAO instrument
    centralwave = {} #in microns
    fpm_diam = {} #in pixels
    flux_zeropt = {}
    spot_ratio = {} #w.r.t. central star
    lenslet_scale = 1.0 #arcseconds per pixel (pixel scale)
    ifs_rotation = 0.0 #degrees CCW from +x axis to zenith
    
    observatory_latitude = 0.0

    #read in MAGAO configuration file and set these static variables
    package_directory = os.path.dirname(os.path.abspath(__file__))
    configfile = package_directory + "/" + "MAGAO.ini"
    config = ConfigParser.ConfigParser()
    try:
        config.read(configfile)
        #get pixel scale
        lenslet_scale = float(config.get("instrument", "ifs_lenslet_scale")) #!
        #get IFS rotation
        ifs_rotation = float(config.get("instrument", "ifs_rotation"))
        bands = ['HA', 'CONT']
        for band in bands:
            centralwave[band] = float(config.get("instrument", "cen_wave_{0}".format(band)))
            fpm_diam[band] = float(config.get("instrument", "fpm_diam_{0}".format(band))) #!
            flux_zeropt[band] = float(config.get("instrument", "zero_pt_flux_{0}".format(band))) #!
        observatory_latitude = float(config.get("observatory", "observatory_lat"))
    except ConfigParser.Error as e:
        print("Error reading MAGAO configuration file: {0}".format(e.message))
        raise e
    
    #########################
   ###    Constructors     ###
    #########################
    def __init__(self, filepaths=None):
        """
        Initialization code for MAGAOData
        """
        super(MAGAOData, self).__init__()
        self._output = None
        if filepaths is None:
            self._input = None
            self._centers = None
            self._filenums = None
            self._filenames = None
            self._PAs = None
            self._wvs = None
            self._wcs = None
            self._IWA = None
            self._OWA = None
            self.spot_flux = None #!
            self.star_flux = None
            self.contrast_scaling = None
            self.prihdrs = None
            self.exthdrs = None
        else:
            self.readdata(filepaths)
    
    ##############################
   ### Instance Required Fields ###
    ##############################
    @property
    def input(self):
        return self._input
    @input.setter
    def input(self, newval):
        self._input = newval
    
    @property
    def centers(self):
        return self._centers
    @centers.setter
    def centers(self, newval):
        self._centers = newval

    @property
    def PAs(self):
        return self._PAs
    @PAs.setter
    def PAs(self, newval):
        self._PAs = newval
    
    @property
    def wvs(self):
        return self._wvs
    @wvs.setter
    def wvs(self, newval):
        self._wvs = newval
    
    @property
    def wcs(self):
        return self._wcs
    @wcs.setter
    def wcs(self, newval):
        self._wcs = newval

    @property
    def IWA(self):
        return self._IWA
    @IWA.setter
    def IWA(self, newval):
        self._IWA = newval

    @property
    def OWA(self):
        return self._OWA
    @OWA.setter
    def OWA(self, newval):
        self._OWA = newval
    
    @property
    def output(self):
        return self._output
    @output.setter
    def output(self, newval):
        self._output = newval

    ###################
   ###    Methods    ###
    ###################
        
    def readdata(self, filepaths):
        """
        Method to open and read a list of MAGAO data
        """
        if isinstance(filepaths, str):
            filepaths = [filepaths]

        data = []
        filenums = []
        filenames = []
        rot_angles = []
        wvs = []
        centers = []
        wcs_hdrs = []
        star_fluxes = []
        spot_fluxes = [] #!
        prihdrs = []
        
        runningSum = 0
        for index, filepath in enumerate(filepaths):
            cube, center, pa, wv, astr_hdrs, filt_band, fpm_band, ppm_band, star_flux, spot_flux, prihdr, exthdr = _magao_process_file(filepath, index)
            
            runningSum = runningSum + 1
            #print("CUBE[0][0]: " + str(cube[0][0][0]))
            data.append(cube)
            centers.append(center)
            star_fluxes.append(star_flux)
            spot_fluxes.append(spot_flux) #!
            rot_angles.append(pa)
            wvs.append(wv)
            filenums.append(np.ones(pa.shape[0]) * index)
            wcs_hdrs.append(astr_hdrs) #!
            prihdrs.append(prihdr)
            filenames.append([filepath for i in range(pa.shape[0])])
                        
            
        #FILENUMS IS 1D LIST WITH LENGTH 68 AFTER FOR LOOP
        data = np.array(data)
        dims = data.shape
        #data = data.reshape([dims[0] * dims[1], dims[2], dims[3]])
        dims2 = data.shape
        #DATA HAS SIZE (68, 450, 450)
        filenums = np.array(filenums)
        filenums.reshape([dims[0]])
        #filenums = np.array(filenums)
        filenames = np.array(filenames).reshape([dims[0]])
        #filenames = np.array(filenames)
        rot_angles = np.array(rot_angles).reshape([dims[0]])
        #rot_angles = np.array(rot_angles)
        wvs = np.array(wvs).reshape([dims[0]])
        #wvs = np.array(wvs)
        #wcs_hdrs = np.array(wcs_hdrs).reshape([dims[0] * dims[1]])
        wcs_hdrs = np.array(wcs_hdrs)
        #centers = np.array(centers).reshape([dims[0] * dims[1], 2])
        dsize = dims[0]
        centers = np.zeros((dsize,2))
        for y in range(dsize):
            for x in range(2):
                centers[y][x] = dims[1]/2
                #centers[y][x] = 224.5
        #centers = np.array(centers)
        #star_fluxes = np.array(star_fluxes).reshape([dims[0] * dims[1]])
        star_fluxes = np.array(star_fluxes)
        #spot_fluxes = np.array(spot_fluxes).reshape([dims[0] * dims[1]]) #!
        spot_fluxes = np.array(spot_fluxes)

        self._input = data
        self._centers = centers
        self._filenums = filenums
        self._filenames = filenames
        self._PAs = rot_angles
        self._wvs = wvs
        self._wcs = None #wvs_hdrs
        self.spot_flux = spot_fluxes
        #self.IWA = MAGAOData.fpm_diam[fpm_band] / 2.0 #!
        #self.IWA = np.ones((68))
        with open('iwa.txt') as f:
            line=f.readline()
        self.IWA = int(line)
        #self.IWA = 1
        self.OWA = 450
        self.star_flux = star_fluxes
        self.contrast_scaling = 1./star_fluxes
        self.prihdrs = prihdrs

    def calibrate_output(self, img, spectral=False, units="contrast"):
        """
        Calibrates the flux of the output of PSF subtracted data.

        Assumes the broadband flux calibration is just multiplication by a single scalar number whereas spectral
        datacubes may have a separate calibration value for each wavelength

        Args:
            img: unclaibrated image.
                 If spectral is not set, this can either be a 2-D or 3-D broadband image
                 where the last two dimensions are [y,x]
                 If specetral is True, this is a 3-D spectral cube with shape [wv,y,x]
            spectral: if True, this is a spectral datacube. Otherwise, it is a broadband image.
            units: currently only support "contrast" w.r.t central star

        Return:
            img: calibrated image of the same shape (this is the same object as the input!!!)
        """
        if units == "contrast":
            if spectral:
                # spectral cube, each slice needs it's own calibration
                numwvs = img.shape[0]
                img *= self.contrast_scaling[:numwvs, None, None]
            else:
                # broadband image
                img *= np.nanmean(self.contrast_scaling)

        return img
        
    def savedata(self, filepath, data, klipparams = None, filetype = None, zaxis = None, center=None, astr_hdr=None,
                 fakePlparams = None,):
        """
        Save data in a GPI-like fashion. Aka, data and header are in the first extension header
        
        Inputs:
        filepath: path to file to output
        data: 2D or 3D data to save
        klipparams: a string of klip parameters
        filetype: filetype of the object (e.g. "KL Mode Cube", "PSF Subtracted Spectral Cube")
        zaxis: a list of values for the zaxis of the datacub (for KL mode cubes currently)
        astr_hdr: wcs astrometry header (None for NIRC2)
        center: center of the image to be saved in the header as the keywords PSFCENTX and PSFCENTY in pixels.
        The first pixel has coordinates (0,0)
        fakePlparams: fake planet params
        
        """
        hdulist = fits.HDUList()
        hdulist.append(fits.PrimaryHDU(header=self.prihdrs[0]))
        hdulist.append(fits.ImageHDU(data=data, name="Sci"))
        
        # save all the files we used in the reduction
        # we'll assume you used all the input files
        # remove duplicates from list
        #print("filenames = " + self._filenames)
        filenames = np.unique(self._filenames)
        nfiles = np.size(filenames)
        hdulist[0].header["DRPNFILE"] = nfiles
        for i, thispath in enumerate(filenames):
            thispath = thispath.replace("\\", '/')
            splited = thispath.split("/")
            fname = splited[-1]
#            matches = re.search('S20[0-9]{6}[SE][0-9]{4}', fname)
            filename = fname#matches.group(0)
            hdulist[0].header["FILE_{0}".format(i)] = filename

        # write out psf subtraction parameters
        # get pyKLIP revision number
        pykliproot = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        # the universal_newline argument is just so python3 returns a string instead of bytes
        # this will probably come to bite me later
        try:
            pyklipver = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=pykliproot, universal_newlines=True).strip()
        except:
            pyklipver = "unknown"
        hdulist[0].header['PSFSUB'] = "pyKLIP"
        hdulist[0].header.add_history("Reduced with pyKLIP using commit {0}".format(pyklipver))
        #if self.creator is None:
        #    hdulist[0].header['CREATOR'] = "pyKLIP-{0}".format(pyklipver)
        #else:
        #    hdulist[0].header['CREATOR'] = self.creator
        #    hdulist[0].header.add_history("Reduced by {0}".self.creator)

        # store commit number for pyklip
        hdulist[0].header['pyklipv'] = pyklipver

        if klipparams is not None:
            hdulist[0].header['PSFPARAM'] = klipparams
            hdulist[0].header.add_history("pyKLIP reduction with parameters {0}".format(klipparams))

        if fakePlparams is not None:
            hdulist[0].header['FAKPLPAR'] = fakePlparams
            hdulist[0].header.add_history("pyKLIP reduction with fake planet injection parameters {0}".format(fakePlparams))

        if filetype is not None:
            hdulist[0].header['FILETYPE'] = filetype

        if zaxis is not None:
            #Writing a KL mode Cube
            if "KL Mode" in filetype:
                hdulist[0].header['CTYPE3'] = 'KLMODES'
                #write them individually
                for i, klmode in enumerate(zaxis):
                    hdulist[0].header['KLMODE{0}'.format(i)] = klmode

        #use the dataset astr hdr if none was passed in
        #if astr_hdr is None:
        #    print self.wcs[0]
        #    astr_hdr = self.wcs[0]
        if astr_hdr is not None:
            #update astro header
            #I don't have a better way doing this so we'll just inject all the values by hand
            astroheader = astr_hdr.to_header()
            exthdr = hdulist[0].header
            exthdr['PC1_1'] = astroheader['PC1_1']
            exthdr['PC2_2'] = astroheader['PC2_2']
            try:
                exthdr['PC1_2'] = astroheader['PC1_2']
                exthdr['PC2_1'] = astroheader['PC2_1']
            except KeyError:
                exthdr['PC1_2'] = 0.0
                exthdr['PC2_1'] = 0.0
            #remove CD values as those are confusing
            exthdr.remove('CD1_1')
            exthdr.remove('CD1_2')
            exthdr.remove('CD2_1')
            exthdr.remove('CD2_2')
            exthdr['CDELT1'] = 1
            exthdr['CDELT2'] = 1

        #use the dataset center if none was passed in
        if center is None:
            center = self.centers[0]
        if center is not None:
            hdulist[0].header.update({'PSFCENTX':center[0],'PSFCENTY':center[1]})
            hdulist[0].header.update({'CRPIX1':center[0],'CRPIX2':center[1]})
            hdulist[0].header.add_history("Image recentered to {0}".format(str(center)))

        hdulist.writeto(filepath, clobber=True)
        hdulist.close()

        
def _magao_process_file(filepath, filetype):
    #filetype == 0 --> HA
    #filetype == 1 --> CONT
    #hdulist = fits.open(filepath)
    #header = hdulist[0].header
    #angle = float(header['ROTOFF'])
    #angle = 90+angle
    #angles = [angle]
    #angles = np.array(angles)
    #hdulist = fits.open(os.getcwd()+"/../HD142527/HD142527/8Apr14/MERGED_long_sets/rotoff_preproc.fits")
    #rotangles = hdulist[0].data
    #rotangles = np.array(rotangles)
    #rotangles = np.zeros((0))
    #print("ROTANGLES = " + str(angle))
    #hdulist.close()
    #hdulist = fits.open(filepath)
    try:
        hdulist = fits.open(filepath)
        header = hdulist[0].header
        angle = float(header['ROTOFF'])
        angle = 90+angle
        angles = [angle]
        angles = np.array(angles)
        cube = hdulist[0].data
        exthdr = None
        prihdr = hdulist[0].header
        
        if filetype == 0:
            filt_band = "H-Alpha"
        else:
            filt_band = "Continuum"
            
        fpm_band = filt_band
        ppm_band = None
            
        wvs = [1.0]
        #center = [[224.5,224.5]]
        datasize = cube.shape[1]
        center = [[datasize/2, datasize/2]]
        
        dims = cube.shape
        x, y = np.meshgrid(np.arange(dims[1], dtype=np.float32), np.arange(dims[0], dtype=np.float32))
        nx = center[0][0] - (x - center[0][0])
        #print("nx is " + str(nx))
        minval = np.min([np.nanmin(cube), 0.0])
        #flipped_cube = ndimage.map_coordinates(np.copy(cube), [y, nx], cval=minval * 5.0)
            
        #star_flux = calc_starflux(flipped_cube, center) #WRITE THIS FUNCTION
        star_flux = [[1]]
        #print("flipped_cube shape is " + str(flipped_cube.shape))
        #cube = flipped_cube.reshape([1, flipped_cube.shape[0], flipped_cube.shape[1]])
        cube.reshape([1, cube.shape[0], cube.shape[1]])
        parang = angles
        astr_hdrs = np.repeat(None, 1)
        spot_fluxes = [[1]] #!

    finally:
        hdulist.close()
        #print("Closing file")
        #fits.close(filepath)
        
    return cube, center, parang, wvs, astr_hdrs, filt_band, fpm_band, ppm_band, star_flux, spot_fluxes, prihdr, exthdr


def calc_starflux(cube, center):
    dims = cube.shape
    y, x = np.meshgrid(np.arange(dims[0]), np.arange(dims[1]))
    g_init = models.Gaussian2D(cube.max(), x_mean=center[0][0], y_mean=center[0][1], x_stddev=5, y_stddev=5, fixed={'x_mean':True,'y_mean':True,'theta':True})
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, y, x, cube)
    return [[g.amplitude]]
    
    
