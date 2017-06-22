__author__ = 'JB'

import os
from copy import copy
from glob import glob
from sys import stdout
import multiprocessing as mp
import itertools

import astropy.io.fits as pyfits
import numpy as np

from pyklip.kpp.utils.kppSuperClass import KPPSuperClass
from pyklip.instruments import GPI
import pyklip.kpp.utils.mathfunc as kppmath
import pyklip.spectra_management as spec


class ShapeOrMF(KPPSuperClass):
    """
    Class calculating the shape or matched filter metric for a GPI 2D image or a 3D cube.


    => The matched filter is a metric allowing one to pull out a known signal from noisy (not correlated) data.
    It is basically a dot product (meaning projection) of the template PSF cube with a stamp cube.

    => The shape metric is a normalized matched filter from which the flux component has been removed.
    It is a sort of pattern recognition.
    The shape value is a dot product normalized by the norm of the vectors

    It is meant to work on GPI campaign data with some automatic file search if inputs are missing.
    For kernels based on satellite spots it is assumed that inputDir will be an epoch folder from the GPI campaign.
    If it is not the case you should provide the PSF manually.

    More generally there is little guarantee of this function to work on non GPI data although this should be improved
    in the future.

    /!\ Caution: Most of the features need testing. Please contact jruffio@stanford.edu if there is a bug.

    Inherit from the Super class KPPSuperClass
    """
    def __init__(self,filename,metric_type,kernel_type,
                 inputDir = None,
                 outputDir = None,
                 PSF_cube_filename = None,
                 mute=None,
                 N_threads=None,
                 overwrite=False,
                 kernel_width = None,
                 sky_aper_radius = None,
                 label = None,
                 add2prefix = None,
                 keepPrefix = None,
                 SpT_file_csv=None):
        """
        Define the general parameters of the metric.

        The idea is that the directories and eventual spectra will be defined later by calling the
        initialize() function. Furthermore initialize() is where files can be read.
        It is done that way to ease the analysis of an entire campaign. ie one object is defined
        for a given metric and then a function calls initialize several time for each target's directory and each
        spectra.


        :param filename: Filename of the file on which to calculate the metric. It should be the complete path unless
                        inputDir is defined.
                        It can include wild characters. The file will be selected using the first output of glob.glob().
        :param metric_type: Either "shape" or "MF".
                        Shape is a matched filter which is also normalized wrt the image stamp.
        :param kernel_type: String defining type of model to be used for the cross correlation:
                    - "PSF": (Default) Use the PSF cube saved as fits with filename PSF_cube_filename.
                            There is no need to define PSF_cube_filename in the case of GPI campaign epochs folders.
                            Indeed it will try to find a PSF in inputDir (make sure there is only one PSF in there)
                            and if it can't find it it will try to generate one based on the raw spectral cubes in
                            inputDir. It uses the cubes listed in the extension header of pyklip images.
                            The PSF is flattened if the data is 2D.
                    - "hat": Define the kernel as a simple aperture photometry with radius kernel_width.
                            Default radius is 2 pixels.
                    - "Gaussian": define the kernel as a symmetric 2D gaussian with width (ie standard deviation) equal
                            to kernel_width. Default value of the width is 1.25.
        :param inputDir: If defined it allows filename to not include the whole path and just the filename.
                        Files will be read from inputDir.
                        Note tat inputDir might be redefined using initialize at any point.
                        If inputDir is None then filename is assumed to have the absolute path.
        :param outputDir: Directory where to create the folder containing the outputs.
                        Note tat inputDir might be redefined using initialize at any point.
                        If outputDir is None:
                            If inputDir is defined: outputDir = inputDir+os.path.sep+"planet_detec_"
        :param PSF_cube_filename: Filename filter for the PSF in inputdir.
                        If inputDir is not defined it should be an absolute path.
                        Default value is "*-original_radial_PSF_cube.fits".
                        Useful only if kernel_type = "PSF"
                        If a PSF cube is not explicitly given and one is read automatically it assumes there is only
                        one PSF cube in this folder.
        :param folderName: Name of the folder containing the outputs. It will be located in outputDir.
                        Default folder name is "default_out" if not spectrum is defined or the name of the spectrum
                        otherwise.
        :param mute: If True prevent printed log outputs.
        :param N_threads: Number of threads to be used for the metrics calculations.
                        If None use mp.cpu_count().
                        If -1 do it sequentially.
        :param overwrite_metric: If True check_existence() will always return False.
        :param kernel_width: Define the width of the Kernel depending on kernel_type. See kernel_width.
        :param sky_aper_radius: Radius of the mask applied on the stamps to calculated the background value.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        :param keepPrefix: Keep the prefix of the input file instead of using the default:
            self.prefix = self.star_name+"_"+self.compact_date+"_"+self.filter+self.add2prefix
        """
        # allocate super class
        super(ShapeOrMF, self).__init__(filename,
                                     inputDir = inputDir,
                                     outputDir = outputDir,
                                     folderName = None,
                                     mute=mute,
                                     N_threads=N_threads,
                                     label=label,
                                     overwrite=overwrite)

        if metric_type != "shape" and metric_type != "MF":
            raise Exception("Bad metric_type. Should be either 'shape' or 'PSF'.")
        else:
            self.metric_type = metric_type

        # Defalt value of kernel_type is "PSF"
        if kernel_type == None:
            self.kernel_type = "PSF"
        else:
            self.kernel_type = kernel_type

        # The default value is defined later
        self.kernel_width = kernel_width


        # Radius of the mask applied on the stamps to calculated the background value.
        if sky_aper_radius is None:
            self.sky_aper_radius = 2.5
        else:
            self.sky_aper_radius = sky_aper_radius

        if self.kernel_type == "PSF":
            self.PSF_cube_filename = PSF_cube_filename

        if add2prefix is not None:
            self.add2prefix = "_"+add2prefix
        else:
            self.add2prefix = ""

        if keepPrefix is not None:
            self.keepPrefix = keepPrefix
        else:
            self.keepPrefix = False

        self.SpT_file_csv=SpT_file_csv


    def spectrum_iter_available(self):
        """
        Indicates wether or not the class is equipped for iterating over spectra. That depends on the number of
        dimensions of the input file. If it is 2D there is no spectrum template so it can't and shouldn't iterate over
        spectra however if the input file is a cube then a spectrum is required.

        In order to iterate over spectra the function new_init_spectrum() can be called.
        spectrum_iter_available is a utility function for campaign data processing.

        :return: True if the file read is 3D.
                False if it is 2D.
                None if no file has been read and there is no way to know.
        """

        if hasattr(self,"is3D"):
            return self.is3D
        else:
            return None

    def init_new_spectrum(self,spectrum):
        """
        Function allowing the reinitialization of the class with a new spectrum without reinitializing what doesn't need
         to be.

        :param spectrum: spectrum path relative to pykliproot + os.path.sep + "spectra" with pykliproot the directory in
                        which pyklip is installed. It that case it should be a spectrum from Mark Marley.
                        Instead of a path it can be a simple ndarray with the right dimension.
                        Or by default it is a completely flat spectrum.

        :return: None
        """
        # Define the output Foldername
        if isinstance(spectrum, str):
            if spectrum != "":
                pykliproot = os.path.dirname(os.path.realpath(spec.__file__))
                self.spectrum_filename = os.path.abspath(glob(os.path.join(pykliproot,"spectra","*",spectrum+".flx"))[0])
                spectrum_name = self.spectrum_filename.split(os.path.sep)
                self.spectrum_name = spectrum_name[len(spectrum_name)-1].split(".")[0]
            else:
                self.spectrum_name = "satSpotSpec"

            # Do the best it can with the spectral information given in inputs.
            if spectrum != "":
                # spectrum_filename is not empty it is assumed to be a valid path.
                if not self.mute:
                    print("Spectrum model: "+self.spectrum_filename)
                # Interpolate the spectrum of the planet based on the given filename
                wv,planet_sp = spec.get_planet_spectrum(self.spectrum_filename,self.filter)

                if hasattr(self, 'sat_spot_spec') \
                        and (self.star_type is not None or self.star_temperature is not None):
                    # Interpolate a spectrum of the star based on its spectral type/temperature
                    wv,star_sp = spec.get_star_spectrum(self.filter,self.star_type,self.star_temperature)
                    # Correct the ideal spectrum given in spectrum_filename for atmospheric and instrumental absorption.
                    self.spectrum_vec = (self.sat_spot_spec/star_sp)*planet_sp
                else:
                    # If the sat spot spectrum or the spectral type of the star is missing it won't correct anything
                    # and keep the ideal spectrum as it is.
                    if not self.mute:
                        print("No star spec or sat spot spec so using sole planet spectrum.")
                    self.spectrum_vec = copy(planet_sp)
            else:
                # If spectrum_filename is an empty string the function takes the sat spot spectrum by default.
                if hasattr(self, 'sat_spot_spec'):
                    if not self.mute:
                        print("Default sat spot spectrum will be used.")
                    self.spectrum_vec = copy(self.sat_spot_spec)
                else:
                    # If the sat spot spectrum is also not given then it just take the band filter spectrum.
                    if not self.mute:
                        print("Using gpi filter "+self.filter+" spectrum. Could find neither sat spot spectrum nor planet spectrum.")
                    wv,self.spectrum_vec = spec.get_gpi_filter(self.filter)

        elif isinstance(spectrum, np.ndarray):
            self.spectrum_vec = spectrum
            self.spectrum_name = "custom"
        else:
            if not self.mute:
                print("Spectrum is not or badly defined so taking flat spectrum")
            self.spectrum_vec = np.ones(self.nl)
            self.spectrum_name = "flat"


        self.folderName = self.spectrum_name+os.path.sep


        for k in range(self.nl):
            self.PSF_cube[k,:,:] *= self.spectrum_vec[k]
        # normalize spectrum with norm 2.
        self.spectrum_vec = self.spectrum_vec / np.sqrt(np.nansum(self.spectrum_vec**2))
        # normalize PSF with norm 2.
        self.PSF_cube = self.PSF_cube/np.sqrt(np.sum(self.PSF_cube**2))

        return None

    def initialize(self, inputDir = None,
                         outputDir = None,
                         spectrum = None,
                         folderName = None,
                         PSF_cube_filename = None,
                         prihdr = None,
                         exthdr = None,
                         star_type = None,
                         star_temperature = None,
                         compact_date = None,
                         label=None):
        """

        Initialize the non general inputs that are needed for the metric calculation and load required files.


        - It reads the input fits file to be reduced.
        - Build or load the PSF if needed
        - load the spectrum if needed (it can be redefined later using self.init_new_spectrum())
        - read headers if they exist


        :param inputDir: If defined it allows filename to not include the whole path and just the filename.
                        Files will be read from inputDir.
                        Note tat inputDir might be redefined using initialize at any point.
                        If inputDir is None then filename is assumed to have the absolute path.
        :param outputDir: Directory where to create the folder containing the outputs.
                        Note tat inputDir might be redefined using initialize at any point.
                        If outputDir is None:
                            If inputDir is defined: outputDir = inputDir+os.path.sep+"planet_detec_"
        :param spectrum: spectrum path relative to pykliproot + os.path.sep + "spectra" with pykliproot the directory in
                        which pyklip is installed. It that case it should be a spectrum from Mark Marley.
                        Instead of a path it can be a simple ndarray with the right dimension.
                        Or by default it is a completely flat spectrum.
        :param folderName: Name of the folder containing the outputs. It will be located in outputDir.
                        Default folder name is "default_out" if not spectrum is defined or the name of the spectrum
                        otherwise.
        :param PSF_cube_filename: Filename filter for the PSF in inputdir.
                        If inputDir is not defined it should be an absolute path.
                        Default value is "*-original_radial_PSF_cube.fits".
                        Useful only if kernel_type = "PSF"
                        If a PSF cube is not explicitly given and one is read automatically it assumes there is only
                        one PSF cube in this folder.
        :param prihdr: User defined primary fits headers in case the file read has none.
        :param exthdr: User defined extension fits headers in case the file read has none.
        :param star_type: String containing the spectral type of the star. 'A5','F4',... Assume type V star. It is ignored
                        of temperature is defined.
        :param star_temperature: Temperature of the star. Overwrite star_type if defined.
        :param compact_date: Define the compact date to be used in the output filenames.
                            If a PSF has to be measured from the satellite spots it will define itself.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        :return:
        """
        if not self.mute:
            print("~~ INITializing "+self.__class__.__name__+" ~~")

        # The super class already read the fits file
        init_out = super(ShapeOrMF, self).initialize(inputDir = inputDir,
                                         outputDir = outputDir,
                                         folderName = folderName,
                                         label=label)

        if compact_date is None:
            self.compact_date = "noDate"
        else:
            self.compact_date = compact_date

        if prihdr is not None:
            self.prihdr = prihdr
            if not self.mute:
                print("Overwriting primary header with user defined prihdr")
        if exthdr is not None:
            self.exthdr = exthdr
            if not self.mute:
                print("Overwriting extension header with user defined exthdr")

        ## Read few fits headers
        # Get center of the image (star position)
        try:
            # Retrieve the center of the image from the fits headers.
            self.center = [self.exthdr['PSFCENTX'], self.exthdr['PSFCENTY']]
        except:
            # If the keywords could not be found the center is defined as the middle of the image
            if not self.mute:
                print("Couldn't find PSFCENTX and PSFCENTY keywords.")
            self.center = [(self.nx-1)/2,(self.ny-1)/2]

        if self.label == "CADI":
            self.center = [140,140]

        # Get the filter of the image
        try:
            # Retrieve the filter used from the fits headers.
            self.filter = self.prihdr['IFSFILT'].split('_')[1]
        except:
            # If the keywords could not be found assume that the filter is H...
            if not self.mute:
                print("Couldn't find IFSFILT keyword. If relevant assuming H.")
            self.filter = "H"

        # Get current star name
        try:
            # OBJECT: keyword in the primary header with the name of the star.
            self.star_name = self.prihdr['OBJECT'].strip().replace (" ", "_")
        except:
            # If the object name could nto be found cal lit unknown_object
            self.star_name = "UNKNOWN_OBJECT"

        # Check if the input file is 2D or 3D
        if hasattr(self, 'nl'): # If the file is a 3D cube
            self.is3D = True
        else:
            self.is3D = False

        # If the Kernel is a PSF build it here
        # if self.kernel_type == "PSF":
        if 1: # Reading the PSF all the time to get access to the sat spot spectrum
            self.PSF_cube_filename = PSF_cube_filename

            if self.PSF_cube_filename is None:
                self.PSF_cube_filename = "*-original_radial_PSF_cube.fits"
                if not self.mute:
                    print("Using default filename for PSF cube: "+self.PSF_cube_filename)

            # Define the actual filename path
            if self.inputDir is None:
                try:
                    self.PSF_cube_path = os.path.abspath(glob(self.PSF_cube_filename)[0])
                except:
                    raise Exception("File "+self.PSF_cube_filename+"doesn't exist.")
            else:
                try:
                    self.PSF_cube_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.PSF_cube_filename)[0])
                except:
                    if not self.mute:
                        print("File "+self.inputDir+os.path.sep+self.PSF_cube_filename+"doesn't exist.")
                        print("Try to generate the PSF from the spectral cube in " + self.inputDir)
                        print("Work only for GPI campaign data.")
                    try:
                        cube_filename_list = []
                        file_index = 0
                        EOList = False
                        while not EOList:
                            try:
                                curr_spdc_filename = self.prihdr["FILE_{0}".format(file_index)]
                                cube_filename_list.append(self.inputDir+os.path.sep+curr_spdc_filename.split(".")[0]+"_spdc_distorcorr.fits")
                            except:
                                EOList = True
                            file_index = file_index +1

                        self.compact_date = curr_spdc_filename.split("S")[1]
                        prefix = self.star_name+"_"+self.compact_date+"_"+self.filter

                        # Load all the data cube into a dataset object
                        dataset = GPI.GPIData(cube_filename_list)
                        if not self.mute:
                            print("Calculating the planet PSF from the satellite spots...")
                        # generate the PSF cube from the satellite spots
                        dataset.generate_psf_cube(20)
                        # Save the original PSF calculated from combining the sat spots
                        dataset.savedata(self.inputDir + os.path.sep + prefix+"-original_PSF_cube.fits", dataset.psfs,
                                                  astr_hdr=dataset.wcs[0], filetype="PSF Spec Cube")
                        # Calculate and save the rotationally invariant psf (ie smeared out/averaged).
                        # Generate a radially averaged PSF cube and save as
                        # self.inputDir + os.path.sep + prefix+"-original_radial_PSF_cube.fits"
                        dataset.get_radial_psf(save = self.inputDir + os.path.sep + prefix)
                        self.PSF_cube_path = self.inputDir + os.path.sep + prefix+"-original_radial_PSF_cube.fits"
                        if not self.mute:
                            print(self.inputDir + os.path.sep + prefix+"-original_radial_PSF_cube.fits"+"and"+\
                                  self.inputDir + os.path.sep + prefix+"-original_PSF_cube.fits"+ " have been saved.")
                    except:
                        raise Exception("File "+self.inputDir+os.path.sep+self.PSF_cube_filename+"doesn't exist and"+
                                        "couldn't generate one. inputDir is probably not a GPI campaign directory or "+\
                                        "the exthdr doesn't include the list of spdc.")

            if not self.mute:
                print("Loading PSF cube: "+self.PSF_cube_path)

            # Load the PSF cube if a file has been found or was just generated
            hdulist = pyfits.open(self.PSF_cube_path)
            self.PSF_cube = hdulist[1].data
            # Extract the spectrum of the satellite spot based of the input PSF
            self.sat_spot_spec = np.nanmax(self.PSF_cube,axis=(1,2))
            # Remove the spectral shape from the psf cube because it is dealt with independently
            for l_id in range(self.PSF_cube.shape[0]):
                self.PSF_cube[l_id,:,:] = self.PSF_cube[l_id,:,:]/self.sat_spot_spec[l_id]


        if self.kernel_type == "PSF" and not self.is3D:
            # flatten the PSF cube if the data is 2D
            self.PSF_flat =  np.nanmean(self.PSF_cube,axis=0)

        # Define the PSF as a gaussian
        if self.kernel_type == "gaussian":
            if self.kernel_width == None:
                self.kernel_width = 1.25
                if not self.mute:
                    print("Default width sigma = {0} used for the gaussian".format(self.kernel_width))

            if not self.mute:
                print("Generate gaussian PSF")
            # Build the grid for PSF stamp.
            ny_PSF = 20 # should be even
            nx_PSF = 20 # should be even
            x_PSF_grid, y_PSF_grid = np.meshgrid(np.arange(0,ny_PSF,1)-ny_PSF/2,np.arange(0,nx_PSF,1)-nx_PSF/2)

            PSF = kppmath.gauss2d(x_PSF_grid, y_PSF_grid,1.0,0.0,0.0,self.kernel_width,self.kernel_width)

            if self.is3D:
                # Duplicate the PSF to get a PSF cube.
                # Caution: No spectral widening implemented here
                self.PSF_cube = np.tile(PSF,(self.nl,1,1))
            else:
                self.PSF_flat = PSF

        # Define the PSF as an aperture or "hat" function
        if self.kernel_type == "hat":
            if self.kernel_width == None:
                self.kernel_width = 2.
                if not self.mute:
                    print("Default radius = {0} used for the hat function".format(self.kernel_width))

            # Build the grid for PSF stamp.
            ny_PSF = 20 # should be even
            nx_PSF = 20 # should be even
            x_PSF_grid, y_PSF_grid = np.meshgrid(np.arange(0,ny_PSF,1)-ny_PSF/2,np.arange(0,nx_PSF,1)-nx_PSF/2)
            # Use aperture for the cross correlation.
            # Calculate the corresponding hat function
            PSF = kppmath.hat(x_PSF_grid, y_PSF_grid, self.kernel_width)

            if self.is3D:
                # Duplicate the PSF to get a PSF cube.
                # Caution: No spectral widening implemented here
                self.PSF_cube = np.tile(PSF,(self.nl,1,1))
            else:
                self.PSF_flat = PSF

        if self.is3D:
            self.nl_PSF, self.ny_PSF, self.nx_PSF = self.PSF_cube.shape
        else:
            self.ny_PSF, self.nx_PSF = self.PSF_flat.shape

        if star_type is None and  star_temperature is None:
            self.star_type = spec.get_specType(self.star_name,self.SpT_file_csv)
            self.star_temperature = star_temperature
        else:
            self.star_type = star_type
            self.star_temperature = star_temperature

        # Load the spectrum here if the data is 3D
        if self.is3D:
            if isinstance(spectrum, str):
                if spectrum != "":
                    pykliproot = os.path.dirname(os.path.realpath(spec.__file__))
                    self.spectrum_filename = os.path.abspath(glob(os.path.join(pykliproot,"spectra","*",spectrum+".flx"))[0])
                    spectrum_name = self.spectrum_filename.split(os.path.sep)
                    self.spectrum_name = spectrum_name[len(spectrum_name)-1].split(".")[0]
                else:
                    self.spectrum_name = "satSpotSpec"

                # Do the best it can with the spectral information given in inputs.
                if spectrum != "":
                    # spectrum_filename is not empty it is assumed to be a valid path.
                    if not self.mute:
                        print("Spectrum model: "+self.spectrum_filename)
                    # Interpolate the spectrum of the planet based on the given filename
                    wv,planet_sp = spec.get_planet_spectrum(self.spectrum_filename,self.filter)

                    if hasattr(self, 'sat_spot_spec') \
                            and (self.star_type is not None or self.star_temperature is not None):
                        # Interpolate a spectrum of the star based on its spectral type/temperature
                        wv,star_sp = spec.get_star_spectrum(self.filter,self.star_type,self.star_temperature)
                        # Correct the ideal spectrum given in spectrum_filename for atmospheric and instrumental absorption.
                        self.spectrum_vec = (self.sat_spot_spec/star_sp)*planet_sp
                    else:
                        # If the sat spot spectrum or the spectral type of the star is missing it won't correct anything
                        # and keep the ideal spectrum as it is.
                        if not self.mute:
                            print("No star spec or sat spot spec so using sole planet spectrum.")
                        self.spectrum_vec = copy(planet_sp)
                else:
                    # If spectrum_filename is an empty string the function takes the sat spot spectrum by default.
                    if hasattr(self, 'sat_spot_spec'):
                        if not self.mute:
                            print("Default sat spot spectrum will be used.")
                        self.spectrum_vec = copy(self.sat_spot_spec)
                    else:
                        # If the sat spot spectrum is also not given then it just take the band filter spectrum.
                        if not self.mute:
                            print("Using gpi filter "+self.filter+" spectrum. Could find neither sat spot spectrum nor planet spectrum.")
                        wv,self.spectrum_vec = spec.get_gpi_filter(self.filter)

            elif isinstance(spectrum, np.ndarray):
                self.spectrum_vec = spectrum
                self.spectrum_name = "custom"
            else:
                if not self.mute:
                    print("Spectrum is not or badly defined so taking flat spectrum")
                self.spectrum_vec = np.ones(self.nl)
                self.spectrum_name = "flat"


            self.folderName = self.spectrum_name+os.path.sep

            for k in range(self.nl):
                self.PSF_cube[k,:,:] = self.PSF_cube[k,:,:]*self.spectrum_vec[k]
            # normalize spectrum with norm 2.
            self.spectrum_vec = self.spectrum_vec / np.sqrt(np.nansum(self.spectrum_vec**2))
            # normalize PSF with norm 2.
            self.PSF_cube = self.PSF_cube / np.sqrt(np.sum(self.PSF_cube**2))
        else: # if 2D data we still want to normalize the PSF
            # normalize PSF with norm 2.
            self.PSF_flat = self.PSF_flat / np.sqrt(np.sum(self.PSF_flat**2))

        # Define the suffix used when saving files
        if self.is3D:
            dim_suffix = "3D"
        else:
            dim_suffix = "2D"
        self.suffix = self.metric_type+dim_suffix+self.kernel_type

        if self.keepPrefix:
            file_ext_ind = os.path.basename(self.filename_path)[::-1].find(".")
            self.prefix = os.path.basename(self.filename_path)[:-(file_ext_ind+1)]
        else:
            self.prefix = self.star_name+"_"+self.compact_date+"_"+self.filter+self.add2prefix


        return init_out

    def check_existence(self):
        """
        Check if this metric has already been calculated for this file.

        If self.overwrite is True then the output will always be False

        :return: Boolean indicating the existence of the metric map.
        """

        file_exist = (len(glob(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')) >= 1)

        if file_exist and not self.mute:
            print("Output already exist.")

        return file_exist and not self.overwrite


    def calculate(self):
        """
        Calculate the metric map.

        :return: shape or matched filter map.
        """

        if not self.mute:
            print("~~ Calculating "+self.__class__.__name__+" with parameters " + self.suffix+" ~~")

        if self.is3D:
            flat_cube = np.nanmean(self.image,axis=0)
        else:
            flat_cube = self.image

        # Get the nans pixels of the flat_cube. We won't bother trying to calculate metrics for those.
        flat_cube_nans = np.where(np.isnan(flat_cube))

        # Remove the very edges of the image. We can't calculate a proper projection of an image stamp onto the PSF if we
        # are too close from the edges of the array.
        flat_cube_mask = np.ones((self.ny,self.nx))
        flat_cube_mask[flat_cube_nans] = np.nan
        flat_cube_noEdges_mask = copy(flat_cube_mask)
        # remove the edges if not already nans
        flat_cube_noEdges_mask[0:self.ny_PSF/2,:] = np.nan
        flat_cube_noEdges_mask[:,0:self.nx_PSF/2] = np.nan
        flat_cube_noEdges_mask[(self.ny-self.ny_PSF/2):self.ny,:] = np.nan
        flat_cube_noEdges_mask[:,(self.nx-self.nx_PSF/2):self.nx] = np.nan
        # Get the pixel coordinates corresponding to non nan pixels and not too close from the edges of the array.
        flat_cube_noNans_noEdges = np.where(np.isnan(flat_cube_noEdges_mask) == 0)

        if not self.mute:
            print("Calculating "+self.metric_type+" for "+self.filename_path)
        metric_map = -np.ones((self.ny,self.nx)) + np.nan

        # Calculate the criterion map.
        # For each pixel calculate the dot product of a stamp around it with the PSF.
        # We use the PSF cube to consider also the spectrum of the planet we are looking for.
        if not self.mute:
            print("Calculate the criterion map. It is done pixel per pixel so it might take a while...")
        stamp_PSF_x_grid, stamp_PSF_y_grid = np.meshgrid(np.arange(0,self.nx_PSF,1)-self.nx_PSF/2,
                                                         np.arange(0,self.ny_PSF,1)-self.ny_PSF/2)
        r_PSF_stamp = (stamp_PSF_x_grid)**2 +(stamp_PSF_y_grid)**2
        where_mask = np.where(r_PSF_stamp < (self.sky_aper_radius**2))
        if self.is3D:
            stamp_PSF_mask = np.ones((self.nl,self.ny_PSF,self.nx_PSF))
            stamp_PSF_mask[:,where_mask] = np.nan
        else:
            stamp_PSF_mask = np.ones((self.ny_PSF,self.nx_PSF))
            stamp_PSF_mask[where_mask] = np.nan


        if self.N_threads > 0:
            pool = mp.Pool(processes=self.N_threads)

            ## cut images in N_threads part
            # get the first and last index of each chunck
            N_pix = flat_cube_noNans_noEdges[0].size
            chunk_size = N_pix/self.N_threads
            N_chunks = N_pix/chunk_size

            # Get the chunks
            chunks_row_indices = []
            chunks_col_indices = []
            for k in range(N_chunks-1):
                chunks_row_indices.append(flat_cube_noNans_noEdges[0][(k*chunk_size):((k+1)*chunk_size)])
                chunks_col_indices.append(flat_cube_noNans_noEdges[1][(k*chunk_size):((k+1)*chunk_size)])
            chunks_row_indices.append(flat_cube_noNans_noEdges[0][((N_chunks-1)*chunk_size):N_pix])
            chunks_col_indices.append(flat_cube_noNans_noEdges[1][((N_chunks-1)*chunk_size):N_pix])

            if self.is3D:
                if self.metric_type == "shape":
                    func = calculate_shape3D_metric_star
                else:
                    func = calculate_MF3D_metric_star
                outputs_list = pool.map(func, itertools.izip(chunks_row_indices,
                                                           chunks_col_indices,
                                                           itertools.repeat(self.image),
                                                           itertools.repeat(self.PSF_cube),
                                                           itertools.repeat(stamp_PSF_mask)))
            else:
                if self.metric_type == "shape":
                    func = calculate_shape2D_metric_star
                else:
                    func = calculate_MF2D_metric_star

                # tpool_outputs = []
                # tpool_outputs += [pool.apply_async(calculate_shape2D_metric,
                #                                     args=(a))
                #                     for a in itertools.izip(chunks_row_indices,
                #                                            chunks_col_indices,
                #                                            itertools.repeat(self.image),
                #                                            itertools.repeat(self.PSF_flat),
                #                                            itertools.repeat(stamp_PSF_mask))]
                # while len(tpool_outputs) > 0:
                #     coucou = tpool_outputs.pop(0)
                #     coucou.wait()
                #     print(coucou.get())
                # print("coucou")

                outputs_list = pool.map(func, itertools.izip(chunks_row_indices,
                                                           chunks_col_indices,
                                                           itertools.repeat(self.image),
                                                           itertools.repeat(self.PSF_flat),
                                                           itertools.repeat(stamp_PSF_mask)))

            for row_indices,col_indices,out in zip(chunks_row_indices,chunks_col_indices,outputs_list):
                metric_map[(row_indices,col_indices)] = out
            pool.close()
        else:
            if self.is3D:
                if self.metric_type == "shape":
                    func = calculate_shape3D_metric
                else:
                    func = calculate_MF3D_metric
                metric_map[flat_cube_noNans_noEdges] = func(flat_cube_noNans_noEdges[0],
                                                           flat_cube_noNans_noEdges[1],
                                                           self.image,
                                                           self.PSF_cube,
                                                           stamp_PSF_mask)
            else:
                if self.metric_type == "shape":
                    func = calculate_shape2D_metric
                else:
                    func = calculate_MF2D_metric
                metric_map[flat_cube_noNans_noEdges] = func(flat_cube_noNans_noEdges[0],
                                                           flat_cube_noNans_noEdges[1],
                                                           self.image,
                                                           self.PSF_flat,
                                                           stamp_PSF_mask)
        self.metricMap = metric_map
        return self.metricMap


    def save(self):
        """
        Save the metric map as a fits file as
        self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits'

        It saves the metric parameters if self.prihdr and self.exthdr were already defined.

        :return: None
        """

        if not os.path.exists(self.outputDir+os.path.sep+self.folderName):
            os.makedirs(self.outputDir+os.path.sep+self.folderName)

        if hasattr(self,"prihdr") and hasattr(self,"exthdr"):
            # Save the parameters as fits keywords
            # MET##### stands for METric
            self.exthdr["MET_TYPE"] = self.metric_type

            self.exthdr["METFILEN"] = self.filename_path
            self.exthdr["METKERTY"] = self.kernel_type
            self.exthdr["METSKYRA"] = self.sky_aper_radius
            self.exthdr["METINDIR"] = self.inputDir
            self.exthdr["METOUTDI"] = self.outputDir
            self.exthdr["METFOLDN"] = self.folderName
            self.exthdr["METCDATE"] = self.compact_date
            self.exthdr["METPREFI"] = self.prefix
            self.exthdr["METSUFFI"] = self.suffix

            # This parameters are not always defined
            if hasattr(self,"spectrum_name"):
                self.exthdr["METSPECN"] = self.spectrum_name
            if hasattr(self,"spectrum_filename"):
                self.exthdr["METSPECF"] = self.spectrum_filename
            if hasattr(self,"PSF_cube_path"):
                self.exthdr["METPSFDI"] = self.PSF_cube_path
            if hasattr(self,"star_type"):
                self.exthdr["METSTTYP"] = self.star_type
            if hasattr(self,"star_temperature"):
                self.exthdr["METSTTEM"] = self.star_temperature
            if hasattr(self,"kernel_width"):
                self.exthdr["METKERWI"] = self.kernel_width

            if not self.mute:
                print("Saving: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')
            hdulist = pyfits.HDUList()
            hdulist.append(pyfits.PrimaryHDU(header=self.prihdr))
            hdulist.append(pyfits.ImageHDU(header=self.exthdr, data=self.metricMap, name=self.suffix))
            hdulist.writeto(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits', clobber=True)
        else:
            hdulist = pyfits.HDUList()
            hdulist.append(pyfits.ImageHDU(data=self.metricMap, name=self.suffix))
            hdulist.append(pyfits.ImageHDU(name=self.suffix))

            hdulist[1].header["MET_TYPE"] = self.metricMap

            hdulist[1].header["METFILEN"] = self.filename_path
            hdulist[1].header["METKERTY"] = self.kernel_type
            hdulist[1].header["METSKYRA"] = self.sky_aper_radius
            hdulist[1].header["METINDIR"] = self.inputDir
            hdulist[1].header["METOUTDI"] = self.outputDir
            hdulist[1].header["METFOLDN"] = self.folderName
            hdulist[1].header["METCDATE"] = self.compact_date
            hdulist[1].header["METPREFI"] = self.prefix
            hdulist[1].header["METSUFFI"] = self.suffix

            # This parameters are not always defined
            if hasattr(self,"spectrum_name"):
                hdulist[1].header["METSPECN"] = self.spectrum_name
            if hasattr(self,"spectrum_filename"):
                hdulist[1].header["METSPECF"] = self.spectrum_filename
            if hasattr(self,"PSF_cube_path"):
                hdulist[1].header["METPSFDI"] = self.PSF_cube_path
            if hasattr(self,"star_type"):
                hdulist[1].header["METSTTYP"] = self.star_type
            if hasattr(self,"star_temperature"):
                hdulist[1].header["METSTTEM"] = self.star_temperature
            if hasattr(self,"kernel_width"):
                hdulist[1].header["METKERWI"] = self.kernel_width

            if not self.mute:
                print("Saving: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')
            hdulist.writeto(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits', clobber=True)

        return None


    def load(self):
        """
        Load the metric map. One should check that it exist first using self.check_existence().

        Define the attribute self.metricMap.

        :return: self.metricMap
        """
        try:
            hdulist = pyfits.open(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')
            self.metricMap = hdulist[1].data
            hdulist.close()
        except:
            hdulist = pyfits.open(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')
            self.metricMap = hdulist[0].data
            hdulist.close()


        return self.metricMap

def calculate_shape2D_metric_star(params):
    """
    Convert `f([1,2])` to `f(1,2)` call.
    It allows one to call calculate_shape3D_metric() with a tuple of parameters.
    """
    return calculate_shape2D_metric(*params)

def calculate_shape2D_metric(row_indices,col_indices,cube,PSF_cube,stamp_PSF_mask, mute = True):
    '''
    Calculate the shape metric on the given 2D image for the pixels targeted by row_indices and col_indices.
    These lists of indices can basically be given from the numpy.where function following the example:
        import numpy as np
        row_indices,col_indices = np.where(np.finite(np.mean(cube,axis=0)))
    By truncating the given lists in small pieces it is then easy to parallelized.

    The shape metric is a normalized matched filter from which the flux component has been removed.
    It is a sort of pattern recognition.
    The shape value is a dot product normalized by the norm of the vectors

    :param row_indices: Row indices list of the pixels where to calculate the metric in cube.
                        Indices should be given from a 2d image.
    :param col_indices: Column indices list of the pixels where to calculate the metric in cube.
                        Indices should be given from a 2d image.
    :param cube: Cube from which one wants the metric map. PSF_cube should be norm-2 normalized.
                PSF_cube /= np.sqrt(np.sum(PSF_cube**2))
    :param PSF_flat: PSF template used for calculated the metric. If ny_PSF,nx_PSF = PSF_cube.shape ny_PSF and nx_PSF
                    are the spatial dimensions of the PSF_cube.
    :param stamp_PSF_mask: 2d mask of size (ny_PSF,nx_PSF) used to mask the central part of a stamp slice. It is used as
                        a type of a high pass filter. Before calculating the metric value of a stamp cube around a given
                        pixel the average value of the surroundings of each slice of that stamp cube will be removed.
                        The pixel used for calculating the average are the one equal to one in the mask.
    :param mute: If True prevent printed log outputs.
    :return: Vector of length row_indices.size with the value of the metric for the corresponding pixels.
    '''

    # Shape of the PSF cube
    ny_PSF,nx_PSF = PSF_cube.shape

    # Number of rows and columns to add around a given pixel in order to extract a stamp.
    row_m = np.floor(ny_PSF/2.0)    # row_minus
    row_p = np.ceil(ny_PSF/2.0)     # row_plus
    col_m = np.floor(nx_PSF/2.0)    # col_minus
    col_p = np.ceil(nx_PSF/2.0)     # col_plus
    # Number of pixels on which the metric has to be computed
    N_it = row_indices.size
    # Define an shape vector full of nans
    shape_map = np.zeros((N_it,)) + np.nan
    # Loop over all pixels (row_indices[id],col_indices[id])
    for id,k,l in zip(range(N_it),row_indices,col_indices):
        if not mute:
            # Print the progress of the function
            stdout.write("\r{0}/{1}".format(id,N_it))
            stdout.flush()

        # Extract stamp cube around the current pixel from the whoel cube
        stamp_cube = copy(cube[(k-row_m):(k+row_p), (l-col_m):(l+col_p)])
        # Remove average value of the surrounding pixels in each slice of the stamp cube
        stamp_cube -= np.nanmean(stamp_cube*stamp_PSF_mask)
        # Dot product of the PSF with stamp cube.
        ampl = np.nansum(PSF_cube*stamp_cube)
        # Normalize the dot product square by the squared norm-2 of the stamp cube.
        # Because we keep the sign shape value is then found in [-1.,1.]
        try:
            shape_map[id] = np.sign(ampl)*ampl**2/np.nansum(stamp_cube**2)
        except:
            # In case ones divide by zero...
            shape_map[id] =  np.nan

    # The shape value here can be seen as a cosine square as it is a normalized squared dot product.
    # Taking the square root to make it a simple cosine.
    return np.sign(shape_map)*np.sqrt(abs(shape_map))


def calculate_shape3D_metric_star(params):
    """
    Convert `f([1,2])` to `f(1,2)` call.
    It allows one to call calculate_shape3D_metric() with a tuple of parameters.
    """
    return calculate_shape3D_metric(*params)

def calculate_shape3D_metric(row_indices,col_indices,cube,PSF_cube,stamp_PSF_mask, mute = True):
    '''
    Calculate the shape metric on the given datacube for the pixels targeted by row_indices and col_indices.
    These lists of indices can basically be given from the numpy.where function following the example:
        import numpy as np
        row_indices,col_indices = np.where(np.finite(np.mean(cube,axis=0)))
    By truncating the given lists in small pieces it is then easy to parallelized.

    The shape metric is a normalized matched filter from which the flux component has been removed.
    It is a sort of pattern recognition.
    The shape value is a dot product normalized by the norm of the vectors

    :param row_indices: Row indices list of the pixels where to calculate the metric in cube.
                        Indices should be given from a 2d image.
    :param col_indices: Column indices list of the pixels where to calculate the metric in cube.
                        Indices should be given from a 2d image.
    :param cube: Cube from which one wants the metric map. PSF_cube should be norm-2 normalized.
                PSF_cube /= np.sqrt(np.sum(PSF_cube**2))
    :param PSF_cube: PSF_cube template used for calculated the metric. If nl,ny_PSF,nx_PSF = PSF_cube.shape, nl is the
                     number of wavelength samples, ny_PSF and nx_PSF are the spatial dimensions of the PSF_cube.
    :param stamp_PSF_mask: 2d mask of size (ny_PSF,nx_PSF) used to mask the central part of a stamp slice. It is used as
                        a type of a high pass filter. Before calculating the metric value of a stamp cube around a given
                        pixel the average value of the surroundings of each slice of that stamp cube will be removed.
                        The pixel used for calculating the average are the one equal to one in the mask.
    :param mute: If True prevent printed log outputs.
    :return: Vector of length row_indices.size with the value of the metric for the corresponding pixels.
    '''

    # Shape of the PSF cube
    nl,ny_PSF,nx_PSF = PSF_cube.shape

    # Number of rows and columns to add around a given pixel in order to extract a stamp.
    row_m = int(np.floor(ny_PSF/2.0))    # row_minus
    row_p = int(np.ceil(ny_PSF/2.0))     # row_plus
    col_m = int(np.floor(nx_PSF/2.0))    # col_minus
    col_p = int(np.ceil(nx_PSF/2.0))     # col_plus

    # Number of pixels on which the metric has to be computed
    N_it = row_indices.size
    # Define an shape vector full of nans
    shape_map = np.zeros((N_it,)) + np.nan
    # Loop over all pixels (row_indices[id],col_indices[id])
    for id,k,l in zip(range(N_it),row_indices,col_indices):
        if not mute:
            # Print the progress of the function
            stdout.write("\r{0}/{1}".format(id,N_it))
            stdout.flush()

        # Extract stamp cube around the current pixel from the whoel cube
        stamp_cube = copy(cube[:,(k-row_m):(k+row_p), (l-col_m):(l+col_p)])
        var_per_wv = np.zeros(nl)
        # Remove average value of the surrounding pixels in each slice of the stamp cube
        for slice_id in range(nl):
            stamp_cube[slice_id,:,:] -= np.nanmean(stamp_cube[slice_id,:,:]*stamp_PSF_mask)
            var_per_wv[slice_id] = np.nanvar(stamp_cube[slice_id,:,:]*stamp_PSF_mask)
        # Normalize the dot product square by the squared norm-2 of the stamp cube.
        # Because we keep the sign shape value is then found in [-1.,1.]
        try:
            # Dot product of the PSF with stamp cube.
            ampl = np.nansum(PSF_cube*stamp_cube/var_per_wv[:,None,None])
            # shape_map[id] = np.sign(ampl)*ampl**2/np.nansum(stamp_cube**2)
            shape_map[id] = ampl/np.sqrt(np.nansum(PSF_cube**2/var_per_wv[:,None,None]))
        except:
            # In case ones divide by zero...
            shape_map[id] =  np.nan

    # The shape value here can be seen as a cosine square as it is a normalized squared dot product.
    # Taking the square root to make it a simple cosine.
    # return np.sign(shape_map)*np.sqrt(abs(shape_map))
    return shape_map


def calculate_MF2D_metric_star(params):
    """
    Convert `f([1,2])` to `f(1,2)` call.
    It allows one to call calculate_MF3D_metric() with a tuple of parameters.
    """
    return calculate_MF2D_metric(*params)

def calculate_MF2D_metric(row_indices,col_indices,cube,PSF_cube,stamp_PSF_mask, mute = True):
    '''
    Calculate the matched filter metric on the given 2D image for the pixels targeted by row_indices and col_indices.
    These lists of indices can basically be given from the numpy.where function following the example:
        import numpy as np
        row_indices,col_indices = np.where(np.finite(np.mean(cube,axis=0)))
    By truncating the given lists in small pieces it is then easy to parallelized.

    The matched filter is a metric allowing one to pull out a known signal from noisy (not correlated) data.
    It is basically a dot product (meaning projection) of the template PSF cube with a stamp cube.

    :param row_indices: Row indices list of the pixels where to calculate the metric in cube.
                        Indices should be given from a 2d image.
    :param col_indices: Column indices list of the pixels where to calculate the metric in cube.
                        Indices should be given from a 2d image.
    :param cube: Cube from which one wants the metric map. PSF_cube should be norm-2 normalized.
                PSF_cube /= np.sqrt(np.sum(PSF_cube**2))
    :param PSF_flat: PSF template used for calculated the metric. If ny_PSF,nx_PSF = PSF_cube.shape ny_PSF and nx_PSF
                    are the spatial dimensions of the PSF_cube.
    :param stamp_PSF_mask: 2d mask of size (ny_PSF,nx_PSF) used to mask the central part of a stamp slice. It is used as
                        a type of a high pass filter. Before calculating the metric value of a stamp cube around a given
                        pixel the average value of the surroundings of each slice of that stamp cube will be removed.
                        The pixel used for calculating the average are the one equal to one in the mask.
    :param mute: If True prevent printed log outputs.
    :return: Vector of length row_indices.size with the value of the metric for the corresponding pixels.
    '''

    # Shape of the PSF cube
    ny_PSF,nx_PSF = PSF_cube.shape

    # Number of rows and columns to add around a given pixel in order to extract a stamp.
    row_m = np.floor(ny_PSF/2.0)    # row_minus
    row_p = np.ceil(ny_PSF/2.0)     # row_plus
    col_m = np.floor(nx_PSF/2.0)    # col_minus
    col_p = np.ceil(nx_PSF/2.0)     # col_plus

    # Number of pixels on which the metric has to be computed
    N_it = row_indices.size
    # Define an matched filter vector full of nans
    matchedFilter_map = np.zeros((N_it)) + np.nan
    # Loop over all pixels (row_indices[id],col_indices[id])
    for id,k,l in zip(range(N_it),row_indices,col_indices):
        if 1:#k == 108 and l == 135:
            if not mute:
                # Print the progress of the function
                stdout.write("\r{0}/{1}".format(id,N_it))
                stdout.flush()

            # Extract stamp cube around the current pixel from the whoel cube
            stamp_cube = copy(cube[(k-row_m):(k+row_p), (l-col_m):(l+col_p)])
            # Remove average value of the surrounding pixels in each slice of the stamp cube
            for slice_id in range(nl):
                stamp_cube[slice_id,:,:] -= np.nanmean(stamp_cube[slice_id,:,:]*stamp_PSF_mask)
                 #print(np.nanmean(stamp_cube[slice_id,:,:]*stamp_PSF_mask))
            # Dot product of the PSF with stamp cube which is the matched filter value.
            matchedFilter_map[id] = np.nansum(PSF_cube*stamp_cube)

    return matchedFilter_map

def calculate_MF3D_metric_star(params):
    """
    Convert `f([1,2])` to `f(1,2)` call.
    It allows one to call calculate_MF3D_metric() with a tuple of parameters.
    """
    return calculate_MF3D_metric(*params)

def calculate_MF3D_metric(row_indices,col_indices,cube,PSF_cube,stamp_PSF_mask, mute = True):
    '''
    Calculate the matched filter metric on the given datacube for the pixels targeted by row_indices and col_indices.
    These lists of indices can basically be given from the numpy.where function following the example:
        import numpy as np
        row_indices,col_indices = np.where(np.finite(np.mean(cube,axis=0)))
    By truncating the given lists in small pieces it is then easy to parallelized.

    The matched filter is a metric allowing one to pull out a known signal from noisy (not correlated) data.
    It is basically a dot product (meaning projection) of the template PSF cube with a stamp cube.

    :param row_indices: Row indices list of the pixels where to calculate the metric in cube.
                        Indices should be given from a 2d image.
    :param col_indices: Column indices list of the pixels where to calculate the metric in cube.
                        Indices should be given from a 2d image.
    :param cube: Cube from which one wants the metric map. PSF_cube should be norm-2 normalized.
                PSF_cube /= np.sqrt(np.sum(PSF_cube**2))
    :param PSF_cube: PSF_cube template used for calculated the metric. If nl,ny_PSF,nx_PSF = PSF_cube.shape, nl is the
                     number of wavelength samples, ny_PSF and nx_PSF are the spatial dimensions of the PSF_cube.
    :param stamp_PSF_mask: 2d mask of size (ny_PSF,nx_PSF) used to mask the central part of a stamp slice. It is used as
                        a type of a high pass filter. Before calculating the metric value of a stamp cube around a given
                        pixel the average value of the surroundings of each slice of that stamp cube will be removed.
                        The pixel used for calculating the average are the one equal to one in the mask.
    :param mute: If True prevent printed log outputs.
    :return: Vector of length row_indices.size with the value of the metric for the corresponding pixels.
    '''

    # Shape of the PSF cube
    nl,ny_PSF,nx_PSF = PSF_cube.shape

    # Number of rows and columns to add around a given pixel in order to extract a stamp.
    row_m = np.floor(ny_PSF/2.0)    # row_minus
    row_p = np.ceil(ny_PSF/2.0)     # row_plus
    col_m = np.floor(nx_PSF/2.0)    # col_minus
    col_p = np.ceil(nx_PSF/2.0)     # col_plus

    # Number of pixels on which the metric has to be computed
    N_it = row_indices.size
    # Define an matched filter vector full of nans
    matchedFilter_map = np.zeros((N_it)) + np.nan
    # Loop over all pixels (row_indices[id],col_indices[id])
    for id,k,l in zip(range(N_it),row_indices,col_indices):
        if 1:#k == 108 and l == 135:
            if not mute:
                # Print the progress of the function
                stdout.write("\r{0}/{1}".format(id,N_it))
                stdout.flush()

            # Extract stamp cube around the current pixel from the whoel cube
            stamp_cube = copy(cube[:,(k-row_m):(k+row_p), (l-col_m):(l+col_p)])
            # Remove average value of the surrounding pixels in each slice of the stamp cube
            for slice_id in range(nl):
                stamp_cube[slice_id,:,:] -= np.nanmean(stamp_cube[slice_id,:,:]*stamp_PSF_mask)
                 #print(np.nanmean(stamp_cube[slice_id,:,:]*stamp_PSF_mask))
            # Dot product of the PSF with stamp cube which is the matched filter value.
            matchedFilter_map[id] = np.nansum(PSF_cube*stamp_cube)

    return matchedFilter_map