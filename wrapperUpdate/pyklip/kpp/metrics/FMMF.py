__author__ = 'JB'

import os
from glob import glob

import astropy.io.fits as pyfits
import numpy as np
import scipy.interpolate as interp

from pyklip.kpp.utils.kppSuperClass import KPPSuperClass
from pyklip.instruments import GPI
import pyklip.spectra_management as spec
import pyklip.fm as fm
import pyklip.fmlib.matchedFilter as mf
import pyklip.klip as klip
import pyklip.kpp.utils.GOI as GOI
import pyklip.kpp.utils.GPIimage as gpiim


class FMMF(KPPSuperClass):
    """
    Forward model matched filter for GPI.

    /!\ Only work on campaign data.

    /!\ Caution: Most of the features need testing. Please contact jruffio@stanford.edu if there is a bug.

    Inherit from the Super class KPPSuperClass
    """
    def __init__(self,
                 filename = None,
                 inputDir = None,
                 outputDir = None,
                 PSF_cube_filename = None,
                 mute=None,
                 N_threads=None,
                 overwrite=False,
                 numbasis = None,
                 maxnumbasis = None,
                 flux_overlap = None,
                 mvt = None,
                 OWA = None,
                 N_pix_sector = None,
                 subsections = None,
                 annuli = None,
                 predefined_sectors = None,
                 label = None,
                 quickTest = False,
                 mute_progression = False,
                 SpT_file_csv = None,
                 fakes_only = None,
                 disable_FM = None,
                 high_pass_filter = None,
                 mvt_noTemplate = None,
                 true_fakes_pos = None):
        """
        Define the general parameters of the metric.

        The idea is that the directories and eventual spectra will be defined later by calling the
        initialize() function. Furthermore initialize() is where files can be read.
        It is done that way to ease the analysis of an entire campaign. ie one object is defined
        for a given metric and then a function calls initialize several time for each target's directory and each
        spectra.


        :param filename: Filename of the file on which to calculate the metric. It should be the complete path unless
                        inputDir is defined.
                        It can include wild characters. The file will be selected using the first output of glob().
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
                        A sequential option is hard coded for debugging purposes only.
        :param overwrite_metric: If True check_existence() will always return False.
        :param kernel_width: Define the width of the Kernel depending on kernel_type. See kernel_width.
        :param sky_aper_radius: Radius of the mask applied on the stamps to calculated the background value.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        :param quickTest: Read only two files (the first and the last) instead of the all sequence
        :param SpT_file_csv: Filename of the table (.csv) containing the spectral type of the stars.
        :param fakes_only:
        :param disable_FM:
        """
        # allocate super class
        super(FMMF, self).__init__(filename,
                                     inputDir = inputDir,
                                     outputDir = outputDir,
                                     folderName = None,
                                     mute=mute,
                                     N_threads=N_threads,
                                     label=label,
                                     overwrite=overwrite)

        # Prevent the class to iterate over all the files matching filename
        self.process_all_files = False
        self.quickTest = quickTest
        self.mute_progression = mute_progression
        self.SpT_file_csv = SpT_file_csv

        if filename is None:
            self.filename = "S*distorcorr.fits"
        else:
            self.filename = filename

        self.fakes_only = fakes_only
        if disable_FM is None:
            self.disable_FM = False
        else:
            self.disable_FM = disable_FM
        if high_pass_filter is None:
            self.high_pass_filter = True
        else:
            self.high_pass_filter = high_pass_filter
        if mvt_noTemplate is None:
            self.mvt_noTemplate = False
        else:
            self.mvt_noTemplate = mvt_noTemplate

        self.PSF_cube_filename = PSF_cube_filename

        self.save_per_sector = None
        self.padding = 10
        self.save_klipped = True
        self.true_fakes_pos = true_fakes_pos

        if numbasis is None:
            self.numbasis = np.array([30])
        else:
            self.numbasis = numbasis

        if maxnumbasis is None:
            self.maxnumbasis = 150
        else:
            self.maxnumbasis = maxnumbasis

        if flux_overlap is None:
            if mvt is None:
                self.flux_overlap = 0.7
                self.mvt = None
            else:
                self.flux_overlap = None
                self.mvt = mvt
        else:
            self.flux_overlap = flux_overlap
            self.mvt = None


        self.OWA = OWA
        self.N_pix_sector = N_pix_sector

        if subsections is None:
            self.subsections = 4
        else:
            self.subsections = subsections

        if annuli is None:
            self.annuli = 5
        else:
            self.annuli = annuli

        if predefined_sectors == "maskEdge":
            self.N_pix_sector = 100
            self.subsections = None
            self.annuli = [(8.698727015558699, 14.326080014734867), (14.326080014734867, 19.953433013911035), (19.953433013911035, 25.580786013087202)]
        if predefined_sectors == "smallSep":
            self.N_pix_sector = 100
            self.subsections = None
            # Define 3 thin annuli (each ~5pix wide) and a big one (~20pix) to cover up to 0.6''
            self.annuli = [(8.698727015558699, 14.326080014734867), (14.326080014734867, 19.953433013911035), (19.953433013911035, 25.580786013087202),(25.580786013087202, 42.46284501061571)]#
            # self.annuli = [(25.580786013087202, 42.46284501061571)]#
            # self.annuli = [(35, 50)]
        if predefined_sectors == "oneAc":
            self.N_pix_sector = 100
            self.subsections = None
            self.annuli = [(8.7, 14.3), (14.3, 20), (20, 25.6),(25.6, 40.5),(40.5,55.5),(55.5,70.8)]
            # self.annuli = [(66,81)]
            # self.annuli = [(55.5,70.8)]
        if predefined_sectors == "smallSepBigSec":
            self.OWA = gpiim.as2pix(0.6)
            self.N_pix_sector = 300
            self.subsections = None
            self.annuli = [(8.698727015558699, 19.953433013911035), (19.953433013911035, 42.46284501061571)]
        if predefined_sectors == "avgSep":
            self.OWA = gpiim.as2pix(0.8)
            self.N_pix_sector = 200
            self.subsections = None
            self.annuli = [(8.698727015558699, 14.326080014734867), (14.326080014734867, 19.953433013911035), (19.953433013911035, 25.580786013087202),(25.580786013087202, 41),(41, 56.5)]
        elif predefined_sectors == "c_Eri":
            self.subsections = [[150./180.*np.pi,190./180.*np.pi]]#None#[[150./180.*np.pi,190./180.*np.pi]]
            self.annuli = [[23,41]]
            # self.N_pix_sector = 100
        elif predefined_sectors == "fluxTest":
            self.subsections = [[(46-10)/180.*np.pi,(46+10)/180.*np.pi]]
            self.annuli = [[36-5,36+5]]
            self.N_pix_sector = None


    def spectrum_iter_available(self):
        """
        Indicates wether or not the class is equipped for iterating over spectra.
        Forward modelling matched filter requires a spectrum so the answer is yes.

        In order to iterate over spectra the function new_init_spectrum() can be called.
        spectrum_iter_available is a utility function for campaign data processing.

        :return: True
        """

        return True

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

        if spectrum is None:
            self.spectrum_name = "t600g100nc"
            # spec_path = "g32ncflx"+os.path.sep+self.spectrum_name+".flx"
        else:
            self.spectrum_name = spectrum
            # self.spectrum_name = spectrum.split(os.path.sep)[-1].split(".")[0]
            # spec_path = spectrum

        self.folderName = self.spectrum_name+os.path.sep

        if self.flux_overlap is not None:
            self.prefix = self.star_name+"_"+self.compact_date+"_"+self.filter+"_"+self.spectrum_name +"_{0:.2f}".format(self.flux_overlap)
        else:
            self.prefix = self.star_name+"_"+self.compact_date+"_"+self.filter+"_"+self.spectrum_name +"_{0:.2f}".format(self.mvt)

        # methane spectral template
        pykliproot = os.path.dirname(os.path.realpath(klip.__file__))
        self.spectrum_filename = os.path.abspath(glob(os.path.join(pykliproot,"spectra","*",self.spectrum_name+".flx"))[0])
        if self.mvt_noTemplate:
            self.spectra_template = None
        else:
            # spectrum_dat = np.loadtxt(self.spectrum_filename )[:160] #skip wavelegnths longer of 10 microns
            # spectrum_wvs = spectrum_dat[:,1]
            # spectrum_fluxes = spectrum_dat[:,3]
            # spectrum_interpolation = interp.interp1d(spectrum_wvs, spectrum_fluxes, kind='cubic')
            # # This spectrum is the one used by klip itself while the matched filter uses the one defined in fm_class.
            # # This should however be the same. I know it's weird but the two codes are seperate
            # self.spectra_template = spectrum_interpolation(self.dataset.wvs)
            self.spectra_template = spec.get_planet_spectrum(self.spectrum_filename,self.dataset.wvs)[1]


        # Build the FM class to do matched filter
        self.fm_class = mf.MatchedFilter(self.dataset.input.shape,self.numbasis, self.dataset.psfs, np.unique(self.dataset.wvs),self.dataset.spot_flux,
                                     spectrallib = [spec.get_planet_spectrum(filename, self.filter)[1] for filename in [self.spectrum_filename ]],
                                     mute = False,
                                     star_type = self.star_type,
                                     filter_name = self.filter,
                                     save_per_sector = self.save_per_sector,
                                     fakes_sepPa_list = self.fakes_sepPa_list,
                                     disable_FM=self.disable_FM,
                                         true_fakes_pos = self.true_fakes_pos)
        return None

    def initialize(self, inputDir = None,
                         outputDir = None,
                         spectrum = None,
                         folderName = None,
                         PSF_cube_filename = None,
                         star_type = "auto",
                         compact_date = None,
                         label=None):
        """

        Initialize the non general inputs that are needed for the metric calculation and load required files.


        - It reads the input fits file to be reduced.
        - Build or load the PSF if needed
        - load the spectrum if needed (it can be redefined later using self.init_new_spectrum())
        - read headers if they exist

        The injected model have a spectrum normalized to unit broadband flux in DN space.


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
                        Spectra are assumed to be in FLUX UNITS.
        :param folderName: Name of the folder containing the outputs. It will be located in outputDir.
                        Default folder name is "default_out" if not spectrum is defined or the name of the spectrum
                        otherwise.
        :param PSF_cube_filename: Filename filter for the PSF in inputdir.
                        If inputDir is not defined it should be an absolute path.
                        Default value is "*-original_radial_PSF_cube.fits".
                        Useful only if kernel_type = "PSF"
                        If a PSF cube is not explicitly given and one is read automatically it assumes there is only
                        one PSF cube in this folder.
        :param star_type: String containing the spectral type of the star. 'A5','F4',... Assume type V star.
                        If "auto" (default), spec.get_specType(self.star_name,self.SpT_file_csv) is called.
                        If self.SpT_file_csv is defined then it reads the spectral type from the table otherwise try to
                        query simbad directly.
                        If None, the spectra are assumed to be in units of contrast.
        :param compact_date: Define the compact date to be used in the output filenames.
                            If a PSF has to be measured from the satellite spots it will define itself.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        :return:
        """
        if not self.mute:
            print("~~ INITializing "+self.__class__.__name__+" ~~")

        # The super class already read the fits file
        init_out = super(FMMF, self).initialize(inputDir = inputDir,
                                         outputDir = outputDir,
                                         folderName = folderName,
                                         label=label)


        if compact_date is None:
            self.compact_date = "noDate"
        else:
            self.compact_date = compact_date

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

        if star_type is "auto":
            self.star_type = spec.get_specType(self.star_name,self.SpT_file_csv)
        else:
            self.star_type = star_type
        if not self.mute:
            print("The star spectral type is: {0}".format(self.star_type))

        # Get the list of spdc files
        filelist = glob(self.inputDir+os.path.sep+self.filename)
        filelist.sort()
        if self.quickTest:
            filelist = [filelist[0],filelist[-1]]

        # Filename of the PSF cube
        if PSF_cube_filename is not None:
            self.PSF_cube_filename = PSF_cube_filename

        if self.PSF_cube_filename == None:
            self.PSF_cube_filename = "*-original_PSF_cube.fits"
            if not self.mute:
                print("Using default filename for PSF cube: "+self.PSF_cube_filename)

        # Define the actual filename path
        if self.inputDir is None:
            try:
                self.PSF_cube_path = os.path.abspath(glob(self.PSF_cube_filename)[0])
                if not self.mute:
                    print("Loading PSF cube: "+self.PSF_cube_path)
                hdulist = pyfits.open(self.PSF_cube_path)
                PSF_cube = hdulist[1].data
                self.dataset = GPI.GPIData(filelist,highpass=self.high_pass_filter,meas_satspot_flux=True,numthreads=None,PSF_cube=PSF_cube)
                self.dataset.psfs = PSF_cube
            except:
                raise Exception("File "+self.PSF_cube_filename+"doesn't exist.")
        else:
            try:
                self.PSF_cube_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.PSF_cube_filename)[0])
                if not self.mute:
                    print("Loading PSF cube: "+self.PSF_cube_path)
                hdulist = pyfits.open(self.PSF_cube_path)
                PSF_cube = hdulist[1].data
                self.dataset = GPI.GPIData(filelist,highpass=self.high_pass_filter,meas_satspot_flux=True,numthreads=None,PSF_cube=PSF_cube)
                self.dataset.psfs = PSF_cube
            except:
                prefix = self.star_name+"_"+self.compact_date+"_"+self.filter

                if not self.mute:
                    print("Calculating the planet PSF from the satellite spots...")
                self.dataset = GPI.GPIData(filelist,highpass=self.high_pass_filter,meas_satspot_flux=False,numthreads=None)
                # generate the PSF cube from the satellite spots
                self.dataset.generate_psf_cube(20,same_wv_only=True)
                # Save the original PSF calculated from combining the sat spots
                self.dataset.savedata(self.inputDir + os.path.sep + prefix+"-original_PSF_cube.fits", self.dataset.psfs,
                                          astr_hdr=self.dataset.wcs[0], filetype="PSF Spec Cube")
                # Calculate and save the rotationally invariant psf (ie smeared out/averaged).
                # Generate a radially averaged PSF cube and save as
                # self.inputDir + os.path.sep + prefix+"-original_radial_PSF_cube.fits"
                self.dataset.get_radial_psf(save = self.inputDir + os.path.sep + prefix)
                self.PSF_cube_path = self.inputDir + os.path.sep + prefix+"-original_radial_PSF_cube.fits"
                self.dataset.spot_flux = GPI.recalculate_sat_spot_fluxes(self.dataset,numthreads=None,residuals=False)
                if not self.mute:
                    print(self.inputDir + os.path.sep + prefix+"-original_radial_PSF_cube.fits"+"and"+\
                          self.inputDir + os.path.sep + prefix+"-original_PSF_cube.fits"+ " have been saved.")


        if spectrum is None:
            self.spectrum_name = "t600g100nc"
            # spec_path = "g32ncflx"+os.path.sep+self.spectrum_name+".flx"
        else:
            self.spectrum_name = spectrum
            # self.spectrum_name = spectrum.split(os.path.sep)[-1].split(".")[0]
            # spec_path = spectrum

        self.folderName = self.spectrum_name+os.path.sep

        if self.flux_overlap is not None:
            self.prefix = self.star_name+"_"+self.compact_date+"_"+self.filter+"_"+self.spectrum_name +"_{0:.2f}".format(self.flux_overlap)
        else:
            self.prefix = self.star_name+"_"+self.compact_date+"_"+self.filter+"_"+self.spectrum_name +"_{0:.2f}".format(self.mvt)

        #todo if self.fakes_only
        # read fakes from headers and give sepPa list to MatchedFilter
        # else give None
        if self.fakes_only:
            sep_list, pa_list = GOI.get_pos_known_objects(self.prihdr,self.exthdr,GOI_list_folder=None,xy = False,pa_sep = True,fakes_only=True)
            self.fakes_sepPa_list = [(gpiim.as2pix(sep),pa) for sep,pa in zip(sep_list, pa_list)]
        else:
            self.fakes_sepPa_list = None

        # methane spectral template
        pykliproot = os.path.dirname(os.path.realpath(klip.__file__))
        self.spectrum_filename = os.path.abspath(glob(os.path.join(pykliproot,"spectra","*",self.spectrum_name+".flx"))[0])

        if self.mvt_noTemplate:
            self.spectra_template = None
        else:
            # spectrum_dat = np.loadtxt(self.spectrum_filename )[:160] #skip wavelegnths longer of 10 microns
            # spectrum_wvs = spectrum_dat[:,1]
            # spectrum_fluxes = spectrum_dat[:,3]
            # spectrum_interpolation = interp.interp1d(spectrum_wvs, spectrum_fluxes, kind='cubic')
            # # This spectrum is the one used by klip itself while the matched filter uses the one defined in fm_class.
            # # This should however be the same. I know it's weird but the two codes are seperate
            # self.spectra_template = spectrum_interpolation(self.dataset.wvs)
            self.spectra_template = spec.get_planet_spectrum(self.spectrum_filename,self.dataset.wvs)[1]


        # Build the FM class to do matched filter
        self.fm_class = mf.MatchedFilter(self.dataset.input.shape,self.numbasis, self.dataset.psfs, np.unique(self.dataset.wvs),self.dataset.spot_flux,
                                     spectrallib = [spec.get_planet_spectrum(self.spectrum_filename, self.filter)[1]],
                                     mute = False,
                                     star_type = self.star_type,
                                     filter_name = self.filter,
                                     save_per_sector = self.save_per_sector,
                                     fakes_sepPa_list = self.fakes_sepPa_list,
                                     disable_FM=self.disable_FM,
                                         true_fakes_pos= self.true_fakes_pos)



        return init_out

    def check_existence(self):
        """
        Check if this metric has already been calculated for this file.

        If self.overwrite is True then the output will always be False

        :return: Boolean indicating the existence of the metric map.
        """

        if self.quickTest:
            susuffix = "QT"
        else:
            susuffix = ""

        if self.disable_FM:
            presuffix = "no"
        else:
            presuffix = ""

        suffix1 = presuffix+"FMMF"+susuffix
        file_exist=(len(glob(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix1+'.fits')) >= 1)

        if file_exist and not self.mute:
            print("Output already exist.")

        return file_exist and not self.overwrite


    def calculate(self):
        """
        Calculate the metric map.

        :return: [self.metric_MF,self.metric_shape,self.sub_imgs].
        """


        # Run KLIP with the forward model matched filter
        sub_imgs, fmout,tmp = fm.klip_parallelized(self.dataset.input, self.dataset.centers, self.dataset.PAs, self.dataset.wvs, self.dataset.IWA, self.fm_class,
                                   numbasis=self.numbasis,
                                   maxnumbasis=self.maxnumbasis,
                                   flux_overlap=self.flux_overlap,
                                   movement=self.mvt,
                                   spectrum=self.spectra_template,
                                   annuli=self.annuli,
                                   subsections=self.subsections,
                                   numthreads=self.N_threads,
                                   padding = self.padding,
                                   N_pix_sector=self.N_pix_sector,
                                   save_klipped = self.save_klipped,
                                   OWA = self.OWA,
                                   mute_progression = self.mute_progression)

        #fmout_shape = (3,self.N_spectra,self.N_numbasis,self.N_frames,self.ny,self.nx)
        self.N_cubes = fmout.shape[3]/37
        # methane_indices = []
        # other_indices = []
        # for k in range(self.N_cubes):
        #     methane_indices.extend(range(k*self.nl,(k+1)*self.nl-19))
        #     other_indices.extend(range((k+1)*self.nl-19,(k+1)*self.nl))
        # self.matched_filter_maps_methane = np.squeeze(np.nansum(fmout[0,:,:,methane_indices,:,:],axis=0))
        # self.matched_filter_maps_other = np.squeeze(np.nansum(fmout[0,:,:,other_indices,:,:],axis=0))
        # self.matched_filter_maps_methane[np.where(self.matched_filter_maps_methane==0)]=np.nan
        # self.matched_filter_maps_other[np.where(self.matched_filter_maps_other ==0)]=np.nan

        fmout[np.where(fmout==0)] = np.nan

        # Retrieve the flux map, mathced filter map and shape map
        # Note: the injected model had a unit broadband flux in DN space.

        #######################
        ## Assuming the noise is constant for all slices of a cube
        metric_shape_nume = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,self.ny,self.nx))
        metric_shape_deno = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,self.ny,self.nx))

        for k in range(self.N_cubes):
            var_per_cube = np.sum(fmout[2,:,:,k*self.nl:(k+1)*self.nl,:,:],axis=2)/self.nl
            metric_shape_nume = metric_shape_nume+np.sum(fmout[0,:,:,k*self.nl:(k+1)*self.nl,:,:]/var_per_cube[:,:,None,:,:],axis=2)
            metric_shape_deno = metric_shape_deno+np.sum(fmout[1,:,:,k*self.nl:(k+1)*self.nl,:,:]/var_per_cube[:,:,None,:,:],axis=2)


        metric_shape = metric_shape_nume/np.sqrt(metric_shape_deno)
        metric_shape[np.where(metric_shape==0)]=np.nan
        metric_shape = np.squeeze(metric_shape)


        # self.metric_pFlux3 = metric_pFlux
        self.metric_shape3 = metric_shape

        #######################
        ## matched filter for each cube and then averaged
        # ppm_band = self.prihdr['APODIZER'].split('_')[1] #to determine sat spot ratios
        # sat_spot_ratio = self.dataset.spot_ratio[ppm_band]
        # metric_pFlux_nume = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,self.ny,self.nx))
        # metric_pFlux_deno = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,self.ny,self.nx))


        metric_CC = np.sum(fmout[0,:,:,:,:,:],axis=2) \
                        / np.sqrt(np.sum(fmout[1,:,:,:,:,:],axis=2))
        metric_CC[np.where(metric_shape==0)]=np.nan
        metric_CC = np.squeeze(metric_CC)
        self.metric_CC2 = metric_CC

        metric_shape = np.sum(fmout[0,:,:,:,:,:]/fmout[2,:,:,:,:,:],axis=2) \
                        / np.sqrt(np.sum(fmout[1,:,:,:,:,:]/fmout[2,:,:,:,:,:],axis=2))
        metric_shape[np.where(metric_shape==0)]=np.nan
        metric_shape = np.squeeze(metric_shape)
        self.metric_shape2 = metric_shape

        metric_pFlux = np.sum(fmout[0,:,:,:,:,:]/fmout[2,:,:,:,:,:],axis=2) \
                        / np.sum(fmout[1,:,:,:,:,:]/fmout[2,:,:,:,:,:],axis=2)
        metric_pFlux[np.where(metric_pFlux==0)]=np.nan
        metric_pFlux = np.squeeze(metric_pFlux)
        self.metric_pFlux2 = metric_pFlux

        # for k in range(self.N_cubes):
        #     sat_spot_flux_for_calib = np.nansum(self.dataset.spot_flux[k*self.nl:(k+1)*self.nl]*self.fm_class.aper_over_peak_ratio)
        #     metric_pFlux_nume = metric_pFlux_nume + np.nansum(fmout[0,:,:,k*self.nl:(k+1)*self.nl,:,:]/fmout[2,:,:,k*self.nl:(k+1)*self.nl,:,:],axis=2) / sat_spot_flux_for_calib * sat_spot_ratio
        #     metric_pFlux_deno = metric_pFlux_deno + np.nansum(fmout[1,:,:,k*self.nl:(k+1)*self.nl,:,:]/fmout[2,:,:,k*self.nl:(k+1)*self.nl,:,:],axis=2)
        # metric_pFlux = metric_pFlux_nume/metric_pFlux_deno
        # # metric_pFlux = metric_pFlux/self.N_cubes
        # metric_pFlux[np.where(metric_pFlux==0)]=np.nan
        # metric_pFlux = np.squeeze(metric_pFlux)
        #
        # # self.metric_MF2 = metric_MF
        # self.metric_pFlux2 = metric_pFlux

        #######################
        ## matched filter for each cube and then averaged
        metric_MF = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,self.ny,self.nx))
        metric_shape = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,self.ny,self.nx))
        for k in range(self.N_cubes):
            metric_MF = metric_MF+np.sum(fmout[0,:,:,k*self.nl:(k+1)*self.nl,:,:],axis=2) \
                            / np.sqrt(np.sum(fmout[1,:,:,k*self.nl:(k+1)*self.nl,:,:],axis=2))
            metric_shape = metric_shape+np.sum(fmout[0,:,:,k*self.nl:(k+1)*self.nl,:,:],axis=2) \
                            / (  np.sqrt(np.sum(fmout[1,:,:,k*self.nl:(k+1)*self.nl,:,:],axis=2)) \
                               * np.sqrt(np.sum(fmout[2,:,:,k*self.nl:(k+1)*self.nl,:,:],axis=2)) )
        metric_MF = metric_MF/self.N_cubes
        metric_MF[np.where(metric_MF==0)]=np.nan
        metric_MF = np.squeeze(metric_MF)
        metric_shape = metric_shape/self.N_cubes
        metric_shape[np.where(metric_shape==0)]=np.nan
        metric_shape = np.squeeze(metric_shape)

        self.metric_MF = metric_MF
        self.metric_shape = metric_shape

        #######################
        ## cubes are coadded and then matched filter is applied
        # fmout0 = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,37,self.ny,self.nx))
        # fmout1 = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,37,self.ny,self.nx))
        # fmout2 = np.zeros((self.fm_class.N_spectra,self.fm_class.N_numbasis,37,self.ny,self.nx))
        # for k in range(self.N_cubes): #JUST RETHINK IT
        #     fmout0 = fmout0+fmout[0,:,:,k*self.nl:(k+1)*self.nl,:,:]
        #     fmout1 = fmout1+fmout[1,:,:,k*self.nl:(k+1)*self.nl,:,:]
        #     fmout2 = fmout0+fmout[2,:,:,k*self.nl:(k+1)*self.nl,:,:]
        # fmout0 = fmout0/37
        # fmout1 = fmout1/37
        # fmout2 = fmout2/37
        # sat_spot_flux_for_calib = 0
        # for k in range(self.N_cubes):
        #     sat_spot_flux_for_calib = sat_spot_flux_for_calib+np.nansum(self.dataset.spot_flux[k*self.nl:(k+1)*self.nl]*self.fm_class.aper_over_peak_ratio)
        # sat_spot_flux_for_calib = sat_spot_flux_for_calib/self.N_cubes
        # metric_pFlux = np.nansum(fmout0,axis=2) \
        #                 / np.nansum(fmout1,axis=2) \
        #                 / sat_spot_flux_for_calib * sat_spot_ratio
        # metric_MF = np.sum(fmout0,axis=2) / np.sqrt(np.sum(fmout1,axis=2))
        # metric_shape = np.sum(fmout0,axis=2) \
        #                 / (  np.sqrt(np.sum(fmout1,axis=2)) \
        #                    * np.sqrt(np.sum(fmout2,axis=2)) )
        #
        # metric_pFlux[np.where(metric_pFlux==0)]=np.nan
        # metric_pFlux = np.squeeze(metric_pFlux)
        # metric_MF[np.where(metric_MF==0)]=np.nan
        # metric_MF = np.squeeze(metric_MF)
        # metric_shape[np.where(metric_shape==0)]=np.nan
        # metric_shape = np.squeeze(metric_shape)

        #self.dataset.psfs



        # # Build the matched filter and shape maps from fmout
        # matched_filter_maps = np.sum(fmout[0,:,:,:,:,:],axis=2)
        # model_square_norm_maps = np.sum(fmout[1,:,:,:,:,:],axis=2)
        # image_square_norm_maps = np.sum(fmout[2,:,:,:,:,:],axis=2)
        # matched_filter_maps[np.where(matched_filter_maps==0)]=np.nan
        # model_square_norm_maps[np.where(model_square_norm_maps==0)]=np.nan
        # image_square_norm_maps[np.where(image_square_norm_maps==0)]=np.nan
        # metric_MF = matched_filter_maps/np.sqrt(model_square_norm_maps)
        # metric_shape = matched_filter_maps/(np.sqrt(model_square_norm_maps)*np.sqrt(image_square_norm_maps))
        # metric_MF = np.squeeze(metric_MF)
        # metric_shape = np.squeeze(metric_shape)

        # Update the wcs headers to indicate North up
        [klip._rotate_wcs_hdr(astr_hdr, angle, flipx=True) for angle, astr_hdr in zip(self.dataset.PAs, self.dataset.wcs)]

        self.sub_imgs = sub_imgs
        self.metricMap = [self.metric_MF,self.metric_shape,self.sub_imgs]


        self.N_cubes = sub_imgs.shape[1]/self.nl
        #print(sub_imgs.shape)
        cubes_list = []
        for k in range(self.N_cubes):
            cubes_list.append(sub_imgs[:,k*self.nl:(k+1)*self.nl,:,:])
        self.final_cube_modes = np.sum(cubes_list,axis = 0)
        self.final_cube_modes[np.where(self.final_cube_modes==0)] = np.nan
        #print(self.final_cube_modes.shape)

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

        # Save the parameters as fits keywords
        extra_exthdr_keywords = [("METNUMBA",str(self.numbasis)),
                                 ("METMAXNB",self.maxnumbasis),
                                 ("METFLXOV",str(self.flux_overlap)),
                                 ("MET_MVT",str(self.mvt)),
                                 ("MET_OWA",str(self.OWA)),
                                 ("METNPIXS",self.N_pix_sector),
                                 ("METSUBSE",str(self.subsections)),
                                 ("METANNUL",str(self.annuli)),
                                 ("METMAXNB",self.maxnumbasis),
                                 ("METFILEN",self.filename_path),
                                 ("METINDIR",self.inputDir),
                                 ("METOUTDI",self.outputDir),
                                 ("METFOLDN",self.folderName),
                                 ("METCDATE",self.compact_date),
                                 ("METSPECN",self.spectrum_name),
                                 ("METSPECF",self.spectrum_filename),
                                 ("METPSFDI",self.PSF_cube_path),
                                 ("METPREFI",self.prefix),
                                 ("METQUICT",self.quickTest)]

        if hasattr(self,"star_type"):
            extra_exthdr_keywords.append(("METSTTYP",self.star_type))

        extra_exthdr_keywords.append(("FMMFVERS",2.0))

        if self.quickTest:
            susuffix = "QT"
        else:
            susuffix = ""

        if self.disable_FM:
            presuffix = "no"
        else:
            presuffix = ""

        # Save the outputs (matched filter, shape map and klipped image) as fits files
        suffix = presuffix+"FMMF"+susuffix
        extra_exthdr_keywords.append(("METSUFFI",suffix))
        self.dataset.savedata(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits',
                         self.metric_shape2,
                         filetype=suffix,
                         astr_hdr=self.dataset.wcs[0], center=self.dataset.centers[0],
                         extra_exthdr_keywords = extra_exthdr_keywords)
        suffix = presuffix+"FMCont"+susuffix
        extra_exthdr_keywords[-1] = ("METSUFFI",suffix)
        self.dataset.savedata(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits',
                         self.metric_pFlux2,
                         filetype=suffix,
                         astr_hdr=self.dataset.wcs[0], center=self.dataset.centers[0],
                         extra_exthdr_keywords = extra_exthdr_keywords)
        suffix = presuffix+"FMCC"+susuffix
        extra_exthdr_keywords[-1] = ("METSUFFI",suffix)
        self.dataset.savedata(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits',
                         self.metric_CC2,
                         filetype=suffix,
                         astr_hdr=self.dataset.wcs[0], center=self.dataset.centers[0],
                         extra_exthdr_keywords = extra_exthdr_keywords)
        for k in range(self.final_cube_modes.shape[0]):
            suffix = "speccube-KL{0}".format(self.numbasis[k])+susuffix
            extra_exthdr_keywords[-1] = ("METSUFFI",suffix)
            self.dataset.savedata(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits',
                             self.final_cube_modes[k],
                             filetype=suffix,
                             astr_hdr=self.dataset.wcs[0], center=self.dataset.centers[0],
                             extra_exthdr_keywords = extra_exthdr_keywords)

        suffix = presuffix+"FMMFv1"+susuffix
        extra_exthdr_keywords[-1] = ("METSUFFI",suffix)
        self.dataset.savedata(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits',
                         self.metric_shape,
                         filetype=suffix,
                         astr_hdr=self.dataset.wcs[0], center=self.dataset.centers[0],
                         extra_exthdr_keywords = extra_exthdr_keywords)
        suffix = presuffix+"FMCCv1"+susuffix
        extra_exthdr_keywords[-1] = ("METSUFFI",suffix)
        self.dataset.savedata(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits',
                         self.metric_MF,
                         filetype=suffix,
                         astr_hdr=self.dataset.wcs[0], center=self.dataset.centers[0],
                         extra_exthdr_keywords = extra_exthdr_keywords)

        suffix = presuffix+"FMMF2"+susuffix
        extra_exthdr_keywords[-1] = ("METSUFFI",suffix)
        self.dataset.savedata(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits',
                         self.metric_shape3,
                         filetype=suffix,
                         astr_hdr=self.dataset.wcs[0], center=self.dataset.centers[0],
                         extra_exthdr_keywords = extra_exthdr_keywords)
        return None


    def load(self):
        """
        Load the metric map. One should check that it exist first using self.check_existence().

        Define the attribute self.metricMap.

        :return: self.metricMap
        """
        # suffix = "FMMF"
        # hdulist = pyfits.open(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits')
        # self.metric_MF = hdulist[1].data
        # hdulist.close()
        # suffix = "FMSH"
        # hdulist = pyfits.open(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits')
        # self.metric_shape = hdulist[1].data
        # hdulist.close()
        # suffix = "speccube-PSFs"
        # hdulist = pyfits.open(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+suffix+'.fits')
        # self.sub_imgs = hdulist[1].data
        # hdulist.close()
        # self.metricMap = [self.metric_MF,self.metric_shape,self.sub_imgs]


        return None
