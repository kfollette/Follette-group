__author__ = 'jruffio'
import os
import astropy.io.fits as pyfits
from glob import glob
import multiprocessing as mp
import numpy as np
from scipy.signal import convolve2d

from pyklip.kpp.utils.kppSuperClass import KPPSuperClass
from pyklip.kpp.stat.statPerPix_utils import *
from pyklip.kpp.stat.statPerPix import StatPerPix
from pyklip.kpp.utils.GOI import *
from pyklip.kpp.utils.GPIimage import *
import pyklip.instruments.GPI as GPI
from pyklip.kpp.metrics.FMMF import FMMF
from pyklip.kpp.metrics.crossCorr import CrossCorr
import pyklip.parallelized as parallelized
from pyklip.kpp.stat.stat import Stat
import pyklip.kpp.utils.mathfunc as kppmath
from pyklip.kpp.metrics.shapeOrMF import ShapeOrMF
from pyklip.kpp.kppPerDir import *

class Contrast(KPPSuperClass):
    """
    Class for SNR calculation.
    """
    def __init__(self,filename,dir_fakes,PSF_cube_filename = None,
                 inputDir = None,
                 outputDir = None,
                 mute=None,
                 N_threads=None,
                 label = None,
                 mask_radius = None,
                 IOWA = None,
                 GOI_list_folder = None,
                 GPI_TSpT_csv = None,
                 overwrite = False,
                 contrast_filename = None,
                 fakes_SNR = None,
                 fakes_spectrum = None,
                 spectrum_name=None,
                 reduce_only = False,
                 plot_only = False,
                 save_contrast=None):
        """

        :param filename: Filename of the file on which to calculate the metric. It should be the complete path unless
                        inputDir is defined.
                        It can include wild characters. The file will be selected using the first output of glob.glob().
        :param mute: If True prevent printed log outputs.
        :param N_threads: Number of threads to be used for the metrics and the probability calculations.
                        If None use mp.cpu_count().
                        If -1 do it sequentially.
                        Note that it is not used for this super class.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        """
        # allocate super class
        super(Contrast, self).__init__(filename,
                                     inputDir = inputDir,
                                     outputDir = outputDir,
                                     folderName = None,
                                     mute=mute,
                                     N_threads=N_threads,
                                     label=label,
                                     overwrite = overwrite)

        if mask_radius is None:
            self.mask_radius = 3
        else:
            self.mask_radius = mask_radius

        if fakes_SNR is None:
            self.fakes_SNR = 10
        else:
            self.fakes_SNR = fakes_SNR

        if fakes_spectrum is None:
            self.fakes_spectrum = "t900g100nc"
        else:
            self.fakes_spectrum = fakes_spectrum

        if spectrum_name is None:
            self.spectrum_name = self.spectrum_name
        else:
            self.spectrum_name = spectrum_name

        self.save_contrast=save_contrast

        self.IOWA = IOWA
        self.Dr = 2
        self.type = "stddev"
        self.suffix = "2Dcontrast"
        self.GOI_list_folder = GOI_list_folder
        self.GPI_TSpT_csv = GPI_TSpT_csv
        self.contrast_filename = contrast_filename
        self.PSF_cube_filename = PSF_cube_filename
        self.dir_fakes = dir_fakes
        self.reduce_only = reduce_only
        self.plot_only = plot_only


    def initialize(self,inputDir = None,
                         outputDir = None,
                         folderName = None,
                         compact_date = None,
                         label = None):
        """
        Initialize the non general inputs that are needed for the metric calculation and load required files.

        For this super class it simply reads the input file including fits headers and store it in self.image.
        One can also overwrite inputDir, outputDir which is basically the point of this function.
        The file is assumed here to be a fits containing a 2D image or a GPI 3D cube (assumes 37 spectral slice).

        Example for inherited classes:
        It can read the PSF cube or define the hat function.
        It can also read the template spectrum in a 3D scenario.
        It could also overwrite this function in case it needs to read multiple files or non fits file.

        :param inputDir: If defined it allows filename to not include the whole path and just the filename.
                        Files will be read from inputDir.
                        Note tat inputDir might be redefined using initialize at any point.
                        If inputDir is None then filename is assumed to have the absolute path.
        :param outputDir: Directory where to create the folder containing the outputs.
                        Note tat inputDir might be redefined using initialize at any point.
                        If outputDir is None:
                            If inputDir is defined: outputDir = inputDir+os.path.sep+"planet_detec_"
        :param folderName: Name of the folder containing the outputs. It will be located in outputDir.
                        Default folder name is "default_out".
                        The convention is to have one folder per spectral template.
                        If the keyword METFOLDN is available in the fits file header then the keyword value is used no
                        matter the input.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".

        :return: None
        """
        if not self.mute:
            print("~~ INITializing "+self.__class__.__name__+" ~~")
        # The super class already read the fits file
        init_out = super(Contrast, self).initialize(inputDir = inputDir,
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

        if self.contrast_filename is not None:
            # Check file existence and define filename_path
            if self.inputDir is None or os.path.isabs(self.contrast_filename):
                try:
                    if len(glob(self.contrast_filename)) == self.N_matching_files:
                        self.contrast_filename_path = os.path.abspath(glob(self.contrast_filename)[self.id_matching_file-1])
                    else:
                        self.contrast_filename_path = os.path.abspath(glob(self.contrast_filename)[0])
                except:
                    raise Exception("File "+self.contrast_filename+"doesn't exist.")
            else:
                try:
                    if len(glob(self.inputDir+os.path.sep+self.contrast_filename)) == self.N_matching_files:
                        self.contrast_filename_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.contrast_filename)[self.id_matching_file-1])
                    else:
                        self.contrast_filename_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.contrast_filename)[0])
                except:
                    raise Exception("File "+self.inputDir+os.path.sep+self.contrast_filename+" doesn't exist.")

        if self.PSF_cube_filename is not None:
            # Check file existence and define filename_path
            if self.inputDir is None or os.path.isabs(self.PSF_cube_filename):
                try:
                    self.PSF_cube_path = os.path.abspath(glob(self.PSF_cube_filename)[self.id_matching_file])
                except:
                    raise Exception("File "+self.PSF_cube_filename+"doesn't exist.")
            else:
                try:
                    self.PSF_cube_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.PSF_cube_filename)[self.id_matching_file])
                except:
                    raise Exception("File "+self.inputDir+os.path.sep+self.PSF_cube_filename+" doesn't exist.")

            # Open the fits file on which the metric will be applied
            hdulist = pyfits.open(self.PSF_cube_path)
            if not self.mute:
                print("Opened: "+self.PSF_cube_path)

            # grab the data and headers
            try:
                self.PSF_cube = hdulist[1].data
            except:
                raise Exception("Couldn't read "+self.PSF_cube_path+". Is it a .fits file with GPI extensions convention?")

        # Get center of the image (star position)
        try:
            # Retrieve the center of the image from the fits headers.
            self.center = [self.exthdr['PSFCENTX'], self.exthdr['PSFCENTY']]
        except:
            # If the keywords could not be found the center is defined as the middle of the image
            if not self.mute:
                print("Couldn't find PSFCENTX and PSFCENTY keywords.")
            self.center = [(self.nx-1)/2,(self.ny-1)/2]


        try:
            self.folderName = self.exthdr["METFOLDN"]+os.path.sep
        except:
            pass

        file_ext_ind = os.path.basename(self.filename_path)[::-1].find(".")
        self.prefix = os.path.basename(self.filename_path)[:-(file_ext_ind+1)]

        return init_out

    def check_existence(self):
        """

        :return: False
        """

        file_exist = (len(glob(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')) >= 1)

        if file_exist and not self.mute:
            print("Output already exist: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')

        if self.overwrite and not self.mute:
            print("Overwriting is turned ON!")

        return file_exist and not self.overwrite


    def calculate(self):
        """

        :param N: Defines the width of the ring by the number of pixels it has to contain
        :return: self.image the imput fits file.
        """
        if not self.mute:
            print("~~ Calculating "+self.__class__.__name__+" with parameters " + self.suffix+" ~~")

        # Inject fakes
        print(self.inputDir)

        throughput_break = 1.1
        contrast_range = [0.2,1.5]

        inputDir_tasks = []
        fakesDir_tasks = []

        overwrite_tmp = False
        resolution = 3.5
        pa_shift_list = [0,30,60]

        if not os.path.exists(self.dir_fakes):
            os.makedirs(self.dir_fakes)

        if not self.plot_only:
            if 0:
                if overwrite_tmp or len(glob(os.path.join(self.inputDir,"pyklip_k150a9s4m1methane_PSFsatSpotFlux","pyklip_k150a9s4m1methane-KL50-speccube.fits"))) == 0:
                    if not self.mute:
                        print("~~ Reducing pyklip no fakes ~~")
                    spdc_glob = glob(self.inputDir+os.path.sep+"S*_spdc_distorcorr.fits")
                    dataset = GPI.GPIData(spdc_glob,highpass=True,meas_satspot_flux=True,numthreads=self.N_threads,PSF_cube = self.PSF_cube)
                    if not os.path.exists(os.path.join(self.inputDir,"pyklip_k150a9s4m1methane_PSFsatSpotFlux")):
                        os.makedirs(os.path.join(self.inputDir,"pyklip_k150a9s4m1methane_PSFsatSpotFlux"))
                    parallelized.klip_dataset(dataset,
                                              outputdir=os.path.join(self.inputDir,"pyklip_k150a9s4m1methane_PSFsatSpotFlux"),
                                              mode="ADI+SDI",
                                              annuli=9,
                                              subsections=4,
                                              movement=1,
                                              numbasis=[20,30,50,150],
                                              spectrum="methane",
                                              fileprefix="pyklip_k150a9s4m1methane",
                                              numthreads=self.N_threads,
                                              calibrate_flux=True)


                with open(self.contrast_filename_path, 'rt') as cvs_contrast:
                    cvs_contrast_reader = csv.reader(filter(lambda row: row[0]!="#",cvs_contrast),delimiter=' ')
                    list_contrast = list(cvs_contrast_reader)
                    contrast_str_arr = np.array(list_contrast, dtype='string')
                    col_names = contrast_str_arr[0]
                    contrast_arr = contrast_str_arr[1::].astype(np.float)
                    sep_samples = contrast_arr[:,0]
                    Ttype_contrast = np.squeeze(contrast_arr[:,np.where("T-Type"==col_names)])
                    Ltype_contrast = np.squeeze(contrast_arr[:,np.where("L-Type"==col_names)])


                if not self.mute:
                    print("~~ Injecting fakes ~~")
                for pa_shift in pa_shift_list:
                    fake_flux_dict = dict(mode = "SNR",SNR=self.fakes_SNR,sep_arr = sep_samples, contrast_arr=Ttype_contrast)
                    fake_position_dict = dict(mode = "spirals",pa_shift=pa_shift)

                    # Inject the fakes
                    spdc_glob = glob(self.inputDir+os.path.sep+"S*_spdc_distorcorr.fits")
                    if overwrite_tmp or len(glob(os.path.join(self.dir_fakes,"S*_spdc_distorcorr_{0}_PA*.fits").format(self.fakes_spectrum))) != 3*len(spdc_glob):
                        if not self.mute:
                            print("~~ Reading dataset ~~")
                        dataset = GPI.GPIData(spdc_glob,highpass=True,meas_satspot_flux=True,numthreads=self.N_threads,PSF_cube = self.PSF_cube)
                        GPI.generate_spdc_with_fakes(dataset,
                                                 fake_position_dict,
                                                 fake_flux_dict,
                                                 outputdir = self.dir_fakes,
                                                 planet_spectrum = self.fakes_spectrum,
                                                 PSF_cube = self.PSF_cube_path,
                                                 star_type = None,
                                                 GOI_list_folder = self.GOI_list_folder,
                                                 mute = False,
                                                 suffix = self.fakes_spectrum+"_PA{0:02d}".format(pa_shift),
                                                 SpT_file_csv = self.GPI_TSpT_csv)

                    # Run pyklip on the fakes
                    if overwrite_tmp or len(glob(os.path.join(self.dir_fakes,"fakes_PA*_k150a9s4m1methane-KL50-speccube.fits"))) != 3:
                        # spdc_glob = glob(self.dir_fakes+os.path.sep+"S*_spdc_distorcorr*_PA{0:02d}.fits".format(pa_shift))
                        # dataset = GPI.GPIData(spdc_glob,highpass=True,meas_satspot_flux=True,numthreads=self.N_threads,PSF_cube = self.PSF_cube)
                        parallelized.klip_dataset(dataset,
                                                  outputdir=self.dir_fakes,
                                                  mode="ADI+SDI",
                                                  annuli=9,
                                                  subsections=4,
                                                  movement=1,
                                                  numbasis=[20,50,150],
                                                  spectrum="methane",
                                                  fileprefix="fakes_PA{0:02d}_k150a9s4m1methane".format(pa_shift),
                                                  numthreads=self.N_threads,
                                                  calibrate_flux=True)


            #############################
            ###### PYKLIP without sky sub MF
            self.ny_PSF = 20 # should be even
            self.nx_PSF = 20 # should be even
            # Define the cross correlation kernel
            pykliproot = os.path.dirname(os.path.realpath(parallelized.__file__))
            planet_spectrum_dir = glob(os.path.join(pykliproot,"spectra","*",self.spectrum_name+".flx"))[0]
            import pyklip.spectra_management as spec
            spectrum = spec.get_planet_spectrum(planet_spectrum_dir,"H")[1]


            #############################
            ###### PYKLIP with cross correlation            #
            # filename = os.path.join("pyklip_k150a9s4m1methane_PSFsatSpotFlux",
            #                                    "*_k150a9s4m1methane-KL50-speccube.fits")
            # inputDir_tasks.append(CrossCorr(filename,kernel_type="gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                            label="k150a9s4m1methane-KL50",mute=self.mute,kernel_width=1.0,
            #                            collapse=True,weights=spectrum,folderName=self.spectrum_name))
            # filename = os.path.join("planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"pyklip_k150a9s4m1methane-KL50-speccube-crossCorrgaussian.fits")
            # inputDir_tasks.append(StatPerPix(filename,
            #                          N_threads=self.N_threads,label="k150a9s4m1methane-KL50",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution))
            #
            #
            # err_list = kppPerDir(self.inputDir,
            #                       inputDir_tasks,
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False,outputDir=self.inputDir)
            # for err in err_list:
            #     print(err)
            #
            #
            # filename = "*_k150a9s4m1methane-KL50-speccube.fits"
            # fakesDir_tasks.append(CrossCorr(filename,kernel_type="gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                                label="k150a9s4m1methane-KL50",mute=self.mute,kernel_width=1.0,
            #                                collapse=True,weights=spectrum,folderName=self.spectrum_name))
            # filename = os.path.join("planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"fakes_PA*_k150a9s4m1methane-KL50-speccube-crossCorrgaussian.fits")
            # filename_noPlanets = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"pyklip_k150a9s4m1methane-KL50-speccube-crossCorrgaussian.fits")
            # fakesDir_tasks.append(StatPerPix(filename,filename_noPlanets=filename_noPlanets,
            #                          N_threads=self.N_threads,label="k150a9s4m1methane-KL50",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution))
            #
            #
            # err_list = kppPerDir(self.dir_fakes,
            #                       fakesDir_tasks,
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False)
            # for err in err_list:
            #     print(err)
            #
            #
            #
            # nofakes_filename = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
            #                                    "pyklip_k150a9s4m1methane-KL50-speccube-crossCorrgaussian.fits")
            # fakes_filename_list = [os.path.join(self.dir_fakes,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
            #                        "fakes_PA{0:02d}_k150a9s4m1methane-KL50-speccube-crossCorrgaussian.fits".format(pa_shift)) for pa_shift in pa_shift_list]
            # fakes_SNR_filename_list = [os.path.join(self.dir_fakes,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
            #                        "fakes_PA{0:02d}_k150a9s4m1methane-KL50-speccube-crossCorrgaussian-SNRPerPixDr2.fits".format(pa_shift)) for pa_shift in pa_shift_list]
            # separation0,contrast_curve0,throughput_tuple0 = calculate_constrat(nofakes_filename,
            #                    fakes_filename_list,
            #                    GOI_list_folder=self.GOI_list_folder,
            #                    mask_radius=self.mask_radius,IOWA=contrast_range,throughput_break=throughput_break,Dr=self.Dr,
            #                    save_dir = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name),
            #                    suffix="pyklip_crossCorr",spec_type=self.spectrum_name,fakes_SNR_filename_list=fakes_SNR_filename_list,resolution=resolution)


            #############################
            ###### PYKLIP with sky sub cross correlation

            # filename = os.path.join("pyklip_k150a9s4m1methane_PSFsatSpotFlux",
            #                                    "*_k150a9s4m1methane-KL50-speccube.fits")
            # pyklip_MFgauss = ShapeOrMF(filename,"MF","gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                            label="k150a9s4m1methane-KL50",mute=self.mute,keepPrefix=True,kernel_width=1.0,
            #                            GPI_TSpT_csv=self.GPI_TSpT_csv)
            # filename = os.path.join("planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"pyklip_k150a9s4m1methane-KL50-speccube-MF3Dgaussian.fits")
            # pyklip_SNR = StatPerPix(filename,
            #                          N_threads=self.N_threads,label="k150a9s4m1methane-KL50",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution)
            # err_list = kppPerDir(self.inputDir,
            #                       [pyklip_MFgauss,pyklip_SNR],
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False,outputDir=self.inputDir)
            # for err in err_list:
            #     print(err)
            #
            # filename = "*_k150a9s4m1methane-KL50-speccube.fits"
            # pyklip_MFgauss = ShapeOrMF(filename,"MF","gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                            label="k150a9s4m1methane-KL50",mute=self.mute,keepPrefix=True,kernel_width=1.0,
            #                            GPI_TSpT_csv=self.GPI_TSpT_csv)
            # filename = os.path.join("planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"fakes_PA*_k150a9s4m1methane-KL50-speccube-MF3Dgaussian.fits")
            # filename_noPlanets = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"pyklip_k150a9s4m1methane-KL50-speccube-MF3Dgaussian.fits")
            # pyklip_SNR = StatPerPix(filename,filename_noPlanets=filename_noPlanets,
            #                          N_threads=self.N_threads,label="k150a9s4m1methane-KL50",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution)
            # err_list = kppPerDir(self.dir_fakes,
            #                       [pyklip_MFgauss,pyklip_SNR],
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False)
            # for err in err_list:
            #     print(err)
            #
            # nofakes_filename = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
            #                                    "pyklip_k150a9s4m1methane-KL50-speccube-MF3Dgaussian.fits")
            # fakes_filename_list = [os.path.join(self.dir_fakes,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
            #                        "fakes_PA{0:02d}_k150a9s4m1methane-KL50-speccube-MF3Dgaussian.fits".format(pa_shift)) for pa_shift in pa_shift_list]
            # fakes_SNR_filename_list = [os.path.join(self.dir_fakes,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
            #                        "fakes_PA{0:02d}_k150a9s4m1methane-KL50-speccube-MF3Dgaussian-SNRPerPixDr2.fits".format(pa_shift)) for pa_shift in pa_shift_list]
            # separation1,contrast_curve1,throughput_tuple1 = calculate_constrat(nofakes_filename,
            #                    fakes_filename_list,
            #                    GOI_list_folder=self.GOI_list_folder,
            #                    mask_radius=self.mask_radius,IOWA=contrast_range,throughput_break=throughput_break,Dr=self.Dr,
            #                    save_dir = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name),
            #                    suffix="pyklip_MFgauss",spec_type=self.spectrum_name,fakes_SNR_filename_list=fakes_SNR_filename_list,resolution=resolution)



            #############################
            ###### PYKLIP with MF
            # filename = os.path.join("pyklip_k150a9s4m1methane_PSFsatSpotFlux",
            #                                    "pyklip_k150a9s4m1methane-KL50-speccube.fits")
            # pyklip_MFgauss = ShapeOrMF(filename,"shape","gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                            label="k150a9s4m1methane-KL50",mute=self.mute,keepPrefix=True,kernel_width=1.0,
            #                            GPI_TSpT_csv=self.GPI_TSpT_csv)
            # filename = os.path.join("planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"pyklip_k150a9s4m1methane-KL50-speccube-shape3Dgaussian.fits")
            # pyklip_SNR = StatPerPix(filename,
            #                          N_threads=self.N_threads,label="k150a9s4m1methane-KL50",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution)
            # err_list = kppPerDir(self.inputDir,
            #                       [pyklip_MFgauss,pyklip_SNR],
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False,outputDir=self.inputDir)
            # for err in err_list:
            #     print(err)
            #
            # filename = "*_k150a9s4m1methane-KL50-speccube.fits"
            # pyklip_SHgauss = ShapeOrMF(filename,"shape","gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                            label="k150a9s4m1methane-KL50",mute=self.mute,keepPrefix=True,kernel_width=1.0,
            #                            GPI_TSpT_csv=self.GPI_TSpT_csv)
            # filename = os.path.join("planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"fakes_PA*_k150a9s4m1methane-KL50-speccube-shape3Dgaussian.fits")
            # filename_noPlanets = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,"pyklip_k150a9s4m1methane-KL50-speccube-shape3Dgaussian.fits")
            # pyklip_SNR = StatPerPix(filename,filename_noPlanets=filename_noPlanets,
            #                          N_threads=self.N_threads,label="k150a9s4m1methane-KL50",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution)
            # err_list = kppPerDir(self.dir_fakes,
            #                       [pyklip_SHgauss,pyklip_SNR],
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False)
            # for err in err_list:
            #     print(err)
            #
            #
            nofakes_filename = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
                                               "pyklip_k150a9s4m1methane-KL50-speccube-shape3Dgaussian.fits")
            fakes_filename_list = [os.path.join(self.dir_fakes,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
                                   "fakes_PA{0:02d}_k150a9s4m1methane-KL50-speccube-shape3Dgaussian.fits".format(pa_shift)) for pa_shift in pa_shift_list]
            fakes_SNR_filename_list = [os.path.join(self.dir_fakes,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name,
                                   "fakes_PA{0:02d}_k150a9s4m1methane-KL50-speccube-shape3Dgaussian-SNRPerPixDr2.fits".format(pa_shift)) for pa_shift in pa_shift_list]
            separation2,contrast_curve2,throughput_tuple2 = calculate_constrat(nofakes_filename,
                               fakes_filename_list,
                               GOI_list_folder=self.GOI_list_folder,
                               mask_radius=self.mask_radius,IOWA=contrast_range,throughput_break=throughput_break,Dr=self.Dr,
                               save_dir = os.path.join(self.inputDir,"planet_detec_k150a9s4m1methane-KL50",self.spectrum_name),
                               suffix="pyklip_SHgauss",spec_type=self.spectrum_name,fakes_SNR_filename_list=fakes_SNR_filename_list,resolution=resolution)

            #############################
            ###### FMMF
            if 0:
                for pa_shift in pa_shift_list:
                    FMMFObj = FMMF(filename = "S*_spdc_distorcorr*_PA{0:02d}.fits".format(pa_shift),
                                    outputDir=None,
                                    N_threads=self.N_threads,
                                    predefined_sectors = "oneAc",#"oneAc",#"HR_4597",#"smallSep",
                                    label = "FMMF_PA{0:02d}".format(pa_shift),
                                    quickTest=False,
                                    overwrite=False,
                                    mute_progression = True,
                                    numbasis=[30],
                                    mvt=0.5,
                                    mvt_noTemplate=False,
                                    SpT_file_csv = self.GPI_TSpT_csv,
                                    fakes_only=True)
                    inputDir = self.dir_fakes
                    kppPerDir(inputDir,[FMMFObj],spec_path_list=[self.spectrum_name],mute_error=False)

            FMMF_metric_list = ["FMMF","FMSH","FMpF"]
            # for FMMF_metric in FMMF_metric_list:
            # #     filename = os.path.join("planet_detec_FMMF",self.spectrum_name,"*0.50-{0}.fits".format(FMMF_metric))
            # #     FMMF_SNR = StatPerPix(filename,
            # #                              N_threads=self.N_threads,label="FMMF",IOWA = None,
            # #                              type="SNR",overwrite=False,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution)
            # #     err_list = kppPerDir(self.inputDir,
            # #                           [FMMF_SNR],
            # #                           spec_path_list=[self.spectrum_name],
            # #                           mute_error = False)
            # #     for err in err_list:
            # #         print(err)
            # #     for pa_shift in pa_shift_list:
            # #         filename = os.path.join("planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,"*0.50-{0}.fits".format(FMMF_metric))
            # #         filename_noPlanets = os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name,"*0.50-{0}.fits".format(FMMF_metric))
            # #         FMMF_SNR = StatPerPix(filename,filename_noPlanets=filename_noPlanets,
            # #                                  N_threads=-1,label="FMMF_PA{0:02d}".format(pa_shift),IOWA = None,
            # #                                  type="SNR",overwrite=False,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution)
            # #         err_list = kppPerDir(self.dir_fakes,
            # #                               [FMMF_SNR],
            # #                               spec_path_list=[self.spectrum_name],
            # #                               mute_error = False)
            # #         for err in err_list:
            # #             print(err)
            # #
            #     nofakes_filename = os.path.join(os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name,
            #                                        "*_0.50-{0}.fits".format(FMMF_metric)))
            #     fakes_filename_list = [os.path.join(self.dir_fakes,"planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,
            #                                            "*_0.50-{0}.fits".format(FMMF_metric).format(pa_shift)) for pa_shift in pa_shift_list]
            #     fakes_SNR_filename_list = [os.path.join(self.dir_fakes,"planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,
            #                                            "*_0.50-{0}-SNRPerPixDr2.fits".format(FMMF_metric).format(pa_shift)) for pa_shift in pa_shift_list]
            #     separation3,contrast_curve3,throughput_tuple3 = calculate_constrat(nofakes_filename,
            #                        fakes_filename_list,
            #                        GOI_list_folder=self.GOI_list_folder,
            #                        mask_radius=self.mask_radius,IOWA=contrast_range,Dr=self.Dr,
            #                        save_dir = os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name),
            #                        suffix=FMMF_metric,spec_type=self.spectrum_name,fakes_SNR_filename_list=fakes_SNR_filename_list,resolution=resolution)


            #############################
            ###### PYKLIP (fm.py) with cross correlation
            # # Define the cross correlation kernel
            # pykliproot = os.path.dirname(os.path.realpath(parallelized.__file__))
            # planet_spectrum_dir = glob(os.path.join(pykliproot,"spectra","*",self.spectrum_name+".flx"))[0]
            # import pyklip.spectra_management as spec
            # spectrum = spec.get_planet_spectrum(planet_spectrum_dir,"H")[1]
            #
            # filename = os.path.join("planet_detec_FMMF",self.spectrum_name,
            #                                    "*_t600g100nc_0.50-speccube-KL30.fits")
            # inputDir_tasks.append(CrossCorr(filename,kernel_type="gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                            label="FMMF",mute=self.mute,kernel_width=1.0,
            #                            collapse=True,weights=spectrum,folderName=self.spectrum_name))
            # filename = os.path.join("planet_detec_FMMF",self.spectrum_name,
            #                                    "*_t600g100nc_0.50-speccube-KL30-crossCorrgaussian.fits")
            # inputDir_tasks.append(StatPerPix(filename,
            #                          N_threads=self.N_threads,label="FMMF",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution))
            #
            #
            # err_list = kppPerDir(self.inputDir,
            #                       inputDir_tasks,
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False)
            # for err in err_list:
            #     print(err)
            #
            #
            # for pa_shift in pa_shift_list:
            #     fakesDir_tasks = []
            #     filename = os.path.join("planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,
            #                                        "*_t600g100nc_0.50-speccube-KL30.fits")
            #     fakesDir_tasks.append(CrossCorr(filename,kernel_type="gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                                    label="FMMF_PA{0:02d}".format(pa_shift),mute=self.mute,kernel_width=1.0,
            #                                    collapse=True,weights=spectrum,folderName=self.spectrum_name))
            #     filename = os.path.join("planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,
            #                                        "*_t600g100nc_0.50-speccube-KL30-crossCorrgaussian.fits")
            #     filename_noPlanets = os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name,"*_t600g100nc_0.50-speccube-KL30-crossCorrgaussian.fits")
            #     fakesDir_tasks.append(StatPerPix(filename,filename_noPlanets=filename_noPlanets,
            #                              N_threads=self.N_threads,label="FMMF_PA{0:02d}".format(pa_shift),IOWA = None,
            #                              type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution))
            #
            #
            #     err_list = kppPerDir(self.dir_fakes,
            #                           fakesDir_tasks,
            #                           spec_path_list=[self.spectrum_name],
            #                           mute_error = False)
            #     for err in err_list:
            #         print(err)
            #
            #
            #
            # nofakes_filename = os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name,
            #                                    "*_t600g100nc_0.50-speccube-KL30-crossCorrgaussian.fits")
            # fakes_filename_list = [os.path.join(self.dir_fakes,"planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,
            #                        "*_t600g100nc_0.50-speccube-KL30-crossCorrgaussian.fits") for pa_shift in pa_shift_list]
            # fakes_SNR_filename_list = [os.path.join(self.dir_fakes,"planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,
            #                        "*_t600g100nc_0.50-speccube-KL30-crossCorrgaussian-SNRPerPixDr2.fits") for pa_shift in pa_shift_list]
            # separation0,contrast_curve0,throughput_tuple0 = calculate_constrat(nofakes_filename,
            #                    fakes_filename_list,
            #                    GOI_list_folder=self.GOI_list_folder,
            #                    mask_radius=self.mask_radius,IOWA=contrast_range,throughput_break=None,Dr=self.Dr,
            #                    save_dir = os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name),
            #                    suffix="fmpyklip_crossCorr",spec_type=self.spectrum_name,fakes_SNR_filename_list=fakes_SNR_filename_list,resolution=resolution)


            #############################
            ###### PYKLIP (fm.py) with MF
            # filename = os.path.join("planet_detec_FMMF",self.spectrum_name,
            #                                    "*_t600g100nc_0.50-speccube-KL30.fits")
            # pyklip_MFgauss = ShapeOrMF(filename,"shape","gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                            label="FMMF",mute=self.mute,keepPrefix=True,kernel_width=1.0,
            #                            GPI_TSpT_csv=self.GPI_TSpT_csv)
            # filename = os.path.join("planet_detec_FMMF",self.spectrum_name,
            #                                    "*_t600g100nc_0.50-speccube-KL30-shape3Dgaussian.fits")
            # pyklip_SNR = StatPerPix(filename,
            #                          N_threads=self.N_threads,label="FMMF",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution)
            # err_list = kppPerDir(self.inputDir,
            #                       [pyklip_MFgauss,pyklip_SNR],
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False)
            # for err in err_list:
            #     print(err)
            #
            # filename = os.path.join("planet_detec_FMMF_PA??",self.spectrum_name,
            #                                    "*_t600g100nc_0.50-speccube-KL30.fits")
            # pyklip_SHgauss = ShapeOrMF(filename,"shape","gaussian",N_threads=self.N_threads,overwrite=overwrite_tmp,
            #                            label="FMMF",mute=self.mute,keepPrefix=True,kernel_width=1.0,
            #                            GPI_TSpT_csv=self.GPI_TSpT_csv)
            # filename = os.path.join("planet_detec_FMMF_PA??",self.spectrum_name,"*_t600g100nc_0.50-speccube-KL30-shape3Dgaussian.fits")
            # filename_noPlanets = os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name,"*_t600g100nc_0.50-speccube-KL30-shape3Dgaussian.fits")
            # pyklip_SNR = StatPerPix(filename,filename_noPlanets=filename_noPlanets,
            #                          N_threads=self.N_threads,label="FMMF",IOWA = None,
            #                          type="SNR",overwrite=overwrite_tmp,GOI_list_folder=self.GOI_list_folder,mute=self.mute,resolution=resolution)
            # err_list = kppPerDir(self.dir_fakes,
            #                       [pyklip_SHgauss,pyklip_SNR],
            #                       spec_path_list=[self.spectrum_name],
            #                       mute_error = False)
            # for err in err_list:
            #     print(err)
            #
            #
            # nofakes_filename = os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name,
            #                                    "*_t600g100nc_0.50-speccube-KL30-shape3Dgaussian.fits")
            # fakes_filename_list = [os.path.join(self.dir_fakes,"planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,
            #                        "*_t600g100nc_0.50-speccube-KL30-shape3Dgaussian.fits") for pa_shift in pa_shift_list]
            # fakes_SNR_filename_list = [os.path.join(self.dir_fakes,"planet_detec_FMMF_PA{0:02d}".format(pa_shift),self.spectrum_name,
            #                        "*_t600g100nc_0.50-speccube-KL30-shape3Dgaussian-SNRPerPixDr2.fits") for pa_shift in pa_shift_list]
            # separation2,contrast_curve2,throughput_tuple2 = calculate_constrat(nofakes_filename,
            #                    fakes_filename_list,
            #                    GOI_list_folder=self.GOI_list_folder,
            #                    mask_radius=self.mask_radius,IOWA=contrast_range,throughput_break=None,Dr=self.Dr,
            #                    save_dir = os.path.join(self.inputDir,"planet_detec_FMMF",self.spectrum_name),
            #                    suffix="fmpyklip_SHgauss",spec_type=self.spectrum_name,fakes_SNR_filename_list=fakes_SNR_filename_list,resolution=resolution)


        # Removing fakes
        spdc_glob = glob(self.dir_fakes+os.path.sep+"S*_spdc_distorcorr*_PA*.fits")
        for filename in spdc_glob:
            print("Removing {0}".format(filename))
            os.remove(filename)

        if not self.reduce_only:

            with open(self.contrast_filename_path, 'rt') as cvs_contrast:
                cvs_contrast_reader = csv.reader(filter(lambda row: row[0]!="#",cvs_contrast),delimiter=' ')
                list_contrast = list(cvs_contrast_reader)
                contrast_str_arr = np.array(list_contrast, dtype='string')
                col_names = contrast_str_arr[0]
                contrast_arr = contrast_str_arr[1::].astype(np.float)
                sep_samples = contrast_arr[:,0]
                Ttype_contrast = np.squeeze(contrast_arr[:,np.where("T-Type"==col_names)])
                Ltype_contrast = np.squeeze(contrast_arr[:,np.where("L-Type"==col_names)])

            import matplotlib.pyplot as plt
            #############################
            ###### FINAL CONTRAST PLOT
            legend_str_list = []
            plt.figure(1,figsize=(8,6))
            plt.plot(sep_samples,Ttype_contrast,"--", color='b', linewidth=3.0)
            legend_str_list.append("Jason T-Type pyklip")
            plt.plot(sep_samples,Ltype_contrast,"--", color='r', linewidth=3.0)
            legend_str_list.append("Jason L-type pyklip")
            # plt.plot(self.pyklip_noSky_stddev_rSamp*0.01413,5*self.pyklip_noSky_1Dstddev/pyklip_noSky_throughput_func(self.pyklip_noSky_stddev_rSamp*0.01413),":", color='b', linewidth=3.0)
            # legend_str_list.append("JB's T-Type pyklip no sky sub")

            suffix_list = ["pyklip_crossCorr","pyklip_MFgauss","pyklip_SHgauss","FMSH","FMMF","FMpF","fmpyklip_crossCorr","fmpyklip_SHgauss"]
            planet_detec_list = ["planet_detec_k150a9s4m1methane-KL50","planet_detec_k150a9s4m1methane-KL50","planet_detec_k150a9s4m1methane-KL50","planet_detec_FMMF","planet_detec_FMMF","planet_detec_FMMF","planet_detec_FMMF","planet_detec_FMMF"]
            linestyle_list = [":","-","-","-.","-.","-.","-.","-."]
            color_list = ["blue","purple","cyan","yellow","orange","red","blue","cyan"]
            for suffix,linestyle,color,planet_detec in zip(suffix_list,linestyle_list,color_list,planet_detec_list):
                with open(os.path.join(self.inputDir,planet_detec,self.spectrum_name,"contrast-"+suffix+'.csv'), 'rt') as cvs_contrast:
                    cvs_contrast_reader = csv.reader(filter(lambda row: row[0]!="#",cvs_contrast),delimiter=' ')
                    list_contrast = list(cvs_contrast_reader)
                    contrast_str_arr = np.array(list_contrast, dtype='string')
                    col_names = contrast_str_arr[0]
                    contrast_arr = contrast_str_arr[1::].astype(np.float)
                    sep_samples = contrast_arr[:,0]
                    Ttype_contrast = np.squeeze(contrast_arr[:,np.where(self.spectrum_name==col_names)])


                plt.plot(sep_samples,Ttype_contrast,linestyle=linestyle, color=color, linewidth=3.0)
                legend_str_list.append("JB's T-Type "+suffix)

            plt.xlabel("Separation (arcsec)", fontsize=20)
            plt.ylabel("Contrast", fontsize=20)
            plt.legend(legend_str_list)
            ax= plt.gca()
            ax.set_yscale('log')
            ax.tick_params(axis='x', labelsize=20)
            ax.tick_params(axis='y', labelsize=20)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            plt.grid(True)
            # plt.ylim([0,1.2])
            # plt.xlim([0,4])
            if not self.mute:
                print("Saving: "+os.path.join(self.save_contrast,self.star_name+"_"+self.compact_date+"_"+self.filter+"_contrast-all.png"))
            plt.savefig(os.path.join(self.save_contrast,self.star_name+"_"+self.compact_date+"_"+self.filter+"_contrast-all.png"), bbox_inches='tight')


            plt.figure(2,figsize=(20,10))
            # suffix_list = ["pyklip_crossCorr","pyklip_MFgauss","pyklip_SHgauss","FMSH","FMMF","FMpF"]
            # linestyle_list = [":","-","-","-.","-.","-."]
            # color_list = ["blue","purple","cyan","yellow","orange","red"]
            title_list = suffix_list
            # suffix_list = ["FMSH","fmpyklip_SHgauss","fmpyklip_crossCorr"]
            # title_list = ["Forward Model Matched Filter","Gaussian Matched Filter","Gaussian Cross Correlation"]
            # planet_detec_list = ["planet_detec_FMMF","planet_detec_FMMF","planet_detec_FMMF"]
            # linestyle_list = ["-","--",":"]
            # color_list = ["#ff9900","#006699","grey"]
            for k,(suffix,linestyle,color,planet_detec,title_str) in enumerate(zip(suffix_list,linestyle_list,color_list,planet_detec_list,title_list)):
                # plt.subplot(2,4,k+1)
                plt.subplot(2,4,k+1)
                with open(os.path.join(self.inputDir,planet_detec,self.spectrum_name,"throughput-"+suffix+'.csv'), 'rt') as cvs_throughput:
                    cvs_reader = csv.reader(filter(lambda row: row[0]!="#",cvs_throughput),delimiter=' ')
                    list_contrast = list(cvs_reader)
                    contrast_str_arr = np.array(list_contrast, dtype='string')
                    col_names = contrast_str_arr[0]
                    contrast_arr = contrast_str_arr[1::].astype(np.float)
                    sep_samples = contrast_arr[:,0]
                    Ttype_throughput = np.squeeze(contrast_arr[:,np.where("throughput"==col_names)])
                    Ttype_throughput_fit = np.squeeze(contrast_arr[:,np.where("fit"==col_names)])


                ax = plt.gca()
                plt.title(title_str,fontsize=20)
                # plt.text(0.5,1.1,title_str,horizontalalignment='center',fontsize=20,transform=ax.transAxes)
                plt.scatter(sep_samples,Ttype_throughput,c=color)
                sep_samples_unique,sep_indices = np.unique(sep_samples,return_index=True)
                plt.plot(sep_samples_unique,Ttype_throughput_fit[sep_indices],linestyle="-",linewidth=5,color=color)
                plt.xlabel("Separation (as)", fontsize=20)
                plt.ylabel("Conversion Factor", fontsize=20)
                plt.xlim([0,1.8])
                plt.xticks([0,0.5,1.0,1.5])
                # plt.xlim([0,1.])
                # plt.ylim([0,None])
                ax= plt.gca()
                plt.grid(True)
                ax.tick_params(axis='x', labelsize=20)
                ax.tick_params(axis='y', labelsize=20)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')
                plt.ticklabel_format(style='sci',axis='y',scilimits=(0,1))

            plt.tight_layout()

            if not self.mute:
                print("Saving: "+os.path.join(self.save_contrast,self.star_name+"_"+self.compact_date+"_"+self.filter+"_throughput-all.png"))
            plt.savefig(os.path.join(self.save_contrast,self.star_name+"_"+self.compact_date+"_"+self.filter+"_throughput-all.png"), bbox_inches='tight')

            plt.figure(3,figsize=(15,10))
            for k,(suffix,linestyle,color,planet_detec) in enumerate(zip(suffix_list,linestyle_list,color_list,planet_detec_list)):
                plt.subplot(2,4,k+1)
                with open(os.path.join(self.inputDir,planet_detec,self.spectrum_name,"contrast-SNR-check-"+suffix+'.csv'), 'rt') as cvs_SNR:
                    cvs_reader = csv.reader(filter(lambda row: row[0]!="#",cvs_SNR),delimiter=' ')
                    list_contrast = list(cvs_reader)
                    contrast_str_arr = np.array(list_contrast, dtype='string')
                    col_names = contrast_str_arr[0]
                    contrast_arr = contrast_str_arr[1::].astype(np.float)
                    sep_samples = contrast_arr[:,0]
                    SNR_from_contrast = np.squeeze(contrast_arr[:,np.where("ContSNR"==col_names)])
                    SNR_fakes = np.squeeze(contrast_arr[:,np.where("PixSNR"==col_names)])


                plt.title(suffix)
                plt.plot([0,20],[0,20],"r-",linewidth=5)
                plt.ylim([0,20])
                plt.xlim([0,20])
                plt.scatter(SNR_from_contrast,SNR_fakes)
                plt.xlabel("SNR from contrast curve", fontsize=20)
                plt.ylabel("SNR from image", fontsize=20)
                ax= plt.gca()
                ax.tick_params(axis='x', labelsize=20)
                ax.tick_params(axis='y', labelsize=20)

            plt.tight_layout()
            if not self.mute:
                print("Saving: "+os.path.join(self.save_contrast,self.star_name+"_"+self.compact_date+"_"+self.filter+"_contrast-SNR-check-all.png"))
            plt.savefig(os.path.join(self.save_contrast,self.star_name+"_"+self.compact_date+"_"+self.filter+"_contrast-SNR-check-all.png"), bbox_inches='tight')

            # plt.show()


        # self.flux_1Dstddev_map = get_image_stat_map(self.image,
        #                                             self.image_without_planet,
        #                                             IOWA = self.IOWA,
        #                                             N = None,
        #                                             centroid = self.center,
        #                                             r_step = self.Dr/2,
        #                                             Dr = self.Dr,
        #                                             type = "stddev",
        #                                             image_wide = None)
        #
        #
        # self.fluxMap_stddev = get_image_stat_map_perPixMasking(self.image,
        #                                                  self.image_without_planet,
        #                                                  mask_radius = self.mask_radius,
        #                                                  IOWA = self.IOWA,
        #                                                  N = self.N,
        #                                                  centroid = self.center,
        #                                                  mute = self.mute,
        #                                                  N_threads = self.N_threads,
        #                                                  Dr= self.Dr,
        #                                                  Dth = None,
        #                                                  type = self.type)

        return None

    def save(self):
        """

        :return: None
        """

        if not os.path.exists(self.outputDir+os.path.sep+self.folderName):
            os.makedirs(self.outputDir+os.path.sep+self.folderName)


        return None

    def load(self):
        """

        :return: None
        """

        return None


def gather_contrasts(base_dir,filename_filter_list,mute = False,epoch_suffix=None):
    """
    Build the multiple combined ROC curve from individual frame ROC curve while making sure they have the same inputs.
    If the folders are organized following the convention below then it will make sure there is a ROC file for each
    filename_filter in each epoch. Otherwise it skips the epoch.

    The folders need to be organized as:
     base_dir/TARGET/autoreduced/EPOCH_Spec/filename_filter

    In the function TARGET and EPOCH are wild characters.

    It looks for all the file matching filename_filter using glob.glob and then add each individual ROC to build the
    master ROC.

    Plot master_N_false_pos vs master_N_true_detec to get a ROC curve.

    :param base_dir: Base directory from which the file search go.
    :param filename_filter: Filename filter with wild characters indicating which files to pick.
    :param mute: If True, mute prints. Default is False.
    :return: threshold_sampling,master_N_false_pos,master_N_true_detec:
        threshold_sampling: The metric sampling. It is the curve parametrization.
        master_N_false_pos: Number of false positives as a function of threshold_sampling
        master_N_true_detec: Number of true positives as a function of threshold_sampling
    """

    N_CONT = len(filename_filter_list)
    sep_samp_list = [[]]*N_CONT
    cont_list = [[]]*N_CONT

    if epoch_suffix is None:
        epoch_suffix = ""

    dirs_to_reduce = os.listdir(base_dir)
    N=0
    for object in dirs_to_reduce:
        if not object.startswith('.'):
            #print(object)

            epochDir_glob = glob(base_dir+object+os.path.sep+"autoreduced"+os.path.sep+"*_*_Spec"+epoch_suffix+os.path.sep)

            for epochDir in epochDir_glob:
                inputDir = os.path.abspath(epochDir)

                file_list = []
                for filename_filter in filename_filter_list:
                    try:
                        print(inputDir+os.path.sep+filename_filter)
                        print(glob(inputDir+os.path.sep+filename_filter)[0])
                        file_list.append(glob(inputDir+os.path.sep+filename_filter)[0])
                        # if not mute:
                        #     print("ROC: {0} in {1}. Adding.".format(filename_filter,inputDir))
                    except:
                        pass
                        # if not mute:
                        #     print("ROC: {0} unvailable in {1}. Skipping".format(filename_filter,inputDir))

                if len(file_list) == N_CONT:
                    print(file_list)
                    N=N+1
                    for index,filename in enumerate(file_list):
                        with open(filename, 'rb') as csvfile:
                            # reader = csv.reader(csvfile, delimiter=' ')
                            # csv_as_list = list(reader)
                            # detec_table_labels = csv_as_list[0]
                            # detec_table = np.array(csv_as_list[1::], dtype='string').astype(np.float)

                            cvs_contrast_reader = csv.reader(filter(lambda row: row[0]!="#",csvfile),delimiter=' ')
                            list_contrast = list(cvs_contrast_reader)
                            # print(list_contrast)
                            contrast_str_arr = np.array(list_contrast, dtype='string')
                            # print(contrast_str_arr)
                            col_names = contrast_str_arr[0]
                            contrast_arr = contrast_str_arr[1::].astype(np.float)
                            sep_samples = contrast_arr[:,0]

                            methane_idx = np.where("T-Type"==col_names)[0]

                            methane_contrast = np.squeeze(contrast_arr[:,methane_idx])

                        try:
                            cont_list[index] = cont_list[index]+methane_contrast
                        except:
                            sep_samp_list[index] = sep_samples
                            cont_list[index] = methane_contrast

    print("N files = {0}".format(N))

    return sep_samp_list,np.array(cont_list)/N

def calculate_constrat(nofakes_filename,fakes_filename_list,
                       GOI_list_folder=None,
                       mask_radius=None,
                       IOWA=None,
                       Dr=None,
                       save_dir = None,
                       suffix=None,
                       spec_type=None,
                       fakes_SNR_filename_list=None,
                       resolution=None,
                       conversion_break = None,
                       linfit = False):
    '''

    :param nofakes_filename:
    :param fakes_filename_list:
    :param GOI_list_folder:
    :param mask_radius:
    :param IOWA:
    :param Dr:
    :return:
    '''

    hdulist = pyfits.open(glob(nofakes_filename)[0])
    metric_image = hdulist[1].data
    exthdr = hdulist[1].header
    prihdr = hdulist[0].header
    center = [exthdr['PSFCENTX'], exthdr['PSFCENTY']]
    if IOWA is None:
        IWA,OWA,inner_mask,outer_mask = get_occ(metric_image, centroid = center)
        IOWA = (IWA,OWA)
        IOWA_as = (pix2as(IWA),pix2as(OWA))

    real_contrast_list = []
    sep_list = []
    pa_list = []
    metric_fakes_val = []
    metric_nofakes_val = []

    for fakes_filename in fakes_filename_list:
        # get conversion for pyklip images
        hdulist = pyfits.open(glob(fakes_filename)[0])
        metric_image_fakes = hdulist[1].data
        exthdr_fakes = hdulist[1].header
        prihdr_fakes = hdulist[0].header

        row_real_object_list,col_real_object_list = get_pos_known_objects(prihdr_fakes,exthdr_fakes,fakes_only=True)
        sep,pa = get_pos_known_objects(prihdr_fakes,exthdr_fakes,pa_sep=True,fakes_only=True)
        sep_list.extend(sep)
        pa_list.extend(pa)
        for fake_id in range(100):
            try:
                real_contrast_list.append(exthdr_fakes["FKCONT{0:02d}".format(fake_id)])
            except:
                continue
        # print([(metric_image_fakes[np.round(row_real_object),np.round(col_real_object)],
        #         np.nanmax(metric_image_fakes[(np.round(row_real_object)-1):(np.round(row_real_object)+2),(np.round(col_real_object)-1):(np.round(col_real_object)+2)]))\
        #                              for row_real_object,col_real_object in zip(row_real_object_list,col_real_object_list)])
        # metric_fakes_val.extend([np.nanmax(metric_image_fakes[(np.round(row_real_object)-1):(np.round(row_real_object)+2),(np.round(col_real_object)-1):(np.round(col_real_object)+2)]) \
        #                              for row_real_object,col_real_object in zip(row_real_object_list,col_real_object_list)])
        metric_fakes_val.extend([metric_image_fakes[int(np.round(row_real_object)),int(np.round(col_real_object))] \
                                     for row_real_object,col_real_object in zip(row_real_object_list,col_real_object_list)])
        metric_nofakes_val.extend([metric_image[int(np.round(row_real_object)),int(np.round(col_real_object))] \
                                     for row_real_object,col_real_object in zip(row_real_object_list,col_real_object_list)])

    metric_fakes_val =  np.array(metric_fakes_val) - np.array(metric_nofakes_val)

    if GOI_list_folder is not None:
        metric_image_without_planet = mask_known_objects(metric_image,prihdr,exthdr,GOI_list_folder, mask_radius = mask_radius)
        # metric_image_without_planet = mask_known_objects(metric_image_fakes,prihdr_fakes,exthdr_fakes,GOI_list_folder, mask_radius = mask_radius)
    else:
        metric_image_without_planet = metric_image

    metric_1Dstddev,metric_stddev_rSamp = get_image_stddev(metric_image_without_planet,
                                                                 IOWA,
                                                                 N = None,
                                                                 centroid = center,
                                                                 r_step = Dr/2,
                                                                 Dr=Dr,
                                                                 resolution=resolution)
    metric_stddev_rSamp = np.array([r_tuple[0] for r_tuple in metric_stddev_rSamp])
    metric_1Dstddev = np.array(metric_1Dstddev)
    from  scipy.interpolate import interp1d
    metric_1Dstddev_func = interp1d(metric_stddev_rSamp,metric_1Dstddev,bounds_error=False, fill_value=np.nan)

    whereNoNans = np.where(np.isfinite(metric_fakes_val))
    metric_fakes_val = metric_fakes_val[whereNoNans]
    sep_list = np.array(sep_list)[whereNoNans]
    pa_list = np.array(pa_list)[whereNoNans]
    real_contrast_list =  np.array(real_contrast_list)[whereNoNans]

    sep_list,pa_list,metric_fakes_val,real_contrast_list = zip(*sorted(zip(sep_list,pa_list,metric_fakes_val,real_contrast_list)))
    metric_fakes_val = np.array(metric_fakes_val)
    sep_list =  np.array(sep_list)
    pa_list =  np.array(pa_list)
    real_contrast_list =  np.array(real_contrast_list)
    if linfit:
        if conversion_break is not None:
            whereInRange = np.where((np.array(sep_list)>IOWA_as[0])*(np.array(sep_list)<conversion_break))
            z1 = np.polyfit(np.array(sep_list)[whereInRange],np.array(real_contrast_list)[whereInRange]/np.array(metric_fakes_val)[whereInRange],1)

            whereInRange = np.where((np.array(sep_list)>conversion_break)*(np.array(sep_list)<IOWA_as[1]))
            z2 = np.polyfit(np.array(sep_list)[whereInRange],np.array(real_contrast_list)[whereInRange]/np.array(metric_fakes_val)[whereInRange],1)

            linfit1 = np.poly1d(z1)
            linfit2 = np.poly1d(z2)
            metric_conversion_func = lambda sep: np.concatenate((linfit1(np.array(sep)[np.where(np.array(sep)<conversion_break)]),
                                                                 linfit2(np.array(sep)[np.where(conversion_break<np.array(sep))])))

        else:
            whereInRange = np.where((np.array(sep_list)>IOWA_as[0])*(np.array(sep_list)<IOWA_as[1]))
            z = np.polyfit(np.array(sep_list)[whereInRange],np.array(real_contrast_list)[whereInRange]/np.array(metric_fakes_val)[whereInRange],1)
            metric_conversion_func = np.poly1d(z)
    else:
        whereInRange = np.where((sep_list>IOWA_as[0])*(sep_list<IOWA_as[1]))
        metric_fakes_in_range = metric_fakes_val[whereInRange]
        cont_in_range = real_contrast_list[whereInRange]
        unique_sep = np.unique(sep_list[whereInRange])
        med_conversion = np.zeros(len(unique_sep))
        std_conversion = np.zeros(len(unique_sep))
        for k,sep_it in enumerate(unique_sep):
            where_sep = np.where(sep_list==sep_it)
            # med_conversion[k] = np.nanmedian(cont_in_range[where_sep]) / np.nanmedian(metric_fakes_in_range[where_sep])
            # std_conversion[k] = np.nanmedian(np.abs(cont_in_range[where_sep]/metric_fakes_in_range[where_sep] - med_conversion[k]))
            med_conversion[k] = np.nanmean(cont_in_range[where_sep]/metric_fakes_in_range[where_sep])
            var_conversion = np.nanvar(cont_in_range[where_sep] / metric_fakes_in_range[where_sep])
            std_conversion[k] = np.sqrt(var_conversion)
        metric_conversion_func = interp1d(unique_sep,med_conversion,bounds_error=False, fill_value=np.nan)
        # std_conversion_func = interp1d(unique_sep,1.4826*std_conversion,bounds_error=False, fill_value=np.nan)
        std_conversion_func = interp1d(unique_sep,std_conversion,bounds_error=False, fill_value=np.nan)


    if 0:
        import matplotlib.pyplot as plt
        plt.figure(2)
        plt.title(suffix)
        plt.plot(sep_list,np.array(metric_fakes_val)/np.array(real_contrast_list),"*")
        plt.plot(sep_list,metric_conversion_func(sep_list),"-")
        plt.xlabel("Separation (arcsec)", fontsize=20)
        plt.ylabel("Throughput (arbritrary units)", fontsize=20)
        ax= plt.gca()
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        plt.show()


    contrast_curve = 5*metric_1Dstddev*metric_conversion_func(pix2as(metric_stddev_rSamp))

    if fakes_SNR_filename_list is not None:
        SNR_real_contrast_list = []
        SNR_sep_list = []
        SNR_fakes=[]

        for fakes_SNR_filename in fakes_SNR_filename_list:
            # get conversion for pyklip images
            hdulist = pyfits.open(glob(fakes_SNR_filename)[0])
            SNR_map_fakes = hdulist[1].data
            exthdr_fakes_SNR = hdulist[1].header
            prihdr_fakes_SNR = hdulist[0].header

            row_real_object_list,col_real_object_list = get_pos_known_objects(prihdr_fakes_SNR,exthdr_fakes_SNR,fakes_only=True)
            sep,pa_real_object_list = get_pos_known_objects(prihdr_fakes_SNR,exthdr_fakes_SNR,pa_sep=True,fakes_only=True)
            SNR_sep_list.extend(sep)
            for fake_id in range(100):
                try:
                    SNR_real_contrast_list.append(exthdr_fakes_SNR["FKCONT{0:02d}".format(fake_id)])
                except:
                    continue
            SNR_fakes.extend([np.nanmax(SNR_map_fakes[(int(np.round(row_real_object))-1):(int(np.round(row_real_object))+2),(int(np.round(col_real_object))-1):(int(np.round(col_real_object))+2)]) \
                                         for row_real_object,col_real_object in zip(row_real_object_list,col_real_object_list)])
            # SNR_fakes.extend([SNR_map_fakes[np.round(row_real_object),np.round(col_real_object)] \
            #                              for row_real_object,col_real_object in zip(row_real_object_list,col_real_object_list)])

        from scipy.interpolate import interp1d
        contrast_curve_interp = interp1d(pix2as(metric_stddev_rSamp),contrast_curve,kind="linear",bounds_error=False)
        SNR_from_contrast = np.array(SNR_real_contrast_list)/(contrast_curve_interp(SNR_sep_list)/5.0)

        with open(os.path.join(save_dir,"contrast-SNR-check-"+suffix+'.csv'), 'w+') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=' ')
            csvwriter.writerows([["Seps","ContSNR","PixSNR","contrast"]])
            csvwriter.writerows(zip(SNR_sep_list,SNR_from_contrast,SNR_fakes,SNR_real_contrast_list))



    if save_dir is not None:
        if suffix is None:
            suffix = "default"

        with open(os.path.join(save_dir,"contrast-"+suffix+'.csv'), 'w+') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=' ')
            csvwriter.writerows([["Seps",spec_type,spec_type+"_met_std",spec_type+"_conv",spec_type+"_conv_std"]])
            contrast_curve[np.where(np.isnan(contrast_curve))] = -1.
            not_neg = np.where(contrast_curve>0)
            csvwriter.writerows(zip(pix2as(metric_stddev_rSamp[not_neg]),
                                    contrast_curve[not_neg],
                                    metric_1Dstddev[not_neg],
                                    metric_conversion_func(pix2as(metric_stddev_rSamp[not_neg])),
                                    std_conversion_func(pix2as(metric_stddev_rSamp[not_neg]))))

        with open(os.path.join(save_dir,"conversion-"+suffix+'.csv'), 'w+') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=' ')
            csvwriter.writerows([["Seps","PA","conversion","fit","metric","contrast","kMAD"]])
            csvwriter.writerows(zip(sep_list,pa_list,
                                    np.array(real_contrast_list)/np.array(metric_fakes_val),
                                    metric_conversion_func(sep_list),
                                    np.array(metric_fakes_val),
                                    np.array(real_contrast_list),
                                    std_conversion_func(sep_list)))

    conversion_tuple = (sep_list,np.array(metric_fakes_val)/np.array(real_contrast_list),np.array(metric_fakes_val),np.array(real_contrast_list))

    return pix2as(metric_stddev_rSamp),contrast_curve,conversion_tuple

