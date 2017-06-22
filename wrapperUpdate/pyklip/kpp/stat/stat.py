__author__ = 'JB'
import os
import astropy.io.fits as pyfits
from glob import glob
import multiprocessing as mp
import numpy as np
from scipy.signal import correlate2d

from pyklip.kpp.utils.kppSuperClass import KPPSuperClass
from pyklip.kpp.stat.stat_utils import *
from pyklip.kpp.utils.GOI import *
import pyklip.kpp.utils.mathfunc as kppmath

class Stat(KPPSuperClass):
    """
    Class for SNR calculation.
    """
    def __init__(self,filename,
                 filename_noPlanets = None,
                 inputDir = None,
                 outputDir = None,
                 folderName = None,
                 mute=None,
                 N_threads=None,
                 label = None,
                 mask_radius = None,
                 IOWA = None,
                 N = None,
                 Dr = None,
                 r_step = None,
                 type = None,
                 rm_edge = None,
                 GOI_list_folder = None,
                 overwrite = False,
                 kernel_type = None,
                 kernel_width = None,
                 image_wide = None,
                 collapse = None,
                 weights = None,
                 nans2zero = None):
        """


        :param filename: Filename of the file on which to calculate the metric. It should be the complete path unless
                        inputDir is defined.
                        It can include wild characters. The file will be selected using the first output of glob.glob().
        :param filename_noPlanets: One should be careful with this one since it requires it should find the same number
                            of files with no signal than normal images when calling glob.glob().
                            Besides one has to check that the ordering of glob outputs are matching for both lists.
        :param mute: If True prevent printed log outputs.
        :param N_threads: Number of threads to be used for the metrics and the probability calculations.
                        If None use mp.cpu_count().
                        If -1 do it sequentially.
                        Note that it is not used for this super class.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        :param collapse: If true and input is 3D then it will collapse the final map.
        :param weights: If not None and collapse is True then a weighted mean is performed using the weights.
        """
        # allocate super class
        super(Stat, self).__init__(filename,
                                     inputDir = inputDir,
                                     outputDir = outputDir,
                                     folderName = folderName,
                                     mute=mute,
                                     N_threads=N_threads,
                                     label=label,
                                     overwrite = overwrite)

        if mask_radius is None:
            self.mask_radius = 7
        else:
            self.mask_radius = mask_radius

        if Dr is None:
            self.Dr = 2
        else:
            self.Dr = Dr


        if r_step is None:
            self.r_step = 2
        else:
            self.r_step = r_step

        if type is None:
            self.type = "SNR"
        else:
            self.type = type

        self.IOWA = IOWA
        self.N = N
        self.rm_edge = rm_edge
        self.GOI_list_folder = GOI_list_folder
        self.filename_noPlanets = filename_noPlanets
        self.image_wide = image_wide

        self.kernel_type = kernel_type
        # The default value is defined later
        self.kernel_width = kernel_width

        if collapse is None:
            self.collapse = False
        else:
            self.collapse = collapse
        self.weights = weights
        if nans2zero is None:
            self.nans2zero = False
        else:
            self.nans2zero = nans2zero

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
        init_out = super(Stat, self).initialize(inputDir = inputDir,
                                         outputDir = outputDir,
                                         folderName = folderName,
                                         label=label)



        if self.filename_noPlanets is not None:# Check file existence and define filename_path
            if self.inputDir is None or os.path.isabs(self.filename_noPlanets):
                try:
                    if len(glob(self.filename_noPlanets)) == self.N_matching_files:
                        self.filename_noPlanets_path = os.path.abspath(glob(self.filename_noPlanets)[self.id_matching_file-1])
                    else:
                        self.filename_noPlanets_path = os.path.abspath(glob(self.filename_noPlanets)[0])
                except:
                    raise Exception("File "+self.filename_noPlanets+"doesn't exist.")
            else:
                try:
                    if len(glob(self.inputDir+os.path.sep+self.filename_noPlanets)) == self.N_matching_files:
                        self.filename_noPlanets_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_noPlanets)[self.id_matching_file-1])
                    else:
                        self.filename_noPlanets_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_noPlanets)[0])
                except:
                    raise Exception("File "+self.inputDir+os.path.sep+self.filename_noPlanets+" doesn't exist.")

            # Open the fits file on which the metric will be applied
            hdulist = pyfits.open(self.filename_noPlanets_path)
            if not self.mute:
                print("Opened: "+self.filename_noPlanets_path)

            # grab the data and headers
            try:
                self.image_noPlanets = hdulist[1].data
                self.exthdr_noPlanets = hdulist[1].header
                self.prihdr_noPlanets = hdulist[0].header
            except:
                # This except was used for datacube not following GPI headers convention.
                if not self.mute:
                    print("Couldn't read the fits file with GPI conventions. Try assuming data in primary and no headers.")
                try:
                    self.image_noPlanets = hdulist[0].data
                except:
                    raise Exception("Couldn't read "+self.filename_noPlanets_path+". Is it a fits?")

        # if self.filename_noSignal is not None:
        #     # Check file existence and define filename_path
        #     if self.inputDir is None:
        #         try:
        #             self.filename_noSignal_path = os.path.abspath(glob(self.filename_noSignal)[self.id_matching_file])
        #             self.N_matching_files = len(glob(self.filename_noSignal))
        #         except:
        #             raise Exception("File "+self.filename_noSignal+"doesn't exist.")
        #     else:
        #         try:
        #             self.filename_noSignal_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_noSignal)[self.id_matching_file])
        #             self.N_matching_files = len(glob(self.inputDir+os.path.sep+self.filename_noSignal))
        #         except:
        #             raise Exception("File "+self.inputDir+os.path.sep+self.filename_noSignal+" doesn't exist.")
        #
        #     # Open the fits file on which the metric will be applied
        #     hdulist = pyfits.open(self.filename_noSignal_path)
        #     if not self.mute:
        #         print("Opened: "+self.filename_noSignal_path)
        #
        #     # grab the data and headers
        #     try:
        #         self.image_noSignal = hdulist[1].data
        #     except:
        #         # This except was used for datacube not following GPI headers convention.
        #         if not self.mute:
        #             print("Couldn't read the fits file with GPI conventions. Try assuming data in primary.")
        #         try:
        #             self.image_noSignal = hdulist[0].data
        #         except:
        #             raise Exception("Couldn't read "+self.filename_noSignal_path+". Is it a fits?")
        # else:
        #     self.image_noSignal = None
        #     self.filename_noSignal_path = None


        # Get center of the image (star position)
        try:
            # Retrieve the center of the image from the fits headers.
            self.center = [self.exthdr['PSFCENTX'], self.exthdr['PSFCENTY']]
        except:
            # If the keywords could not be found the center is defined as the middle of the image
            if not self.mute:
                print("Couldn't find PSFCENTX and PSFCENTY keywords.")
            self.center = [(self.nx-1)/2,(self.ny-1)/2]

        # if self.label == "CADI":
        #     self.center = [140,140]

        try:
            self.folderName = self.exthdr["METFOLDN"]+os.path.sep
        except:
            try:
                self.folderName = self.exthdr["STAFOLDN"]+os.path.sep
            except:
                pass

        file_ext_ind = os.path.basename(self.filename_path)[::-1].find(".")
        self.prefix = os.path.basename(self.filename_path)[:-(file_ext_ind+1)]
        #self.prefix = "".join(os.path.basename(self.filename_path).split(".")[0:-1])
        self.suffix = self.type
        if self.image_wide is None:
            tmp_suffix = ""
            if self.Dr is not None:
                tmp_suffix = tmp_suffix+"Dr"+str(self.Dr)
            elif self.N is not None:
                tmp_suffix = tmp_suffix+"N"+str(self.N)
            if self.r_step is not None:
                tmp_suffix = tmp_suffix+"rs"+str(self.r_step)
            if self.kernel_type is not None:
                tmp_suffix = tmp_suffix+self.kernel_type
                # if self.kernel_width is not None:
                #     tmp_suffix = tmp_suffix+str(self.kernel_width)
        else:
            tmp_suffix = "IW"
        self.suffix = self.suffix+tmp_suffix



        if self.kernel_type is not None:
            self.ny_PSF = 21 # should be even
            self.nx_PSF = 21 # should be even
            # Define the PSF as a gaussian
            if self.kernel_type == "gaussian":
                if self.kernel_width == None:
                    self.kernel_width = 1.25
                    if not self.mute:
                        print("Default width sigma = {0} used for the gaussian".format(self.kernel_width))

                if not self.mute:
                    print("Generate gaussian PSF")
                # Build the grid for PSF stamp.
                x_PSF_grid, y_PSF_grid = np.meshgrid(np.arange(0,self.ny_PSF,1)-self.ny_PSF/2,
                                                     np.arange(0,self.nx_PSF,1)-self.nx_PSF/2)

                self.PSF = kppmath.gauss2d(x_PSF_grid, y_PSF_grid,1.0,0.0,0.0,self.kernel_width,self.kernel_width)

            # Define the PSF as an aperture or "hat" function
            if self.kernel_type == "hat":
                if self.kernel_width == None:
                    self.kernel_width = 1.5
                    if not self.mute:
                        print("Default radius = {0} used for the hat function".format(self.kernel_width))

                # Build the grid for PSF stamp.
                x_PSF_grid, y_PSF_grid = np.meshgrid(np.arange(0,self.ny_PSF,1)-self.ny_PSF/2,
                                                     np.arange(0,self.nx_PSF,1)-self.nx_PSF/2)
                # Use aperture for the cross correlation.
                # Calculate the corresponding hat function
                self.PSF = kppmath.hat(x_PSF_grid, y_PSF_grid, self.kernel_width)

            self.PSF = self.PSF / np.sqrt(np.nansum(self.PSF**2))




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

        if self.rm_edge is not None:
            # Mask out a band of 10 pixels around the edges of the finite pixels of the image.
            IWA,OWA,inner_mask,outer_mask = get_occ(self.image, centroid = self.center)
            conv_kernel = np.ones((self.rm_edge,self.rm_edge))
            wider_mask = correlate2d(outer_mask,conv_kernel,mode="same")
            self.image[np.where(np.isnan(wider_mask))] = np.nan

        # If GOI_list_folder is not None. Mask the known objects from the image that will be used for calculating the
        # PDF. This masked image is given separately to the probability calculation function.
        if self.filename_noPlanets is not None:
            self.image_without_planet = mask_known_objects(self.image_noPlanets,self.prihdr_noPlanets,self.exthdr_noPlanets,self.GOI_list_folder, mask_radius = self.mask_radius)
        else:
            self.image_without_planet = mask_known_objects(self.image,self.prihdr,self.exthdr,self.GOI_list_folder, mask_radius = self.mask_radius)
        # if self.image_noSignal is None:
        #     if self.GOI_list_folder is not None:
        #         self.image_noSignal = mask_known_objects(self.image,self.prihdr,self.exthdr,self.GOI_list_folder, mask_radius = self.mask_radius)
        #     else:
        #         self.image_noSignal = self.image

        if self.collapse:
            if self.weights is not None:
                image_collapsed = np.zeros((self.ny,self.nx))
                image_without_planet_collapsed = np.zeros((self.ny,self.nx))
                for k in range(self.nl):
                    image_collapsed = image_collapsed + self.weights[k]*self.image[k,:,:]
                    image_without_planet_collapsed = image_without_planet_collapsed + self.weights[k]*self.image_without_planet[k,:,:]
                self.image = image_collapsed/np.sum(self.weights)
                self.image_without_planet = image_without_planet_collapsed/np.sum(self.weights)
            else:
                self.image = np.nanmean(self.image,axis=0)
                self.image_without_planet = np.nanmean(self.image_without_planet,axis=0)

        if self.nans2zero:
            self.image = np.nan_to_num(self.image)

        if self.kernel_type is not None:
            # Check if the input file is 2D or 3D
            if np.size(self.image.shape) == 3: # If the file is a 3D cube
                for l_id in np.arange(self.nl):
                    self.image[l_id,:,:] = correlate2d(self.image[l_id,:,:],self.PSF,mode="same")
                    self.image_without_planet[l_id,:,:] = correlate2d(self.image_without_planet[l_id,:,:],self.PSF,mode="same")
            else: # image is 2D
                self.image = correlate2d(self.image,self.PSF,mode="same")
                self.image_without_planet = correlate2d(self.image_without_planet,self.PSF,mode="same")


        if np.size(self.image.shape) == 3:
            # Not tested
            self.stat_cube_map = np.zeros(self.image.shape)
            for k in range(self.nl):
                self.stat_cube_map[k,:,:] = get_image_stat_map(self.image[k,:,:],
                                                               self.image_without_planet[k,:,:],
                                                               IOWA = self.IOWA,
                                                               N = self.N,
                                                               centroid = self.center,
                                                               r_step = self.r_step,
                                                               mute = self.mute,
                                                               Dr= self.Dr,
                                                               type = self.type,
                                                               image_wide=self.image_wide)

            # if self.collapse:
            #     if self.weights is not None:
            #         stat_collapsed_map = np.zeros((self.ny,self.nx))
            #         for k in range(self.nl):
            #             stat_collapsed_map = stat_collapsed_map + self.weights[k]*self.stat_cube_map[k,:,:]
            #         self.stat_cube_map = stat_collapsed_map/np.sum(self.weights)
            #     else:
            #         self.stat_cube_map = np.nanmean(self.stat_cube_map,axis=0)
        elif np.size(self.image.shape) == 2:
            self.stat_cube_map = get_image_stat_map(self.image,
                                                    self.image_without_planet,
                                                    IOWA = self.IOWA,
                                                    N = self.N,
                                                    centroid = self.center,
                                                    r_step = self.r_step,
                                                    mute = self.mute,
                                                    Dr= self.Dr,
                                                    type = self.type,
                                                    image_wide=self.image_wide)
        return self.stat_cube_map


    def save(self):
        """

        :return: None
        """

        if not os.path.exists(self.outputDir+os.path.sep+self.folderName):
            os.makedirs(self.outputDir+os.path.sep+self.folderName)

        if hasattr(self,"prihdr") and hasattr(self,"exthdr"):
            # Save the parameters as fits keywords
            # STA##### stands for STAtistic
            self.exthdr["STA_TYPE"] = self.type

            self.exthdr["STAFILEN"] = self.filename_path
            self.exthdr["STAINDIR"] = self.inputDir
            self.exthdr["STAOUTDI"] = self.outputDir
            self.exthdr["STAFOLDN"] = self.folderName

            self.exthdr["STAMASKR"] = self.mask_radius
            self.exthdr["STA_IOWA"] = str(self.IOWA)
            self.exthdr["STARSTEP"] = str(self.r_step)
            self.exthdr["STA_N"] = self.N
            self.exthdr["STA_DR"] = self.Dr
            self.exthdr["STA_TYPE"] = self.type
            self.exthdr["STARMEDG"] = self.rm_edge
            self.exthdr["STAGOILF"] = self.GOI_list_folder
            self.exthdr["STAKERTY"] = str(self.kernel_type)
            self.exthdr["STAKERWI"] = str(self.kernel_width)
            self.exthdr["STAIMWID"] = str(self.image_wide)

            # This parameters are not always defined
            if hasattr(self,"filename_noSignal_path"):
                self.exthdr["STAFILNS"] = self.filename_noSignal_path

            if not self.mute:
                print("Saving: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')
            hdulist = pyfits.HDUList()
            hdulist.append(pyfits.PrimaryHDU(header=self.prihdr))
            hdulist.append(pyfits.ImageHDU(header=self.exthdr, data=self.stat_cube_map, name=self.suffix))
            hdulist.writeto(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits', clobber=True)
        else:
            hdulist = pyfits.HDUList()
            hdulist.append(pyfits.ImageHDU(data=self.stat_cube_map, name=self.suffix))
            hdulist.append(pyfits.ImageHDU(name=self.suffix))

            hdulist[1].header["STA_TYPE"] = self.type

            hdulist[1].header["STAFILEN"] = self.filename_path
            hdulist[1].header["STAFILNS"] = self.filename_noSignal_path
            hdulist[1].header["STAINDIR"] = self.inputDir
            hdulist[1].header["STAOUTDI"] = self.outputDir
            hdulist[1].header["STAFOLDN"] = self.folderName

            hdulist[1].header["STAMASKR"] = self.mask_radius
            hdulist[1].header["STA_IOWA"] = self.IOWA
            hdulist[1].header["STARSTEP"] = str(self.r_step)
            hdulist[1].header["STA_N"] = self.N
            hdulist[1].header["STA_DR"] = self.Dr
            hdulist[1].header["STA_TYPE"] = self.type
            hdulist[1].header["STARMEDG"] = self.rm_edge
            hdulist[1].header["STAGOILF"] = self.GOI_list_folder
            hdulist[1].header["STAKERTY"] = str(self.kernel_type)
            hdulist[1].header["STAKERWI"] = str(self.kernel_width)
            hdulist[1].header["STAIMWID"] = str(self.image_wide)

            if not self.mute:
                print("Saving: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits')
            hdulist.writeto(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.fits', clobber=True)

        return None

    def load(self):
        """

        :return: None
        """

        return None