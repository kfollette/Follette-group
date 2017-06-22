__author__ = 'JB'
import os
import astropy.io.fits as pyfits
from glob import glob
import multiprocessing as mp
import numpy as np
from scipy.signal import convolve2d
import matplotlib.pyplot as plt
import shutil

from pyklip.kpp.utils.kppSuperClass import KPPSuperClass
from pyklip.kpp.utils.GOI import *

class CADIQuicklook(KPPSuperClass):
    """
    Class for CADI quicklook.
    """
    def __init__(self,
                 inputDir = None,
                 outputDir = None,
                 mute=None,
                 label = None,
                 GOI_list_folder = None,
                 overwrite = False,
                 copy_save = None):
        """


        :param filename: Filename of the file on which to calculate the metric. It should be the complete path unless
                        inputDir is defined.
                        It can include wild characters. The file will be selected using the first output of glob.glob().
        :param filename_noSignal: One should be careful with this one since it requires it should find the same number
                            of files with no signal than normal images when calling glob.glob().
                            Besides one has to check that the ordering of glob outputs are matching for both lists.
        :param mute: If True prevent printed log outputs.
        :param N_threads: Number of threads to be used for the metrics and the probability calculations.
                        If None use mp.cpu_count().
                        If -1 do it sequentially.
                        Note that it is not used for this super class.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        """
        # allocate super class
        super(CADIQuicklook, self).__init__("No filename for CADI Quicklook",
                                     inputDir = inputDir,
                                     outputDir = outputDir,
                                     folderName = None,
                                     mute=mute,
                                     N_threads=None,
                                     label=label,
                                     overwrite = overwrite)

        self.copy_save = copy_save
        self.GOI_list_folder = GOI_list_folder
        self.suffix = "quicklook"


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
        init_out = super(CADIQuicklook, self).initialize(inputDir = inputDir,
                                         outputDir = outputDir,
                                         folderName = folderName,
                                         label = label,
                                         read = False)

        # self.filename_nosdi = "cadi_*_nosdi_hp4.fits"
        # self.filename_sdi = "cadi_*_sdi_hp4.fits"
        # self.filename_nosdiSNR = "planet_detec_CADI"+os.path.sep+"default_out"+os.path.sep+"cadi_*_nosdi_hp4-SNR_Dr2rs2.fits"
        # self.filename_sdiSNR = "planet_detec_CADI"+os.path.sep+"default_out"+os.path.sep+"cadi_*_sdi_hp4-SNR_Dr2rs2.fits"

        self.filename_nosdi = "cadi_*_nosdi_hp4.fits"
        self.filename_sdi = "cadi_*_sdi_hp4.fits"
        self.filename_nosdiSNR = "kpop_CADI"+os.path.sep+"default_out"+os.path.sep+"*_nosdi_hp4-SNRDr2rs2hat.fits"
        self.filename_sdiSNR = "kpop_CADI"+os.path.sep+"default_out"+os.path.sep+"*_sdi_hp4-SNRDr2rs2hat.fits"
        # self.filename_nosdiSNR = "planet_detec_CADI"+os.path.sep+"default_out"+os.path.sep+"*_noSDI-shape2Dhat-SNR_Dr2rs2.fits"
        # self.filename_sdiSNR = "planet_detec_CADI"+os.path.sep+"default_out"+os.path.sep+"*_SDI-shape2Dhat-SNR_Dr2rs2.fits"

        # Check file existence and define filename_path
        if self.inputDir is None:
            try:
                self.filename_nosdi_path = os.path.abspath(glob(self.filename_nosdi)[0])
            except:
                raise Exception("File "+self.filename_nosdi+"doesn't exist.")
            try:
                self.filename_sdi_path = os.path.abspath(glob(self.filename_sdi)[0])
            except:
                raise Exception("File "+self.filename_sdi+"doesn't exist.")
            try:
                self.filename_nosdiSNR_path = os.path.abspath(glob(self.filename_nosdiSNR)[0])
            except:
                raise Exception("File "+self.filename_nosdiSNR+"doesn't exist.")
            try:
                self.filename_sdiSNR_path = os.path.abspath(glob(self.filename_sdiSNR)[0])
            except:
                raise Exception("File "+self.filename_sdiSNR+"doesn't exist.")
        else:
            try:
                self.filename_nosdi_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_nosdi)[0])
            except:
                raise Exception("File "+self.inputDir+os.path.sep+self.filename_nosdi+"doesn't exist.")
            try:
                self.filename_sdi_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_sdi)[0])
            except:
                raise Exception("File "+self.inputDir+os.path.sep+self.filename_sdi+"doesn't exist.")
            try:
                self.filename_nosdiSNR_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_nosdiSNR)[0])
            except:
                raise Exception("File "+self.inputDir+os.path.sep+self.filename_nosdiSNR+"doesn't exist.")
            try:
                self.filename_sdiSNR_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_sdiSNR)[0])
            except:
                raise Exception("File "+self.inputDir+os.path.sep+self.filename_sdiSNR+"doesn't exist.")

        # Define this attribute in case something needs it.
        self.filename_path = self.filename

        # Open the fits file on which the metric will be applied
        hdulist1 = pyfits.open(self.filename_nosdi_path)
        if not self.mute:
            print("Opened: "+self.filename_nosdi_path)
        hdulist2 = pyfits.open(self.filename_sdi_path)
        if not self.mute:
            print("Opened: "+self.filename_sdi_path)
        hdulist3 = pyfits.open(self.filename_nosdiSNR_path)
        if not self.mute:
            print("Opened: "+self.filename_nosdiSNR_path)
        hdulist4 = pyfits.open(self.filename_sdiSNR_path)
        if not self.mute:
            print("Opened: "+self.filename_sdiSNR_path)

        # grab the data and headers
        try:
            self.image_nosdi = hdulist1[1].data
            self.exthdr = hdulist1[1].header
            self.prihdr = hdulist1[0].header
            self.ny,self.nx = self.image_nosdi.shape
            self.image_sdi = hdulist2[1].data
            self.image_nosdiSNR = hdulist3[1].data
            self.image_sdiSNR = hdulist4[1].data
        except:
            # This except was used for datacube not following GPI headers convention.
            if not self.mute:
                print("Couldn't read the fits file with GPI conventions. Try assuming data in primary.")
            try:
                self.image_nosdi = hdulist1[0].data
                self.image_sdi = hdulist2[0].data
                self.image_nosdiSNR = hdulist3[0].data
                self.image_sdiSNR = hdulist4[0].data
            except:
                raise Exception("Couldn't read one of the files.")

        hdulist1.close()
        hdulist2.close()
        hdulist3.close()
        hdulist4.close()

        # Get center of the image (star position)
        try:
            # Retrieve the center of the image from the fits headers.
            self.center = [self.exthdr['PSFCENTX'], self.exthdr['PSFCENTY']]
        except:
            # If the keywords could not be found the center is defined as the middle of the image
            if not self.mute:
                print("Couldn't find PSFCENTX and PSFCENTY keywords.")
            self.center = [(self.nx-1)/2,(self.ny-1)/2]


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
            # If the object name could not be found cal lit unknown_object
            self.star_name = "UNKNOWN_OBJECT"

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
        self.N_cubes = len(cube_filename_list)

        self.folderName = ""
        self.prefix = self.star_name+"_"+self.compact_date+"_"+self.filter

        return init_out

    def check_existence_noInit(self,outputDir = None,folderName = None):
        """

        :return: False
        """
        if self.outputDir is not None and self.folderName is not None:
            glob_filename = glob(self.outputDir+os.path.sep+self.folderName+os.path.sep+"*"+'-'+self.suffix+"_*"+'.png')
        else:
            try:
                glob_filename = glob(outputDir+os.path.sep+folderName+os.path.sep+"*"+'-'+self.suffix+"_*"+'.png')
            except:
                raise Exception("self.outputDir or self.folderName are not defined because initialized hasn't been called. You need to give outputDir and folderName as inputs.")

        file_exist = (len(glob_filename) >= 1)

        if self.overwrite and not self.mute:
            print("Overwriting is turned ON!")

        return file_exist and not self.overwrite

    def check_existence(self):
        """

        :return: False
        """
        glob_filename = glob(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_*"+'.png')
        file_exist = (len(glob_filename) >= 1)

        if file_exist:
            current_N_cubes = int(glob_filename[0].split("_")[-1].split(".")[0])
            if not self.mute:
                print("{0}/{1} cubes quicklook output already exist: ".format(current_N_cubes,self.N_cubes)+glob_filename[0])

            if current_N_cubes < self.N_cubes:
                print("More cubes available so overwriting: {0}".format(self.N_cubes))
                file_exist = False

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

        #x_grid, y_grid = np.meshgrid(np.arange(0,self.nx,1)-self.center[0],np.arange(0,self.ny,1)-self.center[1])
        fig = plt.figure(1,figsize=(8*2,16))
        cmap_name = "viridis"
        fig.suptitle(self.star_name+" "+self.compact_date+" Filter: "+self.filter+" Cubes: {0}".format(self.N_cubes), fontsize=25)
        plt.subplot(2,2,1)
        plt.title("NO SDI",y=1.0)
        #np.log10(self.image_nosdi[::-1,:]-np.nanmin(self.image_nosdi[::-1,:]))
        plt.imshow(self.image_nosdi[::-1,:], interpolation="nearest",cmap=cmap_name)#,extent=[x_grid[0,0],x_grid[0,nx-1],y_grid[0,0],y_grid[ny-1,0]])
        plt.clim(-0.01,0.01)
        plt.subplot(2,2,2)
        plt.title("SDI",y=1.0)
        #np.log10(self.image_sdi[::-1,:]-np.nanmin(self.image_sdi[::-1,:]))
        plt.imshow(self.image_sdi[::-1,:], interpolation="nearest",cmap=cmap_name)#,extent=[x_grid[0,0],x_grid[0,nx-1],y_grid[0,0],y_grid[ny-1,0]])
        plt.clim(-0.01,0.01)
        plt.subplot(2,2,3)
        plt.title("SNR NO SDI",y=1.0)
        plt.imshow(self.image_nosdiSNR[::-1,:], interpolation="nearest",cmap=cmap_name)#,extent=[x_grid[0,0],x_grid[0,nx-1],y_grid[0,0],y_grid[ny-1,0]])
        plt.clim(-5.,5.0)
        plt.subplot(2,2,4)
        plt.title("SNR SDI",y=1.0)
        plt.imshow(self.image_sdiSNR[::-1,:], interpolation="nearest",cmap=cmap_name)#,extent=[x_grid[0,0],x_grid[0,nx-1],y_grid[0,0],y_grid[ny-1,0]])
        plt.clim(-5.,5.0)
        #plt.show()
        #plt.title(self.star_name)

        return None


    def save(self):
        """

        :return: None
        """
        # Removing the old quicklook if there is already one available with less input cubes.
        glob_filename = glob(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_*"+'.png')
        file_exist = (len(glob_filename) >= 1)

        if file_exist and not self.mute:
            current_N_cubes = int(glob_filename[0].split("_")[-1].split(".")[0])
            if current_N_cubes < self.N_cubes:
                print("Removing file: "+glob_filename[0])
                os.remove(glob_filename[0])

        if not os.path.exists(self.outputDir+os.path.sep+self.folderName):
            os.makedirs(self.outputDir+os.path.sep+self.folderName)

        plt.figure(1)
        if not self.mute:
            print("Saving: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_"+str(self.N_cubes)+'.png')
        plt.savefig(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_"+str(self.N_cubes)+'.png',
                    bbox_inches='tight')
        if self.copy_save is not None:
            if not self.mute:
                print(self.copy_save+os.path.sep+'CADI_quicklook_current.png')
            src = self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_"+str(self.N_cubes)+'.png'
            dst = self.copy_save+os.path.sep+'CADI_quicklook_current.png'
            shutil.copyfile(src, dst)


        plt.close(1)

        return None

    def load(self):
        """

        :return: None
        """

        return None