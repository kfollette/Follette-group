__author__ = 'JB'

import os
import astropy.io.fits as pyfits
from glob import glob
import multiprocessing as mp
import numpy as np

class KPPSuperClass(object):
    """
    Super class for all kpop classes (ie FMMF, matched filter, shape, weighted collapse...).
    Has fall-back functions for all metric dependent calls so that each metric class does not need to implement
    functions it doesn't want to.
    It is not a completely empty function and includes features that are probably useful to most inherited class though
    one might decide to overwrite them.
    Here it simply returns the input fits file as read.

    I should remove the option to set output dir in the class definition
    """
    def __init__(self,filename,
                 inputDir = None,
                 outputDir = None,
                 folderName = None,
                 mute=None,
                 N_threads=None,
                 label = None,
                 overwrite = False):
        """
        Define the general parameters of the metric.

        The idea is that the directories and eventual spectra will be defined later by calling the
        initialize() function. Furthermore initialize() is where files can be read.
        It is done that way to ease the analysis of an entire campaign. ie one object is defined
        for a given metric and then a function calls initialize several time for each target's directory and each
        spectra.

        Example for inherited classes:
        It can be the width of the hat function used for the cross correlation in a matched filter for example.

        :param filename: Filename of the file on which to calculate the metric. It should be the complete path unless
                        inputDir is defined.
                        It can include wild characters. The file will be selected using the first output of glob.glob().
        :param inputDir: If defined it allows filename to not include the whole path and just the filename.
                        Files will be read from inputDir.
                        Note tat inputDir might be redefined using initialize at any point.
                        If inputDir is None then filename is assumed to have the absolute path.
        :param outputDir: Directory where to create the folder containing the outputs.
                        Note tat inputDir might be redefined using initialize at any point.
                        If outputDir is None:
                            If inputDir is defined: outputDir = inputDir+os.path.sep+"planet_detec_"+label
                            If inputDir is None: outputDir = "."+os.path.sep+"planet_detec_"+label
        :param folderName: Name of the folder containing the outputs. It will be located in outputDir.
                        Default folder name is "default_out".
                        The convention is to have one folder per spectral template.
                        Usually this folderName should be defined by the class itself and not by the user.
        :param mute: If True prevent printed log outputs.
        :param N_threads: Number of threads to be used for the metrics and the probability calculations.
                        If None use mp.cpu_count().
                        If -1 do it sequentially.
                        Note that it is not used for this super class.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        """

        self.overwrite = overwrite

        # Define a default folderName is the one given is None.
        if folderName is None:
            self.folderName = os.path.sep+"default_out" +os.path.sep
        else:
            self.folderName = folderName+os.path.sep

        self.filename = filename
        # Define the actual filename path
        if inputDir is None: # If None inputDir will be defined in initalize()
            self.inputDir = inputDir
        else:
            self.inputDir = os.path.abspath(inputDir)

        if label is None:
            self.label = "default"
        else:
            self.label = label

        if outputDir is None: # If None outputDir will be defined in initalize()
            self.outputDir = None
        else:
            print("DON'T SET OUTPUTDIR WHEN DEFINING THE CLASS. JB SHOULD REMOVE THIS FEATURE")
            self.outputDir = os.path.abspath(outputDir+os.path.sep+"kpop_"+self.label)

        # Number of threads to be used in case of parallelization.
        if N_threads is None:
            self.N_threads = mp.cpu_count()
        else:
            self.N_threads = N_threads

        self.mute = mute

    def spectrum_iter_available(self):
        """
        Should indicate wether or not the class is equipped for iterating over spectra.
        If the metric requires a spectrum one might one to iterate over several spectra without rereading the input
        files. In order to iterate over spectra the function init_new_spectrum() should be defined.
        spectrum_iter_available is a utility function for campaign data processing to know wether or not spectra the
        metric class should be iterated over different spectra.

        In the case of this super class an therefore by default it returns False.

        :return: False
        """

        return False

    def init_new_spectrum(self,spectrum):
        """
        Function allowing the reinitialization of the class with a new spectrum without reinitializing everything.

        :param spectrum: spectrum path relative to pykliproot + os.path.sep + "spectra" with pykliproot the directory in
                        which pyklip is installed. It that case it should be a spectrum from Mark Marley.
                        Instead of a path it can be a simple ndarray with the right dimension.
                        Or by default it is a completely flat spectrum.

        :return: None
        """

        return None


    def initialize(self,inputDir = None,
                         outputDir = None,
                         folderName = None,
                         label = None,
                         read = True):
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
                        Usually this folderName should be defined by the class itself and not by the user.
        :param label: Define the suffix to the output folder when it is not defined. cf outputDir. Default is "default".
        :param read: If true (default) read the fits file according to inputDir and filename.

        :return: None
        """

        if not hasattr(self,"id_matching_file"):
            self.id_matching_file = 0
        if not hasattr(self,"process_all_files"):
            self.process_all_files = True

        if label is not None:
            self.label = label

        # Define a default folderName is the one given is None.
        if folderName is not None:
            self.folderName = folderName+os.path.sep

        # Define the actual filename path
        if inputDir is None:
            self.inputDir = self.inputDir
        else:
            self.inputDir = os.path.abspath(inputDir)


        if read:
            # Check file existence and define filename_path
            if self.inputDir is None or os.path.isabs(self.filename):
                try:
                    self.filename_path = os.path.abspath(glob(self.filename)[self.id_matching_file])
                    self.N_matching_files = len(glob(self.filename))
                except:
                    raise Exception("File "+self.filename+"doesn't exist.")
            else:
                try:
                    self.filename_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename)[self.id_matching_file])
                    self.N_matching_files = len(glob(self.inputDir+os.path.sep+self.filename))
                except:
                    raise Exception("File "+self.inputDir+os.path.sep+self.filename+" doesn't exist.")

            self.id_matching_file = self.id_matching_file+1

            # Open the fits file on which the metric will be applied
            hdulist = pyfits.open(self.filename_path)
            if not self.mute:
                print("Opened: "+self.filename_path)

            # grab the data and headers
            try:
                self.image = hdulist[1].data
                self.exthdr = hdulist[1].header
                self.prihdr = hdulist[0].header
            except:
                # This except was used for datacube not following GPI headers convention.
                if not self.mute:
                    print("Couldn't read the fits file with GPI conventions. Try assuming data in primary and no headers.")
                try:
                    self.image = hdulist[0].data
                except:
                    raise Exception("Couldn't read "+self.filename_path+". Is it a fits?")

            # Get input cube dimensions
            self.image = np.squeeze(self.image)
            if len(self.image.shape) == 3:
                self.nl,self.ny,self.nx = self.image.shape
                # # Checking that the cube has the 37 spectral slices of a normal GPI cube.
                # if self.nl != 37:
                #     raise Exception("Returning None. Spectral dimension of "+self.filename_path+" is not correct...")
            elif len(self.image.shape) == 2:
                self.ny,self.nx = self.image.shape
            else:
                raise Exception("Returning None. fits file "+self.filename_path+" was not a 2D image or a 3D cube...")

            # If outputDir is None define it as the project directory.
            if outputDir is not None:
                self.outputDir = os.path.abspath(outputDir+os.path.sep+"kpop_"+self.label)
            else: # if self.outputDir is None:
                split_path = np.array(os.path.dirname(self.filename_path).split(os.path.sep))
                if np.sum(word.startswith("planet_detec_") for word in split_path):
                    planet_detec_label = split_path[np.where(["planet_detec" in mystr for mystr in split_path])][-1]
                    self.outputDir = os.path.abspath(self.filename_path.split(planet_detec_label)[0]+planet_detec_label)
                elif np.sum(word.startswith("kpop_") for word in split_path):
                    split_path = np.array(os.path.dirname(self.filename_path).split(os.path.sep))
                    planet_detec_label = split_path[np.where(["kpop" in mystr for mystr in split_path])][-1]
                    self.outputDir = os.path.abspath(self.filename_path.split(planet_detec_label)[0]+planet_detec_label)
                else:
                    self.outputDir = os.path.join(os.path.dirname(self.filename_path),"kpop_"+self.label)

            if self.process_all_files:
                if self.id_matching_file < self.N_matching_files:
                    return True
                else:
                    self.id_matching_file = 0
                    return False
            else:
                return False
        else:
            # If outputDir is None define it as the project directory.
            if outputDir is not None:
                self.outputDir = os.path.abspath(outputDir+os.path.sep+"kpop_"+self.label)

            if self.outputDir is None:
                if self.inputDir is None:
                    self.outputDir = os.path.abspath("."+os.path.sep+"kpop_"+self.label)
                else:
                    self.outputDir = os.path.abspath(self.inputDir+os.path.sep+"kpop_"+self.label)

            return False



    def check_existence(self):
        """
        Check if this metric has already been calculated for this file.

        For this super class it returns False.

        Inherited classes:
        It could check at the output folder if the file with the right extension already exist.

        :return: False
        """

        return False


    def calculate(self):
        """
        Calculate the metric map.

        For this super class it returns the input fits file read in initialize().

        Inherited classes:
        It could check at the output folder if the file with the right extension already exist.

        :return: self.image the imput fits file.
        """

        return self.image


    def save(self):
        """
        Save the metric map as a fits file in self.outputDir+os.path.sep+self.folderName

        For this super class it doesn't do anything.

        Inherited classes:
        It should probably include new fits keywords with the metric parameters before saving the outputs.

        :return: None
        """

        return None


    def load(self):
        """
        Load the metric map if it already exist from self.outputDir+os.path.sep+self.folderName

        For this super class it doesn't do anything.

        :return: None
        """

        return None