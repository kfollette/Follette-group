__author__ = 'JB'
import os
import astropy.io.fits as pyfits
from glob import glob
import multiprocessing as mp
import numpy as np
from scipy.signal import convolve2d

from pyklip.kpp.utils.kppSuperClass import KPPSuperClass
from pyklip.kpp.stat.stat_utils import *
from pyklip.kpp.utils.GOI import *
import pyklip.kpp.utils.mathfunc as kppmath
import pyklip.kpp.utils.GPIimage as gpiim

class Detection(KPPSuperClass):
    """
    Class for SNR calculation.
    """
    def __init__(self,filename,
                 inputDir = None,
                 outputDir = None,
                 mute=None,
                 N_threads=None,
                 label = None,
                 mask_radius = None,
                 threshold = None,
                 maskout_edge = None,
                 overwrite = False,
                 IWA = None,
                 OWA = None):
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
        super(Detection, self).__init__(filename,
                                     inputDir = inputDir,
                                     outputDir = outputDir,
                                     folderName = None,
                                     mute=mute,
                                     N_threads=N_threads,
                                     label=label,
                                     overwrite = overwrite)

        if mask_radius is None:
            self.mask_radius = 4
        else:
            self.mask_radius = mask_radius

        if threshold is None:
            self.threshold = 3
        else:
            self.threshold = threshold

        # If true mask out a band of 10pix at the edge of the image following the nan boundary.
        if maskout_edge is None:
            self.maskout_edge = False
        else:
            self.maskout_edge = maskout_edge

        self.IWA = IWA
        self.OWA = OWA


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
        init_out = super(Detection, self).initialize(inputDir = inputDir,
                                         outputDir = outputDir,
                                         folderName = folderName,
                                         label=label)


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

        try:
            self.folderName = self.exthdr["METFOLDN"]+os.path.sep
        except:
            try:
                self.folderName = self.exthdr["STAFOLDN"]+os.path.sep
            except:
                pass

        file_ext_ind = os.path.basename(self.filename_path)[::-1].find(".")
        self.prefix = os.path.basename(self.filename_path)[:-(file_ext_ind+1)]
        self.suffix = "DetecTh{0}Mr{1}".format(self.threshold,self.mask_radius)

        return init_out

    def check_existence(self):
        """

        :return: False
        """

        file_exist = (len(glob(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.csv')) >= 1)

        if file_exist and not self.mute:
            print("Output already exist: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.csv')

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


        # Make a copy of the criterion map because it will be modified in the following.
        # Local maxima are indeed masked out when checked
        image_cpy = copy(self.image)

        # Build as grids of x,y coordinates.
        # The center is in the middle of the array and the unit is the pixel.
        # If the size of the array is even 2n x 2n the center coordinate in the array is [n,n].
        x_grid, y_grid = np.meshgrid(np.arange(0,self.nx,1)-self.center[0],np.arange(0,self.ny,1)-self.center[1])


        # Definition of the different masks used in the following.
        stamp_size = self.mask_radius * 2 + 2
        # Mask to remove the spots already checked in criterion_map.
        stamp_x_grid, stamp_y_grid = np.meshgrid(np.arange(0,stamp_size,1)-6,np.arange(0,stamp_size,1)-6)
        stamp_mask = np.ones((stamp_size,stamp_size))
        r_stamp = abs((stamp_x_grid) +(stamp_y_grid)*1j)
        stamp_mask[np.where(r_stamp < self.mask_radius)] = np.nan

        # Mask out a band of 10 pixels around the edges of the finite pixels of the image.
        if self.maskout_edge:
            IWA,OWA,inner_mask,outer_mask = get_occ(self.image, centroid = self.center)
            conv_kernel = np.ones((10,10))
            flat_cube_wider_mask = convolve2d(outer_mask,conv_kernel,mode="same")
            image_cpy[np.where(np.isnan(flat_cube_wider_mask))] = np.nan


        # Number of rows and columns to add around a given pixel in order to extract a stamp.
        row_m = int(np.floor(stamp_size/2.0))    # row_minus
        row_p = int(np.ceil(stamp_size/2.0))     # row_plus
        col_m = int(np.floor(stamp_size/2.0))    # col_minus
        col_p = int(np.ceil(stamp_size/2.0))     # col_plus

        # Table containing the list of the local maxima with their info
        # Description by column:
        # 1/ index of the candidate
        # 2/ Value of the maximum
        # 3/ Position angle in degree from North in [0,360]
        # 4/ Separation in pixel
        # 5/ Separation in arcsec
        # 6/ x position in pixel
        # 7/ y position in pixel
        # 8/ row index
        # 9/ column index
        self.candidate_table = []
        self.table_labels = ["index","value","PA","Sep (pix)","Sep (as)","x","y","row","col"]
        ## START WHILE LOOP.
        # Each iteration looks at one local maximum in the criterion map.
        k = 0
        max_val_criter = np.nanmax(image_cpy)
        while max_val_criter >= self.threshold:# and k <= max_attempts:
            k += 1
            # Find the maximum value in the current criterion map. At each iteration the previous maximum is masked out.
            max_val_criter = np.nanmax(image_cpy)
            # Locate the maximum by retrieving its coordinates
            max_ind = np.where( image_cpy == max_val_criter )
            row_id,col_id = max_ind[0][0],max_ind[1][0]
            x_max_pos, y_max_pos = x_grid[row_id,col_id],y_grid[row_id,col_id]
            sep_pix = np.sqrt(x_max_pos**2+y_max_pos**2)
            sep_arcsec = gpiim.pix2as(sep_pix)
            pa = np.mod(np.rad2deg(np.arctan2(-x_max_pos,y_max_pos)),360)


            # Mask the spot around the maximum we just found.
            image_cpy[(row_id-row_m):(row_id+row_p), (col_id-col_m):(col_id+col_p)] *= stamp_mask

            if self.IWA is not None:
                if sep_pix < self.IWA:
                    continue
            if self.OWA is not None:
                if sep_pix > self.OWA:
                    continue

            # Store the current local maximum information in the table
            self.candidate_table.append([k,max_val_criter,pa,sep_pix,sep_arcsec,x_max_pos,y_max_pos,row_id,col_id])
        ## END WHILE LOOP.

        return self.candidate_table


    def save(self):
        """

        :return: None
        """

        if not os.path.exists(self.outputDir+os.path.sep+self.folderName):
            os.makedirs(self.outputDir+os.path.sep+self.folderName)

        if not self.mute:
            print("Saving: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.csv')
        with open(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+'.csv', 'w+') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=';')
            csvwriter.writerows([self.table_labels])
            csvwriter.writerows(self.candidate_table)
        return None

    def load(self):
        """

        :return: None
        """

        return None
