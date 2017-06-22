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

class Quicklook(KPPSuperClass):
    """
    Class for CADI quicklook.
    """
    def __init__(self,filename_proba,filename_detec,
                 inputDir = None,
                 outputDir = None,
                 mute=None,
                 label = None,
                 GOI_list_folder = None,
                 overwrite = False,
                 copy_save = None,
                 SNR = None):
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
        super(Quicklook, self).__init__("No filename for Quicklook",
                                     inputDir = inputDir,
                                     outputDir = outputDir,
                                     folderName = None,
                                     mute=mute,
                                     N_threads=None,
                                     label=label,
                                     overwrite = overwrite)

        self.filename_proba = filename_proba
        self.filename_detec = filename_detec

        self.copy_save = copy_save
        self.GOI_list_folder = GOI_list_folder
        self.suffix = label+"quicklook"
        if SNR is None:
            self.SNR = False
        else:
            self.SNR = SNR


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
        init_out = super(Quicklook, self).initialize(inputDir = inputDir,
                                         outputDir = outputDir,
                                         folderName = folderName,
                                         label = label,
                                         read = False)

        # self.filename_nosdi = "cadi_*_nosdi_hp4.fits"
        # self.filename_sdi = "cadi_*_sdi_hp4.fits"
        # self.filename_nosdiSNR = "planet_detec_CADI"+os.path.sep+"default_out"+os.path.sep+"cadi_*_nosdi_hp4-SNR_Dr2rs2.fits"
        # self.filename_sdiSNR = "planet_detec_CADI"+os.path.sep+"default_out"+os.path.sep+"cadi_*_sdi_hp4-SNR_Dr2rs2.fits"


        # Check file existence and define filename_path
        if self.inputDir is None:
            try:
                self.filename_proba_path = os.path.abspath(glob(self.filename_proba)[0])
            except:
                raise Exception("File "+self.filename_proba+"doesn't exist.")
        else:
            try:
                self.filename_proba_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_proba)[0])
            except:
                raise Exception("File "+self.inputDir+os.path.sep+self.filename_proba+"doesn't exist.")

        # Define this attribute in case something needs it.
        self.filename_path = self.filename

        # Open the fits file on which the metric will be applied
        hdulist1 = pyfits.open(self.filename_proba_path)
        if not self.mute:
            print("Opened: "+self.filename_proba_path)

        # grab the data and headers
        try:
            self.image = hdulist1[1].data
            self.exthdr = hdulist1[1].header
            self.prihdr = hdulist1[0].header
            self.ny,self.nx = self.image.shape
        except:
            # This except was used for datacube not following GPI headers convention.
            if not self.mute:
                print("Couldn't read the fits file with GPI conventions. Try assuming data in primary.")
            try:
                self.image = hdulist1[0].data
            except:
                raise Exception("Couldn't read one of the files.")

        hdulist1.close()

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


        # Check file existence and define filename_path
        if self.inputDir is None:
            try:
                self.filename_detec = os.path.abspath(glob(self.filename_detec)[self.id_matching_file])
                self.N_matching_files = len(glob(self.filename_detec))
            except:
                raise Exception("File "+self.filename_detec+"doesn't exist.")
        else:
            try:
                self.filename_detec_path = os.path.abspath(glob(self.inputDir+os.path.sep+self.filename_detec)[self.id_matching_file])
                self.N_matching_files = len(glob(self.inputDir+os.path.sep+self.filename_detec))
            except:
                raise Exception("File "+self.inputDir+os.path.sep+self.filename_detec+" doesn't exist.")

        # Open the fits file on which the metric will be applied
        with open(self.filename_detec_path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=';')
            csv_as_list = list(reader)
            self.detec_table_labels = csv_as_list[0]
            self.detec_table = np.array(csv_as_list[1::], dtype='string').astype(np.float)
        if not self.mute:
            print("Opened: "+self.filename_detec_path)

        self.N_detec = self.detec_table.shape[0]
        self.val_id = self.detec_table_labels.index("value")
        self.x_id = self.detec_table_labels.index("x")
        self.y_id = self.detec_table_labels.index("y")
        self.row_id = self.detec_table_labels.index("row")
        self.col_id = self.detec_table_labels.index("col")
        self.pa_id = self.detec_table_labels.index("PA")
        self.sep_pix_id = self.detec_table_labels.index("Sep (pix)")
        self.sep_as_id = self.detec_table_labels.index("Sep (as)")

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

        if file_exist and not self.mute:
            print("Quicklook output already exist: "+glob_filename[0])


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

        x_grid, y_grid = np.meshgrid(np.arange(0,self.nx,1)-self.center[0],np.arange(0,self.ny,1)-self.center[1])

        fig = plt.figure(1,figsize=(16,8))
        fig.patch.set_facecolor('black')
        cmap_name = "viridis"
        title_obj = fig.suptitle(self.star_name+" "+self.compact_date+" Filter: "+self.filter+" Cubes: {0}".format(self.N_cubes), fontsize=25)
        plt.setp(title_obj, color='white')
        ax0 = plt.subplot2grid((2,4),(0,1),colspan=2,rowspan=2)
        #plt.title("NO SDI",y=1.0)
        #np.log10(self.image_nosdi[::-1,:]-np.nanmin(self.image_nosdi[::-1,:]))
        plt.imshow(self.image[::-1,:], interpolation="nearest",cmap=cmap_name,extent=[x_grid[0,0],x_grid[0,self.nx-1],y_grid[0,0],y_grid[self.ny-1,0]])
        plt.clim(0,5)
        # Remove box and axes ticks
        ax0.set_axis_off()

        if self.GOI_list_folder is not None:
            sep_GOI_vec,pa_GOI_vec = get_pos_known_objects(self.prihdr,self.exthdr,self.GOI_list_folder,pa_sep = True)
            x_GOI_vec,y_GOI_vec = get_pos_known_objects(self.prihdr,self.exthdr,self.GOI_list_folder,xy = True)
            row_centroid_vec,col_centroid_vec = get_pos_known_objects(self.prihdr,self.exthdr,self.GOI_list_folder,xy = False)
            for x_GOI,y_GOI,sep_GOI,pa_GOI in zip(x_GOI_vec,y_GOI_vec,sep_GOI_vec,pa_GOI_vec):
                circle=plt.Circle((x_GOI,y_GOI),radius=10.,color='r',fill=False)
                ax0.add_artist(circle)
                #print(x_GOI,y_GOI,sep_GOI,pa_GOI)
                ax0.annotate("{0:.2f}ac,{1:.1f}deg".format(sep_GOI,pa_GOI), fontsize=10, color = "red",xy=(float(x_GOI), float(y_GOI)),
                            xycoords='data', xytext=(float(x_GOI)+10, float(y_GOI)-10)
                            )
                # ax0.annotate("{0:.1f}ac,{1:.1f}deg".format(sep_GOI,pa_GOI), fontsize=30, color = "red",xy=(float(x_GOI), float(y_GOI)),
                #             xycoords='data', xytext=(float(x_GOI)+10, float(y_GOI)-10),
                #             arrowprops=dict(arrowstyle="->",
                #                             linewidth = 2.,
                #                             color = 'red')
                #             )

        n_stamp = 31
        row_m = np.floor(n_stamp/2.0)    # row_minus
        row_p = np.ceil(n_stamp/2.0)     # row_plus
        col_m = np.floor(n_stamp/2.0)    # col_minus
        col_p = np.ceil(n_stamp/2.0)     # col_plus

        for k in range(2):
            proba = self.detec_table[k,self.val_id]
            col_pos = self.detec_table[k,self.col_id]
            row_pos = self.detec_table[k,self.row_id]
            curr_x_pos = self.detec_table[k,self.x_id]
            curr_y_pos = self.detec_table[k,self.y_id]
            curr_pa = self.detec_table[k,self.pa_id]
            curr_sep_pix = self.detec_table[k,self.sep_pix_id]
            curr_sep_ac = self.detec_table[k,self.sep_as_id]

            if proba < 4:
                cand_color = "white"
            elif 4 <= proba < 6:
                cand_color = "yellow"
            elif 6 <= proba:
                cand_color = "orange"

            circle=plt.Circle((curr_x_pos,curr_y_pos),radius=7.,color=cand_color,fill=False)
            ax0.add_artist(circle)
            #print(x_GOI,y_GOI,sep_GOI,pa_GOI)
            ax0.annotate("{0}: {1:.2f}ac,{2:.1f}deg".format(k,curr_sep_ac,curr_pa), fontsize=10, color = cand_color,xy=(float(curr_x_pos)+3, float(curr_y_pos)+3),
                        xycoords='data', xytext=(float(curr_x_pos)+10, float(curr_y_pos)+15),
                        arrowprops=dict(arrowstyle="->",
                                        linewidth = 2.,
                                        color = cand_color))

            stamp = self.image[(row_pos-row_m):(row_pos+row_p), (col_pos-col_m):(col_pos+col_p)]

            #print((np.mod(k,2),-int(np.floor(k/2.))))
            if k == 0:
                ax = plt.subplot2grid((2,4),(0,0),colspan=1,rowspan=1)
            elif k == 1:
                ax = plt.subplot2grid((2,4),(0,3),colspan=1,rowspan=1)
            title_obj = plt.title("Candidate {0}".format(k))
            plt.setp(title_obj, color='white')
            if self.SNR == False:
                ax.text(0.,1.3*n_stamp,"False Pos. rate: 10^-{0:.1f}\nSep: {1:.1f}pix, {2:.2f}ac\nPA: {3:.1f}deg".format(proba,curr_sep_pix,curr_sep_ac,curr_pa) ,color=cand_color, fontsize=15)
            else:
                ax.text(0.,1.3*n_stamp,"SNR: {0:.1f}\nSep: {1:.1f}pix, {2:.2f}ac\nPA: {3:.1f}deg".format(proba,curr_sep_pix,curr_sep_ac,curr_pa) ,color=cand_color, fontsize=15)

            plt.imshow(stamp[::-1,:], interpolation="nearest",cmap=cmap_name)#,extent=[x_grid[0,0],x_grid[0,nx-1],y_grid[0,0],y_grid[ny-1,0]])
            # Remove box and axes ticks
            ax.set_axis_off()







        #plt.show()
        #plt.title(self.star_name)

        return None


    def save(self):
        """

        :return: None
        """

        if not os.path.exists(self.outputDir+os.path.sep+self.folderName):
            os.makedirs(self.outputDir+os.path.sep+self.folderName)

        fig = plt.figure(1)
        if not self.mute:
            print("Saving: "+self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_"+str(self.N_cubes)+'.png')
        plt.savefig(self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_"+str(self.N_cubes)+'.png',
                    bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
        if self.copy_save is not None:
            if not self.mute:
                print(self.copy_save+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_"+str(self.N_cubes)+'.png')
            src = self.outputDir+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_"+str(self.N_cubes)+'.png'
            dst = self.copy_save+os.path.sep+self.folderName+os.path.sep+self.prefix+'-'+self.suffix+"_"+str(self.N_cubes)+'.png'
            shutil.copyfile(src, dst)


        plt.close(1)


        return None

    def load(self):
        """

        :return: None
        """

        return None
