import numpy as np
import scipy as sp
import scipy.ndimage
from astropy.io import fits
import sys
import numpy.ma as ma
import math
import os
import glob
import fnmatch
import shutil

def display(keyword, header_keyword):
    '''
    File_Display_Sorter.display:  
    Usage: Prints out the specific header keyword information for an entire array of images: 
    Arguments: keyword = string - file name (ex. 'im_9999.fits' or 'im*')
               header_keyword = string - Specific header subtitle name (ex. 'FILTER' or 'EXPTIME')
    History:
    (v1) 2019-03-28 - SGC: Current version
    '''
    
    image_list = glob.glob(keyword)  #makes an array of fits files according to the keyword given in string form
    for files in image_list:
        header = fits.getheader(files, ignore_missing_end=True)  #Gets header of each fits file
    
        print(header[header_keyword])                   #Print each result for the given header subtitle

def sorter(keyword1, keyword2, header_keyword):
    '''
    File_Display_Sorter.sorter:  
    Usage: Makes a new subdirectory and moves all the corresponding files with a specific keyword into that subdirectory: 
    Arguments: keyword1 = string - file name (ex. 'im_9999.fits' or 'im*')
               keyword2 = string - name for subdirectory, should match the header subtitle iinput (ex. If OBJETCT = 'dark' then keyword2 is 'dark') 
               header_keyword = string - Specific header subtitle name (ex. 'FILTER' or 'EXPTIME')
    History:
    (v1) 2019-03-28 - SGC: Current version
    '''

    image_list = glob.glob(keyword1) #makes an array of fits images according to the keyword given in string form
    
    if not os.path.exists(keyword2):   #makes a subdirectory using keyword2 if the subdirectory does not already exists. 
        os.makedirs(keyword2)
    
    for files in image_list:
        
        header = fits.getheader(files, ignore_missing_end=True)  #gets header of each image
    
        if header[header_keyword] == keyword2:      #if the subtitle header input is equal to keyword2 then print the subtitle 
            print(header[header_keyword])           #header input and move the mathcing file to the corresponding subdirectory
            os.rename(files, keyword2 + '/' + files)

    keyword2_nospace = keyword2.split(' ')[0] + keyword2.split(' ')[1]
    print('Renaming directory to', keyword2_nospace)
    os.rename(keyword2, keyword2_nospace)
            
def exptime(keyword, time, time_file):
    '''
    File_Display_Sorter.sorter:  
    Usage: Makes a new subdirectory for time exposures and moves all the corresponding files with a specific keyword into that subdirectory: 
    Arguments: keyword = string - file name (ex. 'im_9999.fits' or 'im*')
               time = number - exposure time
               time_file = string - exxposure time file name 
    History:
    (v1) 2019-03-28 - SGC: Current version
    '''
    image_list = glob.glob(keyword)     #makes array wuth given keyword
    
    if not os.path.exists(time_file):   #makes a directory for exposure time if not yet existed
        os.makedirs(time_file)
    
    for files in image_list:
        header = fits.getheader(files, ignore_missing_end=True)  #gets header
    
        if header['EXPTIME'] == time:    #Moves the corresponding fies to subdirectory
            print(header['EXPTIME'])
            os.rename(files, time_file + '/' + files)

