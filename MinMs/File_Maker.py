
# coding: utf-8

# In[44]:

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
import itertools
from colorama import Fore, Back, Style

fitsfile = glob.glob('*.fits')
root_source = ['/Users/astrolab/Desktop/Data_backup/']  #The main path for this script, can be changed depending on your own root path

dark_root = 'data_with_raw_calibs/DARK'                 #These are the 4 new paths that the script 
object_root = 'data_with_raw_calibs/OBJECT'             #will make stemming off from the root_source 
sky_root = 'data_with_raw_calibs/SKY'
flat_root = 'data_with_raw_calibs/FLAT'


root = ['data_with_raw_calibs/DARK','data_with_raw_calibs/OBJECT',   #making an array of the new paths
        'data_with_raw_calibs/SKY','data_with_raw_calibs/FLAT']

subject = ['HIP103800','HIP103910','HIP104432','HIP104644',        #arrays of all the target subjects, 
           'HIP106811','HIP108569','HIP108706','HIP109638',        #can be changed for different target names
           'HIP111313','HIP112774','HIP116317','HIP12097',
           'HIP13218','HIP19394','HIP22762','HIP24284',
           'HIP25953','HIP28035','HIP29316','HIP31862','HIP36338',
           'HIP36985','HIP38082','HIP42762','HIP44722','HIP4569',
           'HIP47103','HIP47741','HIP48336','HIP50341','HIP52190',
           'HIP52596','HIP55042','HIP56157','HIP56244','HIP61706',
           'HIP65714','HIP70475','HIP70865','HIP72509','HIP72511',
           'HIP74190','HIP78353','HIP79431','HIP84051','HIP84521',
           'HIP84794', 'HIP91430', 'HIP92871', 'HIP93069', 
           'HIP93101', 'HIP93206']

#def fits_file_manager(fitsfile):
    




# In[40]:

part1 = []                                                     #Makes an array for each HIP target path
for i in subject:
    name_folder = '/Users/astrolab/Desktop/Data_backup/' + i
    part1.append(str(name_folder))
    


# In[41]:

DARK = []                                                      #Makes an array for all dark files                         
for i in part1: 
    folder = i + '/data_with_raw_calibs/DARK'
    DARK.append(str(folder))
#############################################
OBJECT = []                                                    #Makes an array for all object files
for i in part1:
    folder = i + '/data_with_raw_calibs/OBJECT'
    OBJECT.append(str(folder))
#############################################
SKY = []                                                       #Makes an array for all sky files
for i in part1:
    folder = i + '/data_with_raw_calibs/SKY'
    SKY.append(str(folder))
#############################################
FLAT_LAMP = []                                                 #Makes an array for all flat,lamp files
for i in part1:
    folder = i + '/data_with_raw_calibs/FLAT_LAMP'
    FLAT_LAMP.append(str(folder))


# In[50]:

for i in DARK:                                                     #Makes a path for each dark file if it doesn't exist already
    os.path.exists(i)
    if not os.path.exists(i):
        os.makedirs(i)
        
for i in OBJECT:                                                   #Makes a path for each object file if it doesn't exist already
    os.path.exists(i)
    if not os.path.exists(i):
        os.makedirs(i)
        
for i in SKY:                                                      #Makes a path for each sky file if it doesn't exist already
    os.path.exists(i)
    if not os.path.exists(i):
        os.makedirs(i)
        
for i in FLAT_LAMP:                                                #Makes a path for each flat,lamp file if it doesn't exist already
    os.path.exists(i)
    if not os.path.exists(i):
        os.makedirs(i)


