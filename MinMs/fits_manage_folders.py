
# coding: utf-8

# In[ ]:

import glob
import os
import file_manager   #requires file_manager.py

toplevel = os.getcwd() # gets the path to the current directory - this should be the top-level directory containing all the HIP... folders

directories = glob.glob('HIP*/data*') # makes a list of all of everything (should only be directories) that start with HIP

n = len(directories) # number of directories

for folder in directories:
    destination = toplevel + '/' + folder # defines the path to your directory of interest
    os.chdir(destination) # moves into that directory

    print('Now in ' + destination + ', making list of FITS.')
    
    HIP_list = glob.glob('N*.fits')
    
    # run the file_manager function on the current FITS list
    file_manager.file_manager(HIP_list)
    
    print('Ran file_manager in the ' + destination + ' directory')

    os.chdir(toplevel) # move back to the top level folder

print(f'Finished {n} directories. Returned to the top level directory.')

