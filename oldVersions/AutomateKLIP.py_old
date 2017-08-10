
# coding: utf-8

# In[ ]:

#Clare Leonard                                                                
#Version 1.0 - 6/21/17                                                          

import glob
import inspect
import os
#import pyklip.instruments.MAGAO as MAGAO                                       
import MAGAO as MAGAO
#import pyklip.parallelized as parallelized                                     
import parallelized as parallelized
import numpy as np
import sys
#import pyklip.klip as klip                                                     
import klip as klip
from astropy.io import fits
#import SNRMap as SNR   


# In[ ]:

pathToFiles = sys.argv[1]
print("File Path = " + pathToFiles)

annuli2_start = int(sys.argv[5])
annuli2_stop = int(sys.argv[6])
annuli2_inc = int(sys.argv[7])
print("Annuli = " + str(annuli2))

iwa = int(sys.argv[2])
print("IWA = " + str(iwa))

movement2_start = float(sys.argv[8])
movement2_stop = int(sys.argv[9])
movement2_inc = int(sys.argv[10])
print("Movement = " + str(movement2))

outputFileName = sys.argv[4]
print("Output FileName = " + outputFileName)

print("KL Modes = " + str(list(map(int, sys.argv[3].split(",")))))
klmodes = list(map(int, sys.argv[3].split(",")))

subsections2_start = int(sys.argv[11])
subsections2_stop = int(sys.argv[12])
subsections2_inc = int(sys.argv[13])
print("Subsections = " + str(subsections2))

