import SNRMap as snr 
import sys


#Program creates signal to noise ratio map using command line inputs via the SNRMap module

#Input should be in format: python runSNRMap.py med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits None True None 

#First argument is the filepath to the fits file to be analyzed, second argument is parameters for a mask (or None), third argument is a boolean to determine whether the file is saved, fourth argument is desired output file name (or None for default file name 



pathToFile = str(sys.argv[1])
print("File Path = " + pathToFile)

args = 0
if (sys.argv[2] == "None" or sys.argv[2] == "none"):
    planets = None
else:
    ra = list(map(int, sys.argv[2].split(",")))
    pa =list(map(int, sys.argv[3].split(",")))
    wid = list(map(int, sys.argv[4].split(",")))
    planets = (ra,pa,wid)
    print("masking pixels for parameters: " + str(planets))
    args = 2

saveOutput = False
if (sys.argv[3+args] == 'true' or sys.argv[3+args] == 'True'):
    saveOutput = True

    
if (saveOutput):
    if (sys.argv[4+args] == "None" or sys.argv[4+args] == "none"):
        outputName = None
    else:
        outputName = sys.argv[4+args]
        
else:
    outputName = None
        
snr.create_map(pathToFile, 5,planets, saveOutput, outputName)
