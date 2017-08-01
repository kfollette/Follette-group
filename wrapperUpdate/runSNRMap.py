import SNRMap as snr 
import sys


#Program creates signal to noise ratio map using command line inputs via the SNRMap module

#Input should be in format: python runSNRMap.py med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits None True None 

#First argument is the filepath to the fits file to be analyzed, second argument is parameters for a mask (or None), third argument is a boolean to determine whether the file is saved, fourth argument is desired output file name (or None for default file name 



pathToFile = str(sys.argv[1])
print("File Path = " + pathToFile)


if (sys.argv[2] == "None" or sys.argv[2] == "none"):
    planets = None
else:
    planets = (sys.argv[2])
    print("masking pixels for parameters: " + str(planets))

saveOutput = False
if (sys.argv[3] == 'true' or sys.argv[3] == 'True'):
    saveOutput = True

    
if (saveOutput):
    if (sys.argv[4] == "None" or sys.argv[4] == "none"):
        outputName = None
    else:
        outputName = sys.argv[4]
        
else:
    outputName = None
        
snr.create_map(pathToFile, planets, saveOutput, outputName)