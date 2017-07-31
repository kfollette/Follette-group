
#Clare Leonard                                                                
#Version 1.0 - 6/21/17                                                          




##################################################################
#############                                        #############
#############                 IMPORTS                #############
#############                                        #############
##################################################################

import glob
import inspect
import os                                      
import instruments.MagAO as MagAO                                   
import parallelized as parallelized
import numpy as np
import sys                                                   
import klip as klip
from astropy.io import fits
import SNRMap as snr   
import time
import warnings
from astropy.utils.exceptions import AstropyWarning

warnings.filterwarnings('ignore', category=AstropyWarning, append=True)

##################################################################
#############                                        #############
#############               SAVE FILES               #############
#############                                        #############
################################################################## 
 
def writeData(indiv, filepath, filename, annuli, movement, subsections, iwa, klmodes, mask = None, pre = '', suff = ""):    
    #function writes out fits files and writes important information to fits headers
    
    hdu = fits.PrimaryHDU(indiv)
    hdulist = fits.HDUList([hdu])
    prihdr = hdulist[0].header
    
    #creates new strings to add parameter information to file names
    annuli2 = annuli
    movement2 = movement
    subsections2 = subsections
    
    #if program iterates over several parameter values, formats these for fits headers and file names
    if (isinstance(annuli, tuple)):
        annuli = str(annuli[0]) + 'to' + str(annuli[1]) + 'by' + str(annuli[2])
        annuli2 = "a" + str(annuli[0]) + '-' + str(annuli[1]) + 'x' + str(annuli[2])
    if (isinstance(movement, tuple)):
        movement = str(movement[0]) + 'to' + str(movement[1]) + 'by' + str(movement[2])
        movement2 = "m" + str(movement[0]) + '-' + str(movement[1]) + 'x' + str(movement[2])
    if (isinstance(subsections, tuple)):
        subsections = str(subsections[0]) + 'to' + str(subsections[1]) + 'by' + str(subsections[2])
        subsections2 = "s" + str(subsections[0]) + '-' + str(subsections[1]) + '-' + str(subsections[2])
    
    #shortens file path to bottom 4 directories so it will fit in fits header
    pathToFiles_short  = ''
    numdir = 0
    for n in range (len(filepath)):
        pathToFiles_short = filepath[-1-n] + pathToFiles_short
        if (filepath[-1-n] == '/'):
            numdir += 1
        if (numdir >=4 ):
            break
            
    #adds info to fits headers
    prihdr.set('annuli', str(annuli))
    prihdr.set('movement', str(movement))
    prihdr.set('subsctns', str(subsections))
    prihdr.set('klmodes', str(klmodes))
    prihdr.set('filepath', str(pathToFiles_short))
    if (not mask == None):
        rad, pa, wid = mask 
        prihdr.set('mask_rad', str(rad))
        prihdr.set('mask_pa', str(pa))
        prihdr.set('mask_wid', str(wid))
    
    #writes out files
    hdulist.writeto(str(filepath) + "/../" + str(pre) + '_' + filename + "_a" + str(annuli2) + "m" + str(int(movement2)) + "s" + str(subsections2) + "iwa" + str(iwa) + '_' + str(suff) + '_KLmodes-all.fits', clobber=True)




##################################################################
#############                                        #############
#############               GET INPUTS               #############
#############                                        #############
################################################################## 

#value adjusts argument numbering in case of white space in file path 
argnum = 0


pathToFiles = sys.argv[1]

#if filepath doesnt end in sliced, sontinues to add next arguements, helpful in case of whitespace in file path
while (not pathToFiles[-1] == 'd' and not pathToFiles[-1] == '"'):
    argnum += 1
    pathToFiles = pathToFiles + " " + sys.argv[1+argnum]
    
print("File Path = " + pathToFiles)   
print()
 
print("Parameters to explore:")

annuli2_start = int(sys.argv[5+argnum])
annuli2_stop = int(sys.argv[6+argnum])
annuli2_inc = int(sys.argv[7+argnum])
#creates touple for easier eventual string formatting when saving files
annuli = (annuli2_start, annuli2_stop, annuli2_inc)
#if only one parameter is iterated over, makes sure increment is 1 and changes touple to single int
if(annuli2_start == annuli2_stop):
    annuli2_inc = 1;
    annuli = annuli2_start
    print("Annuli = " +str(annuli2_start))    
else:
    print("Annuli: start = %s; end = %s; increment = %s " %(str(annuli2_start), str(annuli2_stop), str(annuli2_inc)))



movement2_start = float(sys.argv[8+argnum])
movement2_stop = float(sys.argv[9+argnum])
movement2_inc = float(sys.argv[10+argnum])
#creates touple for easier eventual string formatting when saving files
movement = (movement2_start, movement2_stop, movement2_inc)
#if parameter is not set to change, makes sure increment is 1 and changes touple to single int
if(movement2_start == movement2_stop):
    movement2_inc = 1;
    movement = movement2_start
    print("Movement = " +str(movement2_start))
    
else:
    print("Movement: start = %s; end = %s; increment = %s " %(str(movement2_start), str(movement2_stop), str(movement2_inc)))


subsections2_start = int(sys.argv[11+argnum])
subsections2_stop = int(sys.argv[12+argnum])
subsections2_inc = int(sys.argv[13+argnum])
#creates touple for easier eventual string formatting when saving files
subsections = (subsections2_start, subsections2_stop, subsections2_inc)
#if parameter is not set to change, makes sure increment is 1 and changes touple to single int
if(subsections2_start == subsections2_stop):
    subsections2_inc = 1;
    subsections = subsections2_start
    print("Subsections = " +str(subsections2_start))    
else:
    print("Subsections: start = %s; end = %s; increment = %s " %(str(subsections2_start), str(subsections2_stop), str(subsections2_inc)))

print()
    
iwa = int(sys.argv[2+argnum])
print("IWA = " + str(iwa))

print("KL Modes = " + str(list(map(int, sys.argv[3+argnum].split(",")))))
klmodes = list(map(int, sys.argv[3+argnum].split(",")))

print()
print('Planet mask parameters:')

print("Radius = " + str(list(map(int, sys.argv[14+argnum].split(",")))))
ra = list(map(int, sys.argv[14+argnum].split(",")))

print("Position Angle = " + str(list(map(int, sys.argv[15+argnum].split(",")))))
pa = list(map(int, sys.argv[15+argnum].split(",")))

print("Mask width (radial, angular): = " + str(list(map(int, sys.argv[16+argnum].split(",")))))
wid = list(map(int, sys.argv[16+argnum].split(",")))

#object to hold mask parameters for snr map 
mask = (ra, pa, wid)
                
print()

outputFileName = sys.argv[4+argnum]
print("Output FileName = " + outputFileName)


saveSNR = False
if (sys.argv[17+argnum] == 'true' or sys.argv[17+argnum] == 'True'):
    saveSNR = True    

print()





##################################################################
#############                                        #############
#############            PERFORM AUTOMATION          #############
#############                                        #############
################################################################## 

print("reading: " + pathToFiles + "/*.fits")

print()

filelist = glob.glob(pathToFiles + '/*.fits')
dataset = MagAO.MagAOData(filelist)
dataset.IWA = iwa

#creates cube to eventually hold average SNR data
snrCube = np.zeros((len(klmodes),int((annuli2_stop-annuli2_start)/annuli2_inc+1),int((movement2_stop-movement2_start)/movement2_inc+1)))

#creates cube to hold all snr maps 
snrMapCube5d = np.zeros((len(klmodes),int((annuli2_stop-annuli2_start)/annuli2_inc+1),int((movement2_stop-movement2_start)/movement2_inc+1), 450, 450))


#loop over annuli, movement, and subsection parameters

#keeps track of number of annuli values that have been tested, used for indexing
acount = 0

for a in range(annuli2_start, annuli2_stop+1, annuli2_inc):
    
    #keeps track of number of movement values that have been tested, used for indexing
    mcount = 0
    
    for m in np.arrange(movement2_start, movement2_stop+1, movement2_inc):
        
        for s in range(subsections2_start, subsections2_stop+1, subsections2_inc):
            print("Starting KLIP for parameters:")
            print("annuli = %d; movement = %d; subections = %d" %(a, m,s))

            #run klip for given parameters
            parallelized.klip_dataset(dataset, outputdir="", fileprefix=outputFileName, annuli=a, subsections=s, movement=m, numbasis=klmodes, calibrate_flux=True, mode="ADI")
 
            #cube to hold median combinations of klipped images
            cube = np.zeros((len(klmodes),450,450))
            
            #flips images
            dataset.output = dataset.output[:,:,:,::-1]
             
            #keeps track of number of KL mode values that have been tested, used for indexing
            kcount = 0
            
            #iterates over kl modes
            for k in klmodes:
                
                #takes median combination of cube made with given number of KL modes
                isolatedKL = np.nanmedian(dataset.output[kcount,:,:,:], axis=0)
                                
                #adds median image to cube 
                cube[kcount,:,:] = isolatedKL
               
                
                #makes SNR map 
                snrmap = snr.create_map(isolatedKL, planets = mask, saveOutput = False)
                
                #adds SNR map to 5d cube 
                if (saveSNR):
                    snrMapCube5d[kcount,acount,mcount,:,:] = snrmap 
                
                #gets highest pixel value in snr map in the location of the planet 
                planetSNRs = []
                
                for x in range (len(ra)):
                    planetSNRs.append(snr.getPlanet(snrmap, ra[x], pa[x], wid[1]))
                
                planetSNR = np.mean(planetSNRs)
                
                #add planet snr value to snrCube
                snrCube[kcount,acount,mcount] = planetSNR
                kcount+=1
                                          
            #write median combination cube to disk 
            print("Writing median image combinations to " + pathToFiles + "/../")
            writeData(cube, pathToFiles, outputFileName, a, m, s, iwa, klmodes, mask = None, pre = 'med')
            print()
            
        mcount+=1
    acount+=1

if (saveSNR):
    print("Writing 4d SNR data to " + pathToFiles + "/../")
    print()
    #writes SNR maps to 4d cubes 
    for x in range (len(klmodes)):
        snr4d = snrMapCube5d[x,:,:,:,:] 
        writeData(snr4d, pathToFiles, outputFileName, annuli, movement, subsections, iwa, klmodes, mask = mask, pre = 'SNRCube')
         
print("Writing average SNR values to " + pathToFiles + "/../")    
#write snr cube to disk 
writeData(snrCube, pathToFiles, outputFileName, annuli, movement, subsections, iwa, klmodes, mask = mask, pre = 'paramexplore')

print()
print("KLIP automation complete")  
            

