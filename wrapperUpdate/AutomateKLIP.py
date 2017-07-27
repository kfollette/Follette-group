
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



##################################################################
#############                                        #############
#############               GET INPUTS               #############
#############                                        #############
################################################################## 

#value adjusts argument numbering in case of white space in file path 
argnum = 0


pathToFiles = sys.argv[1]

#if filepath doesnt end in sliced, sontinues to add next arguements, helpful iin case of whitespace in file path
while (not pathToFiles[-1] == 'd' and not pathToFiles[-1] == '"'):
    argnum += 1
    pathToFiles = pathToFiles + " " + sys.argv[1+argnum]
    
print("File Path = " + pathToFiles)

print()

print("Parameters to explore:")

annuli2_start = int(sys.argv[5+argnum])
annuli2_stop = int(sys.argv[6+argnum])
annuli2_inc = int(sys.argv[7+argnum])

if(annuli2_start == annuli2_stop):
    annuli2_inc = 1;
    print("Annuli = " +str(annuli2_start))
    
else:
    print("Annuli: start = %s; end = %s; increment = %s " %(str(annuli2_start), str(annuli2_stop), str(annuli2_inc)))



movement2_start = int(sys.argv[8+argnum])
movement2_stop = int(sys.argv[9+argnum])
movement2_inc = int(sys.argv[10+argnum])


if(movement2_start == movement2_stop):
    movement2_inc = 1;
    print("Movement = " +str(movement2_start))
    
else:
    print("Movement: start = %s; end = %s; increment = %s " %(str(movement2_start), str(movement2_stop), str(movement2_inc)))


subsections2_start = int(sys.argv[11+argnum])
subsections2_stop = int(sys.argv[12+argnum])
subsections2_inc = int(sys.argv[13+argnum])


if(subsections2_start == subsections2_stop):
    subsections2_inc = 1;
    print("Subsections = " +str(subsections2_start))
    
else:
    print("Subsections: start = %s; end = %s; increment = %s " %(str(subsections2_start), str(subsections2_stop), str(subsections2_inc)))

print()
    
iwa = int(sys.argv[2+argnum])
print("IWA = " + str(iwa))

print()

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

print()

outputFileName = sys.argv[4+argnum]
#outputFileName = nameOutput(pathToFiles)
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


snrCube = np.zeros((len(klmodes),int((annuli2_stop-annuli2_start)/annuli2_inc+1),int((movement2_stop-movement2_start)/movement2_inc+1)))

snrMapCube5d = np.zeros((len(klmodes),int((annuli2_stop-annuli2_start)/annuli2_inc+1),int((movement2_stop-movement2_start)/movement2_inc+1), 450, 450))


#loop over annuli, movement, and subsection parameters

print("running klip for parameters:")


#keeps track of number of annuli values that have been tested, used for indexing
acount = 0

for a in range(annuli2_start, annuli2_stop+1, annuli2_inc):
    
    #keeps track of number of movement values that have been tested, used for indexing
    mcount = 0
    
    for m in range(movement2_start, movement2_stop+1, movement2_inc):
        
        for s in range(subsections2_start, subsections2_stop+1, subsections2_inc):
            print("annuli = %d; movement = %d; subections = %d" %(a, m,s))
            #sys.stdout.flush()

            #run klip for given parameters
            parallelized.klip_dataset(dataset, outputdir="", fileprefix=outputFileName, annuli=a, subsections=s, movement=m, numbasis=klmodes, calibrate_flux=True, mode="ADI")

            
            #cube to hold median combinations of klipped images
            cube = np.zeros((len(klmodes),450,450))
            
            #keeps track of number of KL mode values that have been tested, used for indexing
            kcount = 0
            
            #flips images
            dataset.output = dataset.output[:,:,:,::-1]
            
            
            #iterates over kl modes
            for k in klmodes:
                
                #takes median combination of cube made with given number of KL modes
                isolatedKL = np.nanmedian(dataset.output[kcount,:,:,:], axis=0)
                
                #put together output name
                outputNameSNR = outputFileName + "_a" + str(a) + "m" + str(m) + "KL" + str(k) + "_SNRMap.fits"
                
                #object to hold mask parameters for snr map 
                mask = (ra, pa, wid)
                
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
                
                #adds median image to cube 
                cube[kcount,:,:] = isolatedKL
                
                #add planet snr value to snrCube
                snrCube[kcount,acount,mcount] = planetSNR
                kcount+=1
                                          
            #write median combination cube to disk 
            fits.writeto(pathToFiles + '/../' + 'med_'+ outputFileName + "_a" + str(a) + "m" + str(m) + "s" + str(subsections2) + "iwa" + str(iwa) +'-KLmodes-all.fits', cube, clobber=True)
        mcount+=1
    acount+=1

if (saveSNR):
    #writes SNR maps to 4d cubes 
    for x in range (len(klmodes)):
        snr4d = snrMapCube5d[x,:,:,:,:] 
        fits.writeto(pathToFiles + '/../' + outputFileName + "_a" + str(annuli2_start) + "-" + str(annuli2_stop) + "x" + str(annuli2_inc) + "m" + str(movement2_start) + "-" + str(movement2_stop) + "x" + str(movement2_inc) + "_" + str(klmodes[x]) + 'KLmodes_SNRMaps.fits', snr4d)  
    
#write snr cube to disk 
fits.writeto(pathToFiles + '/../' + outputFileName + "_paramexplore_a" + str(annuli2_start) + "-" + str(annuli2_stop) + "x" + str(annuli2_inc) + "m" + str(movement2_start) + "-" + str(movement2_stop) + "x" + str(movement2_inc) + '-KLmodes-all.fits', snrCube)     
            

