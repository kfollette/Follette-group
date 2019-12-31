#Clare Leonard                                                                
#Version 1.0 - 6/21/17

#Kate Follette
#Version 2.0 - 9/25/19
#Various Modifications.
# 1) Sped up with 3D SNR maps instead of loop through KL modes.
# 2) Added additional parameters to a 5th cube dimension. Now separate maps for peak SNR derived with noise = Standard Deviation,
#peak SNR with noise = median, and the total SNR under the mask for those two cases.
# 3) Runtime printed to terminal
# 4) cleaned up file reading and writing to use fits.getdata and fits.getheader rather than hdulist stuff
# 5) added comments, removed redundancies

##################################################################
#############                                        #############
#############                 IMPORTS                #############
#############                                        #############
##################################################################

import glob
import inspect
import os
import pyklip.instruments.MagAO as MagAO
import pyklip.parallelized as parallelized
import numpy as np
import sys
import pyklip.klip as klip
from astropy.io import fits
import SNRMap_new as snr
import time
import warnings
from astropy.utils.exceptions import AstropyWarning

warnings.filterwarnings('ignore', category=AstropyWarning, append=True)

##################################################################
#############                                        #############
#############               SAVE FILES               #############
#############                                        #############
################################################################## 
 

def writeData(im, prihdr, allParams = False, snrmap = False, pre = ''):
    #function writes out fits files with important info captured in fits headers
    
    if (allParams):
    #for parameter explorer cube output - capture full range of parameter values
        annuli_fname = annuli
        annuli_head = annuli
        movement_fname = movement
        movement_head = movement
        subsections_fname = subsections
        subsections_head = subsections
       
        #if program iterates over several parameter values, formats these for fits headers and file names
        if (isinstance(annuli, tuple)):
            annuli_fname = str(annuli[0]) + '-' + str(annuli[1]) + 'x' + str(annuli[2])
            annuli_head = str(annuli[0]) + 'to' + str(annuli[1]) + 'by' + str(annuli[2])  
        if (isinstance(movement, tuple)):
            movement_fname = str(movement[0]) + '-' + str(movement[1]) + 'x' + str(movement[2])
            movement_head = str(movement[0]) + 'to' + str(movement[1]) + 'by' + str(movement[2])
        if (isinstance(subsections, tuple)):
            subsections_head = str(subsections[0]) + 'to' + str(subsections[1]) + 'by' + str(subsections[2])
            subsections_fname = str(subsections[0]) + '-' + str(subsections[1]) + '-' + str(subsections[2])
    else:
        #for individual images and SNR maps, capture the single parameter values used
        annuli_head = a
        movement_head = m
        subsections_head = s
        annuli_fname = a
        movement_fname = m
        subsections_fname = s

    #shortens file path to bottom 4 directories so it will fit in fits header
    try:
        pathToFiles_short = '/'.join(pathToFiles.split(os.path.sep)[-4:])
    except:
        pathToFiles_short = pathToFiles
            
    #adds info to fits headers
    prihdr['ANNULI']=str(annuli_head)
    prihdr['MOVEMENT']=str(movement_head)
    prihdr['SUBSCTNS']=str(subsections_head)
    prihdr['IWA'] = str(iwa)
    prihdr['KLMODES']=str(klmodes)
    prihdr['FILEPATH']=str(pathToFiles_short)
 
    if(snrmap):
        rad, pa, wid = mask 
        prihdr['MASK_RAD']=str(rad)
        prihdr['MASK_PA']=str(pa)
        prihdr['MASK_WID']=str(wid)
        prihdr['SNRSMTH']=str(_smooth)
        prihdr['SNRFWHM']=str(FWHM)

    if(allParams):
        prihdr["SLICE1"]="average planet peak value under mask in standard deviation noise map"
        prihdr["SLICE2"] = "average planet peak value under mask in median absolute value noise map"
        prihdr["SLICE3"] = "average value of positive pixels under mask in standard deviation noise map"
        prihdr["SLICE4"] = "average value of positive pixels under mask in median absolute value noise map"
        prihdr["SLICE5"] = "total number of pixels >5sigma outside of mask in standard deviation noise map"
        prihdr["SLICE6"] = "total number of pixels >5sigma outside of mask in median absolute value noise map"
        prihdr["SLICE7"] = "total number of pixels >5sigma outside of mask and at similar radius in standard deviation noise map"
        prihdr["SLICE8"] = "total number of pixels >5sigma outside of mask and at similar radius in median absolute value noise map"

    #suff = ''

    #writes out files
    fits.writeto(str(pathToFiles) + "_klip/" + str(pre)  + outputFileName + "_a" + str(annuli_fname) + "m" + str(
        movement_fname) + "s" + str(subsections_fname) + "iwa" + str(iwa) + suff + '_klmodes-all.fits', im, prihdr, overwrite=True)


##################################################################
#############                                        #############
#############       GET INPUTS  FROM GUI             #############
#############                                        #############
################################################################## 

#value adjusts argument numbering in case of white space in file path 
argnum = 0

pathToFiles = sys.argv[1]

#if the file path has white space in it, recognizes the end of the filepath by the phrase '%finish'
#If the phrase '%finish' does not occur, leaves pathToFiles as the first argument
try:
    while (not pathToFiles[-7:] == '%finish'):
        argnum += 1
        pathToFiles = pathToFiles + " " + sys.argv[1+argnum]
    pathToFiles = pathToFiles[:-7]
        
except:
    pathToFiles = sys.argv[1]
    argnum = 0
    
print("File Path = " + pathToFiles)   
print()

#create directory to save ouput to
if not os.path.exists(pathToFiles + "_klip"):
    os.makedirs(pathToFiles + "_klip")
    os.chmod(pathToFiles + "_klip", 0o777)

print("Parameters to explore:")

annuli_start = int(sys.argv[5+argnum])
annuli_stop = int(sys.argv[6+argnum])
annuli_inc = int(sys.argv[7+argnum])
#creates touple for easier eventual string formatting when saving files
annuli = (annuli_start, annuli_stop, annuli_inc)
#if only one parameter is iterated over, makes sure increment is 1 and changes touple to single int
if(annuli_start == annuli_stop):
    annuli_inc = 1;
    annuli = annuli_start
    print("Annuli = " +str(annuli_start))    
else:
    print("Annuli: start = %s; end = %s; increment = %s " %(str(annuli_start), str(annuli_stop), str(annuli_inc)))

movement_start = float(sys.argv[8+argnum])
movement_stop = float(sys.argv[9+argnum])
movement_inc = float(sys.argv[10+argnum])
#creates touple for easier eventual string formatting when saving files
movement = (movement_start, movement_stop, movement_inc)
#if parameter is not set to change, makes sure increment is 1 and changes touple to single int
if(movement_start == movement_stop):
    movement_inc = 1;
    movement = movement_start
    print("Movement = " +str(movement_start))
    
else:
    print("Movement: start = %s; end = %s; increment = %s " %(
        str(movement_start), str(movement_stop), str(movement_inc)))

subsections_start = int(sys.argv[11+argnum])
subsections_stop = int(sys.argv[12+argnum])
subsections_inc = int(sys.argv[13+argnum])
#creates touple for easier eventual string formatting when saving files
subsections = (subsections_start, subsections_stop, subsections_inc)
#if parameter is not set to change, makes sure increment is 1 and changes touple to single int
if(subsections_start == subsections_stop):
    subsections_inc = 1;
    subsections = subsections_start
    print("Subsections = " +str(subsections_start))    
else:
    print("Subsections: start = %s; end = %s; increment = %s " %(
        str(subsections_start), str(subsections_stop), str(subsections_inc)))

print()
    
iwa = int(sys.argv[2+argnum])
print("IWA = " + str(iwa))

print("KL Modes = " + str(list(map(int, sys.argv[3+argnum].split(",")))))
klmodes = list(map(int, sys.argv[3+argnum].split(",")))

FWHM = float(sys.argv[14+argnum])
print("Star FWHM = " + str(FWHM))

_smooth = float(sys.argv[15+argnum])
print("Smoothing Value = " + str(_smooth))

print()
print('Planet mask parameters:')

ra = list(map(int, sys.argv[16+argnum].split(",")))
print("Radius = " + str(ra))

pa = list(map(int, sys.argv[17+argnum].split(",")))
print("Position Angle = " + str(pa))

#catch in case don't specify enough ras
if len(ra) != len(pa):
    print("list of separations is not equal in length to list of position angles. Duplicating to match.")
    ra=np.repeat(ra,len(pa))

wid = list(map(int, sys.argv[18+argnum].split(",")))
print("Mask width (radial, angular): = " + str(wid))

#object to hold mask parameters for snr map 
mask = (ra, pa, wid)
                
print()

outputFileName = sys.argv[4+argnum]
print("Output FileName = " + outputFileName)

saveSNR = False
if (sys.argv[19+argnum] == 'true' or sys.argv[19+argnum] == 'True'):
    saveSNR = True    

suff = ''    
singleAnn = False
if (sys.argv[20+argnum] == 'true' or sys.argv[20+argnum] == 'True'):
    singleAnn = True   
    suff = '_min-annuli'

highpass = False
if (sys.argv[21+argnum] == 'true' or sys.argv[21+argnum] == 'True'):
    highpass = True   
    suff += '_highpass'
    
print()

##################################################################
#############                                        #############
#############            PERFORM AUTOMATION          #############
#############                                        #############
################################################################## 

print("reading: " + pathToFiles + "/*.fits")

start_time = time.time()
print("start clock time is", time.time())

start_process_time = time.process_time()
print("start process time is", time.process_time())

print()

#grab generic header from a generic single image
hdr = fits.getheader(pathToFiles + '/sliced_1.fits')
#erase values that change through image cube
del hdr['ROTOFF']
try:
    del hdr['GSTPEAK']
except:
    print('not a saturated dataset')
del hdr['STARPEAK']


#reads in files
filelist = glob.glob(pathToFiles + '/*.fits')
dataset = MagAO.MagAOData(filelist)

#set iwa and owa
dataset.IWA = iwa
xDim = dataset._input.shape[2]
yDim = dataset._input.shape[1]
owa = min(xDim,yDim)/2

#creates cube to eventually hold parameter explorer data
PECube = np.zeros((8,int((subsections_stop-subsections_start)/subsections_inc+1), len(klmodes),
                    int((annuli_stop-annuli_start)/annuli_inc+1),
                    int((movement_stop-movement_start)/movement_inc+1)))



###BEGIN LOOPS OVER ANNULI, MOVEMENT AND SUBSECTION PARAMETERS

#keeps track of number of annuli values that have been tested, used for indexing
acount = 0

for a in range(annuli_start, annuli_stop+1, annuli_inc):

    #size of annular zones
    dr = float(owa-iwa)/a
    #creates list of zone radii
    all_bounds = [dr*rad+iwa for rad in range(a+1)]
    print('annuli bounds are', all_bounds)
    numAnn = a
    
    if(singleAnn):
        #find maximum annulus boundary radius that is still inside innermost planet injection radius
        lowBound = max([b for b in all_bounds if (min(ra)>b)])
        #find minimum exterior boundary radius that is outside outermost planet injection radius
        upBound = min([b for b in all_bounds if (max(ra)<b)])
        #list of zone boundaries for planets between the two bounds
        all_bounds = [b for b in all_bounds if (b>=lowBound and b<=upBound)]
        numAnn = int(round((upBound-lowBound)/dr))
        #reset iwa and owa to correspond to annulus
        dataset.IWA = lowBound
        dataset.OWA = upBound

    #check to see if any planets fall very close to a zone boundary 
    if (len( [b for b in all_bounds for r in ra if(b <= r+FWHM/2 and b >= r-FWHM/2)] ) == 0):

        #keeps track of number of movement values that have been tested, used for indexing
        mcount = 0

        for m in np.arange(movement_start, movement_stop+1, movement_inc):

            scount = 0

            for s in range(subsections_start, subsections_stop+1, subsections_inc):
                
                if(singleAnn):
                    print("Parameters: movement = %s; subections = %d" %(m,s))
                    print("Running for %d annuli, equivalent to single annulus of width %s pixels" %(annuli_start+acount, dr))
                else:
                    print("Parameters: annuli = %d; movement = %s; subections = %d" %(a, m,s))

                #creates cube to hold snr maps 
                #snrMapCube = np.zeros((2,len(klmodes),yDim,xDim))

                runKLIP = True

                if (os.path.isfile(str(pathToFiles) + "_klip/med_" + outputFileName + "_a" + str(a) + "m" + str(m) + "s" + str(s) + "iwa" + str(iwa) + suff + '_klmodes-all.fits')):
                    print("match")
                    incube = fits.getdata(str(pathToFiles) + "_klip/med_" + outputFileName + "_a" + str(a) + "m" + str(m) + "s" + str(s) + "iwa" + str(iwa) + suff + '_klmodes-all.fits')
                    head = fits.getheader(str(pathToFiles) + "_klip/med_" + outputFileName + "_a" + str(a) + "m" + str(m) + "s" + str(s) + "iwa" + str(iwa) + suff + '_klmodes-all.fits')
                    klmodes2 = head['KLMODES'][1:-1]
                    klmodes2 = list(map(int, klmodes2.split(",")))
                    print(klmodes, klmodes2)

                    if (len([k for k in klmodes if not k in klmodes2]) == 0):
                        print("Found KLIP processed images for same parameters saved to disk. Reading in data.")
                        #don't re-run KLIP
                        runKLIP = False

                if (runKLIP):
                    print("Starting KLIP")
                    #run klip for given parameters
                    parallelized.klip_dataset(dataset, outputdir=(pathToFiles + "_klip/"),
                                              fileprefix=outputFileName, annuli=numAnn, subsections=s, movement=m,
                                              numbasis=klmodes, calibrate_flux=True, mode="ADI", highpass = highpass, time_collapse='median')

                    #collapse in time dimension
                    incube = np.nanmedian(dataset.output, axis=1)
                    #truncates wavelength dimension, which we don't use
                    incube = incube[:,0,:,:]
                    #print('check: input image shape goes from', dataset.output.shape, 'to', incube.shape)

                #list of noise calculation methods
                methods = ['stddev', 'med']

                # makes SNR map
                snrmaps, peaksnr, snrsums, snrspurious= snr.create_map(incube, FWHM, smooth=_smooth, planets=mask, saveOutput=False)
                print(snrspurious.shape)

                #klmode index
                kcount = 0
                # iterates over kl modes
                for k in klmodes:
                    for methodctr in np.arange(2):
                        #loops over planets specified and returns their SNRs
                        planetSNRs = [snr.getPlanet_peak(snrmaps[methodctr,kcount,:,:], ra[x], pa[x], int(FWHM / 2) + 1) for x in range(len(ra))]
                        #print("planet SNRs are", planetSNRs, 'for', methods[methodctr])
                        planetSNR = np.nanmean(planetSNRs)
                        #print("average planet SNR is", planetSNR, 'for', methods[methodctr])

                        #adds peak values from getPlanet_peak to PE cube
                        PECube[methodctr,scount,kcount,acount,mcount] = planetSNR
                    kcount+=1
                    # adds sums under mask from snr.create_map to PE cube
                    PECube[2:4, scount, :, acount, mcount] = snrsums
                    PECube[4:6, scount, :, acount, mcount] = snrspurious[:,:,0]
                    PECube[6:8, scount, :, acount, mcount] = snrspurious[:,:,1]

                if(runKLIP) and np.median(planetSNR)>3:
                    #write median combination cube to disk 
                    print("Writing median image combinations to " + pathToFiles + "_klip/")
                    writeData(incube, hdr, pre = 'med_')

                if (saveSNR) and np.median(planetSNR)>3:
                    print("Writing SNR maps to " + pathToFiles + "_klip/")
                    writeData(snrmaps, hdr, snrmap = True, pre = 'snrmap_')
                print()

                scount+=1
            mcount+=1

    else: 
        print("Planet near annulus boundary; skipping KLIP for annuli = " + str(a))
        print()
        #assign a unique value as a flag for these cases in the parameter explorer map
        PECube[:,:,:,acount,:] = -1000
                
    acount+=1

         
print("Writing parameter explorer file to " + pathToFiles + "_klip/")
#write parameter explorer cube to disk
writeData(PECube, hdr, allParams = True, snrmap = True, pre = 'paramexplore_')



print()
print("KLIP automation complete")

print("end clock time is", time.time())
print("end process time is", time.process_time())
print("total clock runtime: ", time.time()- start_time)
print("total process runtime:", time.process_time()-start_process_time)
