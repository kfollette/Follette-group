import glob
import inspect
import os
import pyklip.instruments.MagAO as MagAO
import pyklip.parallelized as parallelized
import astropy.wcs as wcs
import numpy as np
import sys
import pyklip.klip as klip
from astropy.io import fits
import SNRMap_new as snr
import time
import warnings
import pyklip.fakes as fakes
import multiprocessing as mp
import pdb
import pyklip.instruments.utils.wcsgen as wcsgen
from astropy.utils.exceptions import AstropyWarning
from tqdm import tqdm

warnings.filterwarnings('ignore', category=AstropyWarning, append=True)


def explore_params(path_to_files, outfile_name, iwa, klmodes, annuli_start, annuli_stop, movement_start, 
    movement_stop, FWHM, ra, pa, wid, annuli_inc=1, movement_inc=1, subsections_start=False, subsections_stop=False, subsections_inc=False,  
    smooth=False, input_contrast=False, time_collapse='median', highpass = True, owa=False,
    saveSNR = True, singleAnn = False, boundary=False, verbose = False, snrsmt = False,
    calibrate_flux=False, pickup=True, submean=False):

    """
    I need a docstring

    MODIFICATION NOTES
    2/9/22 - KBF added submean keyword

    """

    #default is 1 subsection
    if subsections_start == False:
        if (subsections_stop != False) or (subsections_inc != False):
            print("must set subsections_start, subsections_stop, and subsections_inc together")
            return()
        subsections_start = 1
        subsections_stop = 1
        subsections_inc = 1

    #pre-klip smooth off = smoothing value of 0
    if smooth == False:
        smooth=0.0
    
    if verbose is True:
        print(f"File Path = {path_to_files}")   
        print()
        print(f"Output Filename = {outfile_name}")
        print("Parameters to explore:")
        print(f"Annuli: start = {annuli_start}; end = {annuli_stop}; increment = {annuli_inc}")
        print(f"Subsections: start = {subsections_start}; end = {subsections_stop}; increment = {subsections_inc} ")
        print(f"Movement: start = {movement_start}; end = {movement_stop}; increment = {movement_inc} ")
        print(f"IWA = {iwa}, KL Modes = {klmodes}, FWHM = {FWHM}, Smoothing Value = {smooth}")
        print()
        print("Planet Parameters")
        print(f"Radius= {ra}, Position Angle = {pa}, Mask Width = {wid}, Input Contrast - {input_contrast}") #, X Positions = {x_positions}, Y Positions = {y_positions} ")
        print()
        print("reading: " + path_to_files + "/*.fits")

    # create directory to save ouput to
    if not os.path.exists(path_to_files + "_klip"):
        os.makedirs(path_to_files + "_klip")
    
    # create tuples for easier eventual string formatting when saving files
    annuli = (annuli_start, annuli_stop, annuli_inc)
    movement = (movement_start, movement_stop, movement_inc)
    subsections = (subsections_start, subsections_stop, subsections_inc)

    # if only one parameter is iterated over, makes sure increment is 1 and changes touple to single int
    if(annuli_start == annuli_stop):
        annuli_inc = 1
        annuli = annuli_start

    # if parameter is not set to change, makes sure increment is 1 and changes touple to single int
    if(movement_start == movement_stop):
        movement_inc = 1
        movement = movement_start

    # if parameter is not set to change, makes sure increment is 1 and changes touple to single int
    if(subsections_start == subsections_stop):
        subsections_inc = 1
        subsections = subsections_start

    # check that position angle and radius lists have the same number of elements
    if len(ra) != len(pa):
        print("List of separations is not equal in length to list of position angles. Duplicating to match.")
        ra=np.repeat(ra,len(pa))

    # object to hold mask parameters for snr map 
    mask = (ra, pa, wid)

    nplanets = len(ra)
    if verbose is True:
        print(nplanets, "planets with separations ", ra, "and PAs ", pa)
    
    # Add suffix to filenames depending on user-specified values
    suff = ''    
    if singleAnn is True:
        suff += '_min-annuli'
    
    if highpass is True:
        suff += '_highpass'

    if type(highpass)!=bool:
        suff+= '_hp'+str(highpass)

    if verbose is True:
    
        print("Reading: " + path_to_files + "/*.fits")
        
        start_time = time.time()
        print("Start clock time is", time.time())
        
        start_process_time = time.process_time()
        print("Start process time is", time.process_time())
        
    # grab generic header from a generic single image
    hdr = fits.getheader(path_to_files + '/sliced_1.fits')

    # erase values that change through image cube
    del hdr['ROTOFF']
    try:
        del hdr['GSTPEAK']
    except:
        print('NOT a saturated dataset')
    del hdr['STARPEAK']
    
    # reads in files
    filelist = glob.glob(path_to_files + '/*.fits')

    # Get star peaks
    starpeak = []
    for i in np.arange(len(filelist)):
        head = fits.getheader(filelist[i])
        starpeak.append(head["STARPEAK"])

    dataset = MagAO.MagAOData(filelist)
    #make a clean copy of the dataset that will be pulled each time (parallelized modifies dataset.input object)
    dataset_input_clean = np.copy(dataset.input)

    palist = sorted(dataset._PAs)
    palist_clean = [pa if (pa < 360) else pa-360 for pa in palist]
    palist_clean_sorted = sorted(palist_clean)
    totrot = palist_clean_sorted[-1]-palist_clean_sorted[0]

    if verbose is True:
        print(f"total rotation for this dataset is {totrot} degrees")

    # set IWA and OWA
    dataset.IWA = iwa

    if owa is False:
        xDim = dataset._input.shape[2]
        yDim = dataset._input.shape[1]
        dataset.OWA = min(xDim,yDim)/2
        owa = dataset.OWA
    else:
        dataset.OWA = owa

    # Make function to write out data 
    def writeData(im, prihdr, annuli, movement, subsections, snrmap = False, pre = ''):
        #function writes out fits files with important info captured in fits headers
        
        #if program iterates over several parameter values, formats these for fits headers and file names
        if (isinstance(annuli, tuple)):
            annuli_fname = str(annuli[0]) + '-' + str(annuli[1]) + 'x' + str(annuli[2])
            annuli_head = str(annuli[0]) + 'to' + str(annuli[1]) + 'by' + str(annuli[2])  
        else: 
            annuli_fname = annuli
            annuli_head = annuli

        if (isinstance(movement, tuple)):
            movement_fname = str(movement[0]) + '-' + str(movement[1]) + 'x' + str(movement[2])
            movement_head = str(movement[0]) + 'to' + str(movement[1]) + 'by' + str(movement[2])
        else: 
            movement_head = movement
            movement_fname = movement

        if (isinstance(subsections, tuple)):
            subsections_head = str(subsections[0]) + 'to' + str(subsections[1]) + 'by' + str(subsections[2])
            subsections_fname = str(subsections[0]) + '-' + str(subsections[1]) + '-' + str(subsections[2])
        else:
            subsections_head = subsections
            subsections_fname = subsections


        #shortens file path to bottom 4 directories so it will fit in fits header
        try:
            path_to_files_short = '/'.join(path_to_files.split(os.path.sep)[-4:])
        except:
            path_to_files_short = path_to_files
                
        #adds info to fits headers
        prihdr['ANNULI']=str(annuli_head)
        prihdr['MOVEMENT']=str(movement_head)
        prihdr['SUBSCTNS']=str(subsections_head)
        prihdr['IWA'] = str(iwa)
        prihdr['KLMODES']=str(klmodes)
        prihdr['FILEPATH']=str(path_to_files_short)
        prihdr['OWA']=str(dataset.OWA)
        prihdr['TIMECOLL']=str(time_collapse)
        prihdr['CALIBFLUX']=str(calibrate_flux)
        prihdr["HIGHPASS"]=str(highpass)

    
        if(snrmap):
            rad, pa, wid = mask 
            prihdr['MASK_RAD']=str(rad)
            prihdr['MASK_PA']=str(pa)
            prihdr['MASK_WID']=str(wid)
            prihdr['SNRSMTH']=str(smooth)
            prihdr['SNRFWHM']=str(FWHM)

        if isinstance(annuli, tuple):
            prihdr["SLICE1"]="planet peak value under mask in standard deviation noise map"
            prihdr["SLICE2"] = "planet peak value under mask in median absolute value noise map"
            prihdr["SLICE3"] = "average value of positive pixels under mask in standard deviation noise map"
            prihdr["SLICE4"] = "average value of positive pixels under mask in median absolute value noise map"
            prihdr["SLICE5"] = "total number of pixels >5sigma outside of mask in standard deviation noise map"
            prihdr["SLICE6"] = "total number of pixels >5sigma outside of mask in median absolute value noise map"
            prihdr["SLICE7"] = "total number of pixels >5sigma outside of mask and at similar radius in standard deviation noise map"
            prihdr["SLICE8"] = "total number of pixels >5sigma outside of mask and at similar radius in median absolute value noise map"
            prihdr["SLICE9"] = "calibrated contrast value of planet/s at a given separation"

        #writes out files
        fits.writeto(str(path_to_files) + "_klip/" + str(pre)  + outfile_name + "_a" + str(annuli_fname) + "m" + str(
            movement_fname) + "s" + str(subsections_fname) + "iwa" + str(iwa) + suff + '-KLmodes-all.fits', im, prihdr, overwrite=True)

        return


    # create cube to eventually hold parameter explorer data
    PECube = np.zeros((9,int((subsections_stop-subsections_start)/subsections_inc+1), len(klmodes), int(nplanets),
                        int((annuli_stop-annuli_start)/annuli_inc+1),
                        int((movement_stop-movement_start)/movement_inc+1)))
    
    #create a file to store paritally completed output
    tmpfile = str(path_to_files) + "_klip/" + "tmp.fits"

    #if restarting from PE that was interrupted
    if (pickup ==True) and (os.path.exists(tmpfile)):
        pe_past = fits.getdata(tmpfile)
    else:
        pe_past = np.zeros(PECube.shape)

    # BEGIN LOOPS OVER ANNULI, MOVEMENT AND SUBSECTION PARAMETERS
    
    # used for indexing: keeps track of number of annuli values that have been tested
    acount = 0
    
    for a in range(annuli_start, annuli_stop+1, annuli_inc):
    
        # calculate size of annular zones
        dr = float(owa-iwa)/a

        # creates list of zone radii
        all_bounds = [dr*rad+iwa for rad in range(a+1)]

        planet_annuli = [a for a in all_bounds if (a<ra[-1]+dr) and (a>ra[0])]
        nplanet_anns = len(planet_annuli)
   
        ann_cen_rad = [ a - dr/2 for a in planet_annuli ]

        if verbose is True:
            print("planets span ", nplanet_anns, "annular zones for annuli = ", a)

        # print('annuli bounds are', all_bounds)
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
    
        #if boundary keyword is set, check to see if any planets are too close to annuli boundaries
        if boundary != False:
            #is planet within +/- number set as boundary pixels
            if not (len( [b for b in all_bounds for r in ra if(b <= r+boundary and b >= r-boundary)] ) == 0):
                print([b for b in all_bounds for r in ra if(b <= r+boundary and b >= r-boundary)])
                print("A planet is near annulus boundary; skipping KLIP for annuli = " + str(a))
                #assign a unique value as a flag for these cases in the parameter explorer map
                PECube[:,:,:,:,acount,:] = np.nan
                #break out of annuli loop before KLIPing
                acount=1
                continue

        # used for indexing: keeps track of number of movement values that have been tested
        mcount = 0
    
        for m in tqdm(np.arange(movement_start, movement_stop+1, movement_inc)):

            #figure out whether there is enough range in the innermost annulus. Skip cases where not enough rotational space.

            if np.arctan(m/(iwa+dr/2))*180/np.pi>totrot:
                if verbose is True:
                    print("movement", m, "=" "%5.1f" % (np.arctan(m/ann_cen_rad[0])*180/np.pi), 
                        "deg. for innermost annulus. Only ", "%5.1f" % (totrot), 
                        "available. skipping this movement/annuli combo") 
                PECube[:,:,:,:,acount,mcount] = np.nan
                mcount+=1
                continue

            else:
                scount = 0
        
                for s in range(subsections_start, subsections_stop+1, subsections_inc):

                    runKLIP = True
                    skip = False

                    klipstr = "_a" + str(a) + "m" + str(m) + "s" + str(s) + "iwa" + str(iwa) 
                    fname  = str(path_to_files) + "_klip/" + outfile_name + klipstr+ suff + '-KLmodes-all.fits'

                    if pickup==True:
                        if np.sum(pe_past[:,scount,:,:,acount,mcount])!=0:
                            PECube[:,scount,:,:,acount,mcount]=pe_past[:,scount,:,:,acount,mcount]
                            skip = True

                    if skip!=True:
                        if verbose is True:  
                            if(singleAnn):
                                print("Parameters: movement = %s; subections = %d" %(m,s))
                                print("Running for %d annuli, equivalent to single annulus of width %s pixels" %(annuli_start+acount, dr))
                            else:
                                print("Parameters: annuli = %d; movement = %s; subections = %d" %(a, m,s))
                
                            # create cube to hold snr maps 
                            #snrMapCube = np.zeros((2,len(klmodes),yDim,xDim))
                        
                        if os.path.isfile(fname):
                            print(outfile_name+klipstr+suff, fname)
                            incube = fits.getdata(fname)
                            head = fits.getheader(fname)
                            klmodes2 = head['KLMODES'][1:-1]
                            klmodes2 = list(map(int, klmodes2.split(",")))
            
                            if (len([k for k in klmodes if not k in klmodes2]) == 0):
                                if verbose is True:
                                    print("Found KLIP processed images for same parameters saved to disk. Reading in data.")
                                #don't re-run KLIP
                                runKLIP = False
            
                        if (runKLIP):
                            if verbose is True:
                                print("Starting KLIP")
                            #run klip for given parameters
                            #read in a fresh copy of dataset so no compound highpass filtering
                            dataset.input = dataset_input_clean
                            parallelized.klip_dataset(dataset, outputdir=(path_to_files + "_klip/"), fileprefix=outfile_name+klipstr+suff, 
                                annuli=numAnn, subsections=s, movement=m, numbasis=klmodes, calibrate_flux=calibrate_flux, 
                                mode="ADI", highpass = highpass, time_collapse=time_collapse, verbose = verbose)

                            #read in the final image and header
                            print(outfile_name+klipstr+suff, fname)
                            #read in file that was created so can add to header
                            incube = fits.getdata(fname)
                            head = fits.getheader(fname)

                            #add KLMODES keyword to header
                            #this also has the effect of giving the file a single header instead of pyklip's double
                            head["KLMODES"]=str(klmodes)
                            fits.writeto(fname, incube, head, overwrite=True)
                        
                        if input_contrast is not False:
                            dataset_copy = np.copy(incube)
                        
                            # Find planet x and y positions from pa and sep
                            x_positions = [r*np.cos((np.radians(p+90)))+ dataset.centers[0][0] for r, p in zip(ra, pa)]
                            y_positions = [r*np.sin((np.radians(p+90)))+ dataset.centers[0][0] for r, p in zip(ra, pa)]
                        
                            # Loop through kl modes
                            cont_meas = np.zeros((len(klmodes), 1))
                            for k in range(len(klmodes)):
                            
                                dataset_contunits = dataset_copy[k]/np.median(starpeak)
                                                           
                                if runKLIP is False:
                                    w = wcsgen.generate_wcs(parangs = 0, center = dataset.centers[0])
                                else:
                                    w = dataset.output_wcs[0]
                                
                                # Retrieve flux of injected planet
                                planet_fluxes = []
                                for sep, p in zip(ra, pa):
                                    fake_flux = fakes.retrieve_planet_flux(dataset_contunits, dataset.centers[0], w, sep, p, searchrad=7)
                                    planet_fluxes.append(fake_flux)
                          
                                # Calculate the throughput
                                tpt = np.array(planet_fluxes)/np.array(input_contrast)
                                
                                # Create an array with the indices are that of KL mode frame with index 2
                                ydat, xdat = np.indices(dataset_contunits.shape)
                           
                                # Mask the planets
                                for x, y in zip(x_positions, y_positions):

                                    # Create an array with the indices are that of KL mode frame with index 2
                                    distance_from_star = np.sqrt((xdat - x) ** 2 + (ydat - y) ** 2)

                                    # Mask
                                    dataset_contunits[np.where(distance_from_star <= 2 * FWHM)] = np.nan
                                    masked_cube = dataset_contunits

                                # Measure the raw contrast
                                contrast_seps, contrast = klip.meas_contrast(dat=masked_cube, iwa=iwa, owa=dataset.OWA, resolution=(7), center=dataset.centers[0], low_pass_filter=True)

                                # Find the contrast to be used 
                                use_contrast = np.interp(np.median(ra), contrast_seps, contrast)
                                
                                # Calibrate the contrast
                                cal_contrast = use_contrast/np.median(tpt)
                                cont_meas[k] = -cal_contrast
                                        
                        # makes SNR map
                        snrmaps, peaksnr, snrsums, snrspurious= snr.create_map(fname, FWHM, smooth=snrsmt, planets=mask, saveOutput=False, sigma = 5, checkmask=False, submean=submean, verbose = verbose)

                        PECube[0:2, scount, :, :, acount, mcount] = peaksnr
                        PECube[2:4, scount, :, :, acount, mcount] = snrsums
                        PECube[4:6, scount, :, :, acount, mcount] = snrspurious[:,:,None,0]
                        PECube[6:8, scount, :, :, acount, mcount] = snrspurious[:,:,None,1]
                        PECube[8, scount, :, :, acount, mcount] = cont_meas

                        pe_past[:,scount,:,:,acount,mcount]=PECube[:,scount,:,:,acount,mcount]
                        fits.writeto(tmpfile,pe_past, overwrite=True)

                        if(runKLIP) and np.nanmedian(peaksnr)>3:
                            writeData(incube, hdr, a, m, s)
                        if verbose is True:
                            print("Median peak SNR > 3. Writing median image combinations to " + path_to_files + "_klip/")
                            
                        if saveSNR is True:
                            writeData(snrmaps, hdr, a, m, s, snrmap = True, pre = 'snrmap_')
                            if verbose is True:
                                print("Writing SNR maps to " + path_to_files + "_klip/")
        
                
                    scount+=1
                mcount+=1                
        acount+=1
    
    if verbose is True:        
        print("Writing parameter explorer file to " + path_to_files + "_klip/")

    #write parameter explorer cube to disk
    writeData(PECube, hdr, annuli, movement, subsections, snrmap = True, pre = 'paramexplore_')

    if verbose is True: 
        print("KLIP automation complete")    
        print("End clock time is", time.time())
        print("End process time is", time.process_time())
        print("Total clock runtime: ", time.time()- start_time)
        print("Total process runtime:", time.process_time()-start_process_time)

    return(PECube)


   