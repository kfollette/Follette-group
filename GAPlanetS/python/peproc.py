import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import pandas as pd
import gaplanets as gp
from importlib import reload
import re
from scipy import ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import glob
import pyklip.klip as klip
import pyklip.instruments.MagAO as MagAO
import pyklip.parallelized as parallelized
import SNRMap_new as snr
import os
import pdb
import pickle
import textwrap
import threading
import itertools as it
import progressbar
from time import sleep
from datetime import date


def collapse_planets(pename, pedir='./', outdir='proc/', writestr=False, writefiles=True, snrthresh=False, oldpe=False, separate_planets=False):
    """
    Averages over the planet dimension of a parameter explorer file

    REQUIRED INPUTS:
    pename:     name of paramexplore file

    OPTIONAL INPUTS
    pedir:      directory holding parameter explorer file
    writestr:   filename prefix for saved file  
                If not specified, preserves name of parameter explorer
    snrthresh:  SNR threshhold for peak pixel under mask. All lower values will be masked.
 
   RETURNS:
    pecube:     parameter explorer matrix with planet dimension collapsed to 1
    writename:  name of output file

    """

    if writestr == False:
        #use the parameter explorer name for this file as well, minus the '_highpass_klmodes-all.fits'
        writestr = pename[:-17]

    # read in image and header
    pecube = fits.getdata(pedir + pename)
    pehead = fits.getheader(pedir + pename)
    dims = pecube.shape
    nplanets=dims[3]
    pehead["NPLANET"]=nplanets
    pehead["PLSEP"]=separate_planets

    #if snrthresh set, find values where peak pixel SNRs (slices 1/2) are below threshhold and replace with nans
    if (snrthresh != False):
        #stdev maps
        ind = np.where(pecube[0,:,:,:,:,:]<snrthresh)
        lowsnr_mask=np.ones(dims[1:])
        lowsnr_mask[ind]=np.nan
        #apply to peak (sl 0) and avg under mask (sl 2)
        for sl in [0,2]:
            pecube[sl,:,:,:,:]*=lowsnr_mask
        
        #absmed maps
        ind_absmed = np.where(pecube[1,:,:,:,:,:]<snrthresh)
        lowsnr_mask_absmed=np.ones(dims[1:])
        lowsnr_mask_absmed[ind_absmed]=np.nan
        #apply to peak (sl 1) and avg under mask (sl 3)
        for sl in [1,3]:
            pecube[sl,:,:,:,:]*=lowsnr_mask_absmed

    #normalize each planet by dividing by its maximum SNR across all movement and annuli
    #does NOT average over KL modes or metric
    #planet loop
    #for i in np.arange(dims[-1]):
        #klloop
        #for j in np.arange(dims[2]):
            #slice loop
            #for k in np.arange(4):
                #normalize by dividing by max for this planet, klmode, slice (snr slices only)
                #pecube[k,:,j,:,:,i]/=np.nanmax(pecube[k,:,j,:,:,i])
    
    #sort of a sloppy fix to make sure spurious pixel metric is NaN where no refims
    #make a mask from one of the SNR metric slices, which have proper NaNs
    if np.ndim(pecube) == 5:
        nanmask = np.where(np.isnan(pecube[3,:,:,:,:])==True)
    elif np.ndim(pecube) == 6:
        nanmask = np.where(np.isnan(pecube[3,:,:,:,:,:])==True)
    #mask with ones everywhere except these locations
    msk = np.ones(dims[1:])
    msk[nanmask]=np.nan

    #for spurious pixel and contrast metrics, multiply by this mask
    for k in np.arange(4,9):
        pecube[k,:,:,:,:,:]*=msk
   
    #create directory to save in if doesn't yet exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #collapse in planet dimension or extract 
    #mean preserves nans so all planets have to be recovered
    #separates planets
    if separate_planets==True: 
        npldim=nplanets
        writename=writestr+'_sepplanets.fits'
    
    #collapses planet dimension
    else: 
        writename = writestr+'_planetcollapse.fits'
        #quick hack for old PEs where planet dimension was 6th
        if oldpe==True:
            pecube=np.mean(pecube,axis=5, keepdims=True)
        else:
            pecube=np.mean(pecube,axis=3, keepdims=True)
        npldim=1
    if writefiles==True:
        fits.writeto(outdir+writename, pecube, pehead, overwrite=True)

    return (pecube, pehead, writename, npldim)

def trimkl(pename, kllist, pedir='./', outdir='proc/', writestr=False, hdr=False, writename=False, writefiles=True):
    """
    Reads in a parameter explorer file and collapses it in the KLmode dimension (axis 3)

    REQUIRED INPUTS:
    pename:     name of paramexplore file
    kllist:     list of KL modes you want to collapse over

    OPTIONAL INPUTS:
    pedir:      directory holding parameter explorer file
    writestr:   filename prefix for saved file 
                If not specified, preserves name of parameter explorer

    RETURNS:

    """

    #if provided string, read in file
    if type(pename) is str:
        if writestr == False:
            writestr = pename[:-5]

        writename=writestr + '_trimmedkls.fits'

        # read in image and header
        klcube_raw = fits.getdata(pedir + pename)
        head = fits.getheader(pedir + pename)

    #if provided array, carry over
    else:
        klcube_raw = pename
        head = hdr
        writename = writename + '_trimmedkls.fits'

    dims=list(klcube_raw.shape)

    #create an array with same dims except for kl modes
    
    dims[2]=len(kllist)
    klkeep = np.zeros(dims)
        
    # pull KL modes of parameter explorer from the header
    allkl = list(map(int, head['KLMODES'][1:-1].split(",")))

    # find only the KL slices in the list of modes to be collapsed and make an array
    i = 0;
    j = 0;
    keepind = []
    for kl in allkl:
        if kl in kllist:
            keepind.append(i)
            #print('keeping kl', kl, "with index", i)
            #only keep KL modes matching kllist
            klkeep[:, :,j,:,:,:] = klcube_raw[:,:,i,:,:,:]
            j += 1
        i += 1

    # replacing -250 (value when planet too close to annulus) with 0 for ease of display
    klkeep[klkeep == -1000] = np.nan

    # update header to reflect KL modes used
    head["KLMODES"] = str(kllist)

    if writefiles==True:
        fits.writeto(outdir+writename, klkeep, head, overwrite=True)

    #return trimmed cubes
    return (klkeep, head, writename)


def filter_nan_gaussian_conserving(arr, sigma):
    """Apply a gaussian filter to an array with nans.

    Intensity is only shifted between not-nan pixels and is hence conserved.
    The intensity redistribution with respect to each single point
    is done by the weights of available pixels according
    to a gaussian distribution.
    All nans in arr, stay nans in gauss.
    """
    nan_msk = np.isnan(arr)

    loss = np.zeros(arr.shape)
    loss[nan_msk] = 1
    loss = ndimage.gaussian_filter(
            loss, sigma=sigma, mode='constant', cval=1)

    gauss = arr.copy()
    gauss[nan_msk] = 0
    gauss = ndimage.gaussian_filter(
            gauss, sigma=sigma, mode='constant', cval=0)
    gauss[nan_msk] = np.nan

    gauss += loss * arr

    return gauss

def find_best_new(pename, kllist, pedir='./', writestr=False, writefiles=True, weights=[1,1,1,1,1,1], outdir='proc/', snrthresh=False,
    oldpe=False, debug=False, smt=3, snrmeth='all',separate_planets=False, separate_kls=False, verbose=False, maxy=25, maxx=25):
    """
    collapses parameter explorer file and extracts the optimal parameter value

    REQUIRED INPUTS:
    pename:     name of paramexplore file
    kllist:     list of KL modes you want to collapse over
   
    OPTIONAL INPUTS:
    pedir:      directory holding parameter explorer file
    weights:    weights for parameter quality metric. Metrics are listed in order below. 
                [peak SNR, peak SNR neighbors, avg SNR, avg SNR neighbors, stdev, stdev neighbors, spurious pixels]
    smt:        FWHM of gaussian smoothing kernel for 'neighbor quality' metrics
    snrmeth:    one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value),
                and 'all' to average the two


    RETURNS:


    """

    #EXTRACT PLANETS OR COLLAPSE

    if verbose==True:
        if separate_planets==False:
            print("COLLAPSING IN PLANET DIMENSION")
        else:
            print("EXTRACTING PLANETS SEPARATELY")

    pecube, pehead, plwritename, npldim = collapse_planets(pename, pedir=pedir, outdir=outdir, snrthresh=snrthresh, oldpe=oldpe, writestr=writestr, writefiles=writefiles, separate_planets=separate_planets)

    #EXTRACT KL MODES OR COLLAPSE
    if verbose==True:
        print("EXTRACTING ONLY KL MODES SPECIFIED")
    if writefiles==True:
        kltrim, head, writename = trimkl(plwritename, kllist, pedir=outdir, outdir=outdir)
    else:
        kltrim, head, writename = trimkl(pecube, kllist, pedir=outdir, outdir=outdir, hdr=pehead, writename=plwritename, writefiles=writefiles)

    if writestr==False:
        writestr=writename[:-5]

    # if collapsing, make mean and stdev arrays
    if separate_kls==False:
        if verbose==True:
            print("COLLAPSING IN KL DIMENSION")
        #stdevkl = np.std(kltrim, axis=2, keepdims=True)
        #sumkl = np.sum(kltrim, axis=2, keepdims=True)
        #overwrite kltrim with average
        kltrim= np.mean(kltrim, axis=2, keepdims=True)
        #print(stdevkl.shape)

        #grab header
        #head = fits.getheader(outdir+writename)
        head["KLCOLL"]='True'

        # write arrays
        fits.writeto(outdir+ writestr + '_avgkl.fits', kltrim, head, overwrite=True)
        #fits.writeto(outdir+ writestr + '_stdevkl.fits', stdevkl, head, overwrite=True)
        #fits.writeto(outdir+writestr + '_sumkl.fits', sumkl, head, overwrite=True)

    ##EXTRACT SNR SLICES
    #pull the SNR map slices with the matching snrmeth value
    if snrmeth == "absmed":
        slice = 1
    if snrmeth == "stdev":
        slice = 0  

    dims = kltrim.shape

    #extract the appropriate slices
    if verbose==True:
        print("EXTRACTING", snrmeth, "SNR SLICES FROM PE CUBE")
    if snrmeth != 'all':
        #for single method cubes, slices are SNR peak, avg SNR, total >thresh pixels, >thresh pixels inside CR
        kltrim_snr=kltrim[slice::2,:,:,:,:,:]

        #if snrmeth = absmed, grab contrast slice too
        if slice==1:
            dmns = np.array(kltrim_snr.shape)
            dmns[0]+=1
            kltrim_plusone = np.zeros(dmns)
            kltrim_plusone[0:-1,:,:,:,:,:]=kltrim_snr
            #fill last one with contrast cube dimension, which is not dependent on snr method
            kltrim_plusone[-1,:,:,:,:,:]=kltrim[-1,:,:,:,:,:]
            kltrim_snr = kltrim_plusone

    #if snrmeth is 'all', average the slices
    else:
        #average each of the first 4 pairs of metrics
        avgdim = list(kltrim.shape)
        avgdim[0]=5
        kltrim_snr=np.zeros(avgdim)
        for i in np.arange(4):
            kltrim_snr[i,:,:,:,:,:]=np.average(kltrim[2*i:2*i+1,:,:,:,:,:], axis=0)
        kltrim_snr[-1,:,:,:,:,:]=kltrim[-1,:,:,:,:,:]

    #set up kllist argument for loop
    if len(kllist)>1 and separate_kls==False:
        klloop = 1
        kllist=[kllist]
        #stdev_valid=True
    else:
        if verbose==True: 
            print("EXTRACTING SLICES FOR KL MODES", kllist, 'SEPARATELY')
        klloop = len(kllist)
        #stdev_valid=False

    #set up cubes for planet, kl loop

    #extracts parameter explorer inputs from the filename (assumes was autonamed by parameter explorer)
    namelist = pename.split('_')
    params = [s for s in namelist if s[0] == 'a']
    params = params[0]
    params = re.split('a|-|x|m|s|iwa', params)
    ymin = int(params[1])
    ymax = int(params[2])
    ystep = int(params[3])
    xmin = int(float(params[4]))
    xmax = int(float(params[5]))
    xstep = int(float(params[6]))
    nstepx = int((xmax - xmin) / xstep) + 1
    nstepy = int((ymax - ymin) / ystep) + 1

    #sets up arrays for filling
    nq_snr = np.zeros([nstepy, nstepx])
    #nq_stdev = np.zeros([nstepy, nstepx])
    agg_cube=np.zeros([klloop,npldim,nstepy,nstepx])
    agg=np.zeros([nstepy,nstepx])
    ann_val = np.zeros([klloop,npldim])
    movm_val = np.zeros([klloop,npldim])
    ##Note - hard coded for 1 subsection. 
    metric_cube = np.zeros([7, 1, klloop, npldim, nstepy, nstepx])
    qual_cube = np.zeros([6, 1, klloop, npldim, nstepy, nstepx])

    for p in np.arange(npldim):

        for k in np.arange(klloop):
                
            #finds locations of peaks
            #note - currently hard coded for 1 subsection (second index)
            maxind = np.where(kltrim_snr[0,0,k,p,:,:] == np.nanmax(kltrim_snr[0,0,k,p,:,:]))
            maxind_umask = np.where(kltrim_snr[1,0,k,p,:,:] == np.nanmax(kltrim_snr[1,0,k,p,:,:]))

            #translates to x and y coordinates
            xcoord = maxind[0][0]
            ycoord = maxind[1][0]
            #print("peak value for peak SNR is at coordinates:", xcoord, ycoord)

            #translates to x and y coordinates for umask
            xcoord_umask = maxind_umask[0][0]
            ycoord_umask = maxind_umask[1][0]
            #print("peak value for avg SNR under mask is at coordinates:", xcoord, ycoord)

            #normalize the SNR metric to range from 0 to 1
            #note - hard-coded for susections = 1
            snr_norm = (kltrim_snr[0,0,k,p,:,:]-np.nanmin(kltrim_snr[0,0,k,p,:,:]))/ np.nanmax(kltrim_snr[0,0,k,p,:,:])
            snr_norm_umask = (kltrim_snr[1,0,k,p,:,:]-np.nanmin(kltrim_snr[1,0,k,p,:,:])) / np.nanmax(kltrim_snr[1,0,k,p,:,:])

            #if stdev_valid==True:
                #normalize standard deviations across KL modes. Low values = good
                # divide by SNR so is Stdev in SNR as fraction of SNR itself
                #stdev_norm_cube = stdevkl[0:2,0,k,p,:,:] / kltrim_snr[0:2,0,k,p,:,:]
                #first slice is peak
                #stdev_norm = 1 - ((stdev_norm_cube[0,:,:]-np.nanmin(stdev_norm_cube[0,:,:]))/np.nanmax(stdev_norm_cube[0,:,:]))
                #second slice is under mask
                #stdev_norm_umask = 1 - ((stdev_norm_cube[1,:,:]-np.nanmin(stdev_norm_cube[1,:,:]))/np.nanmax(stdev_norm_cube[1,:,:]))

            #if extracting separately, fill these arrays with nans 
            #else:
                #stdev_norm = np.zeros([nstepy,nstepx])*np.nan
                #stdev_norm_umask = np.zeros([nstepy,nstepx])*np.nan

            #spurious pixels metrics - pulling slice 3 (in between IWA and CR only)
            spurpix = kltrim_snr[3,0,k,p,:,:]
            
            #make a contrast metric
            #returned quantity is -1*contrast. turn back into contrast
            #and log so that higher = better
            logcontrast = np.log10(-1*kltrim_snr[4,0,k,p,:,:])
            #filter out unphysical contrasts
            logcontrast[logcontrast>0]=np.nan
            #now take absolute value - smaller is better
            logcontrast = np.abs(logcontrast)
            #and subtract the minimum so goes min=0 to max
            logcontrast = logcontrast - np.nanmin(logcontrast)
            #now divide by the max so goes 0-->1
            contrast = logcontrast/np.nanmax(logcontrast)
            
            #compute neighbor quality metrics by smoothing with Gaussian
            sig=smt
            nq_snr = filter_nan_gaussian_conserving(snr_norm,sig)
            #nq_stdev = filter_nan_gaussian_conserving(stdev_norm,sig)
            nq_snr_umask = filter_nan_gaussian_conserving(snr_norm_umask,sig)
            #nq_stdev_umask = filter_nan_gaussian_conserving(stdev_norm_umask,sig)

            #normalizes neighbor quality
            nq_snr /= np.nanmax(nq_snr)
            #nq_stdev /= np.nanmax(nq_stdev)
            nq_snr_umask /= np.nanmax(nq_snr_umask)
            #nq_stdev_umask /= np.nanmax(nq_stdev_umask)

            if debug==True:
                #make a cube of all these metrics for sanity checking
                qual_cube[0,:,k,p,:,:]=snr_norm
                qual_cube[1,:,k,p,:,:]=snr_norm_umask
                #qual_cube[2,:,k,p,:,:]=stdev_norm
                #qual_cube[3,:,k,p,:,:]=stdev_norm_umask
                qual_cube[2,:,k,p,:,:]=spurpix
                qual_cube[3,:,k,p,:,:]=nq_snr
                qual_cube[4,:,k,p,:,:]=nq_snr_umask
                #qual_cube[7,:,k,p,:,:]=nq_stdev
                #qual_cube[8,:,k,p,:,:]=nq_stdev_umask
                qual_cube[5,:,k,p,:,:]=contrast

            #average under mask and peak pixel estimates
            #snr_norm_combo = (snr_norm_avg + snr_norm_avg_umask) / 2.
            #nq_snr_combo = (nq_snr + nq_snr_umask) / 2.
            #stdev_norm_combo = (stdev_norm_avg + stdev_norm_avg_umask) / 2.
            #nq_stdev_combo = (nq_stdev + nq_stdev_umask) / 2.

            #write out the cubes being used for the final metric
            metric_cube[0,:,k,p,:,:]= snr_norm
            metric_cube[1,:,k,p,:,:]= nq_snr
            metric_cube[2,:,k,p,:,:]= snr_norm_umask
            metric_cube[3,:,k,p,:,:]= nq_snr_umask

            #going to use the under mask values for the final answer - probably more stable
            #metric_cube[4,:,k,p,:,:]= stdev_norm_umask
            #metric_cube[5,:,k,p,:,:]= nq_stdev_umask

            #spurious pixel metric = 1 if no spurious pixels and 0 if max number for this dataset
            if np.nanmax(spurpix)>0:
                spurpix_norm = 1-(spurpix-np.nanmin(spurpix))/np.nanmax(spurpix)
            else: #edge case - no spurious pixels in any image
                spurpix_norm= 1+spurpix
             
            metric_cube[4,:,k,p,:,:]= spurpix_norm
            metric_cube[5,:,k,p,:,:]= contrast 

            metriclist = (snr_norm, nq_snr, snr_norm_umask, nq_snr_umask, spurpix_norm, contrast)
            #metriclist = (snr_norm, nq_snr, snr_norm_umask, nq_snr_umask, stdev_norm_umask, nq_stdev_umask, spurpix_norm, contrast)
            
            #make sure weights for stdev slices are 0 if only 1 kl mode or extracting separately
            #if stdev_valid==False:
                #weights[4]=0
                #weights[5]=0

            #calculate an aggregate parameter quality metric by summing the individual metrics * their weights
            agg=np.zeros((nstepy,nstepx))
            for metricind in np.arange(len(metriclist)):
                #only sum if non-zero weight (otherwise nans in unweighted metrics will propagate)
                if weights[metricind]>0:
                    agg+=weights[metricind]*metriclist[metricind]

            metric_cube[6,:,k,p,:,:]=agg
            agg_cube[k,p,:,:]=agg

            ##find location or peak of parameter quality metric and print info
            maxx=int(maxx)
            maxy=int(maxy)
            ind = np.where(agg[:maxy,:maxx] == np.nanmax(agg[:maxy,:maxx]))

            if agg[ind[0],ind[1]].shape[0]>1:
                print("the optimal solution for this choice of parameters/weights is not unique")
                return()

            #extract metric scores for this location    
            #metric_scores = [snr_norm[ind][0], nq_snr[ind][0], snr_norm_umask[ind][0], nq_snr_umask[ind][0], \
            #stdev_norm_umask[ind][0], nq_stdev_umask[ind][0], spurpix_norm[ind][0], contrast[ind][0], agg[ind][0]]
            metric_scores = [snr_norm[ind][0], nq_snr[ind][0], snr_norm_umask[ind][0], nq_snr_umask[ind][0], spurpix_norm[ind][0], contrast[ind][0], agg[ind][0]]


            #translate to annuli and movement values
            ann_val[k,p]= ymin + ind[0][0] * ystep
            movm_val[k,p] = xmin + ind[1][0] * xstep
            
            if separate_planets==False:
                plno = 'all'
            else:
                plno = p+1
            if verbose==True:
                print('peak for planet =', plno, 'klmode = ', kllist[k], 'is at', ind[0][0], ind[1][0], 'corresponding to annuli', ann_val[k][p], ' and movement', movm_val[k][p])
            #print('SNR value for fake planets (avg of SNR methods and planets) is', avgSNR)
            #print('metric scores for (snr peak, snr peak neigbors, snr umask, snr umask neighbors, stdev, stdev neighbors, spurious pix, contrast, agg) are:', metric_scores)



    if debug==True:
        fits.writeto(outdir+writename[:-5]+'_paramqual_fullcube.fits', qual_cube, overwrite=True)

    metric_fname = writename[:-5]+'_paramqual_metrics.fits'

    if writefiles == True:
        fits.writeto(outdir+metric_fname, metric_cube, overwrite=True)
    
    return metric_cube, agg_cube, ann_val, movm_val, metric_scores, metric_fname


def collapse_pes(pedir='./', kllist=[5,10,20,50], wts = [1,1,1,1,1,1], mode='Line', 
                snrmeth='stdev', snrthresh=False, outdir='proc/', xname='', 
                datadir='../', header=True, smt=3, savefig=True,
                hpval=None, collmode=None, owa=None, oldpe=False, calflux=False, 
                separate_planets=False, separate_kls=False):
    """
    Collapses ALL parameter explorer files in a given directory according to the specified combination of KL modes,
    metric weights, and SNR computation method and runs the corresponding KLIP reductions for that set of parameters

    OPTIONAL INPUTS:
    pedir:      directory holding parameter explorer files
    kllist:     kl modes to keep. will average over them if more than one is specified.
    weights:    weights for parameter quality metric. Metrics are listed in order below. 
                [peak SNR, peak SNR neighbors, avg SNR, avg SNR neighbors, stdev, stdev neighbors, spurious pixels]
    mode:       what type of images to run the KLIP reductions on. Options are Line, Cont, or SDI (which will generate
                both and subtract them)
    snrmeth:    one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value),
                and 'all' to average the two
    snrthresh:  SNR threshhold for peak pixel under mask. All lower values will be masked.
    outdir:     location to store output files
    xname:      additional string to add to all filenames 
    datadir:    path to the group's dropbox on your machine (note you can selective sync only relevant dq_cuts subdirectories)
    hpval:      value of the highpass keyword for KLIP. can be True/False or value. 
    collmode:   value for time_collapse keyword for KLIP. Options are 'median', 'mean', 'weighted-mean'
    oldpe:      set to True for older parameter explorers where planet dimension is axis=5 rather than axis=3
    header:     if True, will ignore specified values for owa, collmode, hpval, and calflux and get them from PE header

    RETURNS:
    d:          a dictionary with all relevant parameter explorer info, KLIP parameters, and KLIPed images


    """


    #set up file name
    klliststr=str(kllist).replace(',','_')[1:-1]
    klliststr=klliststr.replace(" ","")
    wtstr = re.sub(r'\W+', '', str(wts))
    xstr = '_'+wtstr+'_'+snrmeth+'_'+klliststr+xname

    #find all the parameter explorer files in the specified directory
    rundir = os.getcwd()
    os.chdir(pedir)
    flist = glob.glob('paramexplore*KLmodes-all.fits')
    
    os.chdir(rundir)

    #sort alphabetically so the same objects are together
    flist = sorted(flist)

    #set up dictionary
    npes=len(flist)
    d={}

    #store some global values for this dictionary
    d["dict_name"]=xname
    d["kllist"]=kllist
    d["weights"]=wts
    d["snrmeth"]=snrmeth
    d["snr_threshhold"]=snrthresh
    d["mode"]=mode
    d["datadir"]=datadir
    d["pestring"]=xstr

    #create directory to save in if doesn't yet exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #if highpass, collapse_mode, owa are not lists, repeat the entered value
    if not isinstance(hpval, list):
        hpval = np.repeat(hpval, npes)
    if not isinstance(collmode, list):
        collmode = np.repeat(collmode, npes)
    if not isinstance(owa, list):
        owa=np.repeat(owa, npes)
    if not isinstance(calflux, list):
        calflux=np.repeat(calflux,npes)

    #loop over number of PEs found 
    for i in np.arange(npes):
        # creates a dictionary with parameters needed to run KLIp for all images in a dataset
        #store PE file
        d["pe{0}".format(i+1)] = fits.getdata(pedir+flist[i])
        print('collapsing PE #', i+1, 'of', npes, ': ', flist[i])

        #store PE name
        pename = flist[i]
        d["pe{0}name".format(i+1)] = pename

        #extract IWA and store
        split=flist[i].split('_')
        s = re.findall(r'iwa(.*?)_', flist[i])

        d['pe{0}iwa'.format(i+1)]=int(s[0])

        head = fits.getheader(pedir+flist[i])
        file = fits.getdata(pedir+flist[i])
        nplanets = file.shape[3]
        #print(file.shape, nplanets)

        #if header keyword set, pull other KLIP values from header
        if header == True:
            hpval[i] = float(head["HIGHPASS"])*2
            collmode[i] = head["TIMECOLL"]
            owa[i] = head["OWA"]
            calflux[i] = [head["CALIBFLUX"]=="True"][0]
        
        #record in dictionary
        d["pe{0}hpval".format(i+1)]=hpval[i]
        d["pe{0}colmode".format(i+1)]=collmode[i]
        d["pe{0}owa".format(i+1)]=owa[i]
        d["pe{0}calflux".format(i+1)]=calflux[i]

        #separate out path, cut, and prefix
        #special case to match naming convention when more than one dataset in night
        #find the index in the list where the default naming string (_ax-yxz...)for the parameter explorer begins
        last=[val for val in split if val[0]=='a'][-1]
        lastind = split.index(last)
        dsetname = ''
        pfx = ''
        for ind in np.arange(1,lastind):
            dsetname+= ' '+str(split[ind])
            pfx+='_'+str(split[ind])

        #the only modification to obj_date folder names is the addition of long or short
        if 'long' in split or 'short' in split:
            d['pe{0}fpath'.format(i+1)]=datadir+split[1]+'/'+split[2]+'_'+split[3]+'/'
            d["pe{0}cut".format(i+1)]=split[5]
            #d["pe{0}pfx".format(i+1)]= split[1]+'_'+split[2]+'_'+split[3]+'_'+split[4]+'_'+split[5]
        else:
            d["pe{0}cut".format(i+1)]=split[4]
            d['pe{0}fpath'.format(i+1)]=datadir+split[1]+'/'+split[2]+'/'
            #d["pe{0}pfx".format(i+1)]= split[1]+'_'+split[2]+'_'+split[3]+'_'+split[4]
        
        #store dataset name for plot titles later
        #remove first underscore
        pfx=pfx[1:]
        d["pe{0}pfx".format(i+1)]=pfx
        d["pe{0}dset".format(i+1)] = dsetname
        writename = d["pe{0}pfx".format(i+1)]+xstr

        print(pename)

        ## extract best parameters according to weights and methods
        metric_cube, agg_cube, ann_val, movm_val, metric_scores, metric_fname = find_best_new(pename, kllist, pedir=pedir, outdir=outdir, \
            writestr=writename, weights=wts, snrmeth=snrmeth, smt=smt, snrthresh=snrthresh, separate_planets=separate_planets, separate_kls=separate_kls)

        d["pe{0}ann".format(i+1)]=ann_val
        d["pe{0}movm".format(i+1)]=movm_val
        d["pe{0}agg".format(i+1)]=agg_cube
        nkldim = ann_val.shape[0]
        npldim = ann_val.shape[1]
        strklip = [['']*npldim]*nkldim

        #save visualization of the metrics for this PE explorer
        if savefig==True:
            paramexplore_fig(pename, kllist, pedir=pedir, outdir=outdir, weights=wts, snrmeth=snrmeth, smt=smt)
    
        #define image input direcotries for KLIP based on PE filename
        haindir=d["pe{0}fpath".format(i+1)]+'dq_cuts/'+'Line_'+d["pe{0}cut".format(i+1)]+'_sliced/'
        contindir=d["pe{0}fpath".format(i+1)]+'dq_cuts/'+'Cont_'+d["pe{0}cut".format(i+1)]+'_sliced/'


        for kl in np.arange(nkldim):
            for pl in np.arange(npldim):        
                
                #set up some unique naming info
                prefix=d["pe{0}pfx".format(i+1)]+xstr
                iwa =d["pe{0}iwa".format(i+1)]
                strklip[kl][pl] = '_a'+str(ann_val[kl][pl])+'m'+str(movm_val[kl][pl])+'iwa'+str(iwa)
                #print(strklip[kl][pl])
                
                #mode loop. need to do everything twice if SDI. Otherwise, choose the correct dataset.
                if mode=='SDI':
                    modes=['Cont', 'Line']
                else:
                    runmode=mode
                    modes=['single']
                
                for j in np.arange(len(modes)):
                    #print(haindir, contindir)
                    if mode=='SDI':
                        runmode=modes[j]
                    if runmode=='Cont':
                        filelist = glob.glob(contindir + "sliced*.fits")
                    if runmode=='Line':
                        #default prefix has Cont in it, since we ran the PE on fake planets in cont images
                        prefix = prefix.replace('Cont','Line') 
                        filelist = glob.glob(haindir + "sliced*.fits")
                
                    #run KLIP with the optimum values for each PE

                    #check if file with this name already exists
                    if os.path.exists("{out}/{pre}-KLmodes-all.fits".format(out=outdir, pre=prefix+strklip[kl][pl])):
                        print("This file already exists. I am NOT re-running KLIP, but just reading the existing image in. Check and make sure you weren't intending to change the name")
                    #otherwise, run KLIP
                    else:
                        pdb.set_trace()
                        dataset = MagAO.MagAOData(filelist) 
                        dataset.IWA=iwa
                        dataset.OWA=float(owa[i])

                        print('running KLIP. highpass = ', hpval[i], ', calflux = ', calflux[i], ', time collapse = ', collmode[i], ', OWA = ', owa[i], 'prefix =', prefix, 'first file =', filelist[0])
                        parallelized.klip_dataset(dataset, outputdir=outdir, fileprefix=prefix+strklip[kl][pl], 
                                algo='klip', annuli=int(ann_val[kl,pl]), subsections=1, movement=movm_val[kl,pl],
                                numbasis=kllist, maxnumbasis=100, calibrate_flux=calflux[i], mode="ADI", highpass=float(hpval[i]), 
                                save_aligned=False, time_collapse=collmode[i])
                        #fits.writeto('test'+str(i)+'.fits', dataset.output, overwrite=True)

                    #pull output image 
                    klim = fits.getdata("{out}/{pre}-KLmodes-all.fits".format(out=outdir, pre=prefix+strklip[kl][pl]))
                    #print(klim.shape)
                    if kl==0 and pl==0:
                        klcube = np.zeros([nkldim, npldim, klim.shape[0],klim.shape[1],klim.shape[2]])
                    klcube[kl,pl,:,:,:]=klim

                    #and store in dictionary (overwrites at present until end of loop - needs work)
                    d["pe{0}strklip".format(i+1)]=strklip
                    if runmode=='Cont':
                        d["pe{0}contklipim".format(i+1)]=klcube
                    if runmode=='Line':
                        d["pe{0}haklipim".format(i+1)]=klcube
    return(d)

def paramexplore_fig(pename, kllist, pedir='proc/', outdir='proc/', writestr=False, weights=[1,1,1,1,1,1], 
    snrmeth='stdev', smt=3, separate_planets=False, separate_kls=False):
    
    metric_cube, agg_cube, ann_val, movm_val, metric_scores, metric_fname = find_best_new(pename, kllist, pedir=pedir, 
        outdir=outdir, writestr=writestr, weights=weights, snrmeth=snrmeth, smt=smt, separate_planets=separate_planets, 
        separate_kls=separate_kls)

    if writestr == False:
        writestr = pename[:-17]

    namelist = pename.split('_')
    params = [s for s in namelist if s[0] == 'a']
    params = params[0]
    params = re.split('a|-|x|m|s|iwa|hp', params)
    ymin = int(params[1])
    ymax = int(params[2])
    ystep = int(params[3])
    xmin = int(float(params[4]))
    xmax = int(float(params[5]))
    xstep = int(float(params[6]))
    nstepx = (xmax - xmin) / xstep
    nstepy = (ymax - ymin) / ystep

    fig_xdim = nstepx*0.25
    fig_ydim = nstepy

    fig = plt.figure(figsize=(fig_ydim,fig_xdim))
    gs = fig.add_gridspec(2, 4) # ,6)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    #ax5 = fig.add_subplot(gs[0,2])
    #ax6 = fig.add_subplot(gs[0,3])
    ax7 = fig.add_subplot(gs[0,2])
    ax8 = fig.add_subplot(gs[1,2])
    ax9 = fig.add_subplot(gs[:,3:])

    plt.setp((ax1, ax2, ax3, ax4, ax7, ax8, ax9), xticks=np.arange(0, nstepx + 1, 2), xticklabels=np.arange(xmin, xmax + 1, 2),
                 yticks=np.arange(0, nstepy + 1, 2), yticklabels=np.arange(ymin, ymax + 1, 2))

    #if extracted kl or planets separately, make figs separately
    nkldim = ann_val.shape[0]
    npldim = ann_val.shape[1]

    for kl in np.arange(nkldim):
        for pl in np.arange(npldim):

            im1 = ax1.imshow(metric_cube[0,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
            ax1.set_xlabel("movement parameter")
            ax1.set_ylabel("annuli parameter")
            ax1.set_title("Peak SNR: Weight = "+str(weights[0]))
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im1, cax=cax, orientation='vertical')

            im2 = ax2.imshow(metric_cube[1,0,kl,pl,:,:], origin='lower',cmap='magma', vmin=0, vmax=1)
            ax2.set_xlabel("movement parameter")
            ax2.set_ylabel("annuli parameter")
            ax2.set_title("Peak SNR Neighbor Quality: Weight = "+str(weights[1]))
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im2, cax=cax, orientation='vertical')

            im3 = ax3.imshow(metric_cube[2,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
            ax3.set_xlabel("movement parameter")
            ax3.set_ylabel("annuli parameter")
            ax3.set_title("Avg SNR Under Mask: Weight = "+str(weights[2]))
            divider = make_axes_locatable(ax3)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im3, cax=cax, orientation='vertical')

            im4 = ax4.imshow(metric_cube[3,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
            ax4.set_xlabel("movement parameter")
            ax4.set_ylabel("annuli parameter")
            ax4.set_title("Avg SNR Neighbor Quality: Weight = "+str(weights[3]))
            divider = make_axes_locatable(ax4)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im4, cax=cax, orientation='vertical')

            
            #im5 = ax5.imshow(metric_cube[4,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
            #ax5.set_xlabel("movement parameter")
            #ax5.set_ylabel("annuli parameter")
            #ax5.set_title("Stdev Across KL: Weight = "+str(weights[4]))
            #divider = make_axes_locatable(ax5)
            #cax = divider.append_axes('right', size='5%', pad=0.05)
            #plt.colorbar(im5, cax=cax, orientation='vertical')

            #im6 = ax6.imshow(metric_cube[5,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
            #ax6.set_xlabel("movement parameter")
            #ax6.set_ylabel("annuli parameter")
            #ax6.set_title("Stdev Neighbor Quality: Weight = "+str(weights[5]))
            #divider = make_axes_locatable(ax6)
            #cax = divider.append_axes('right', size='5%', pad=0.05)
            #plt.colorbar(im6, cax=cax, orientation='vertical')

            im7 = ax7.imshow(metric_cube[4,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
            ax7.set_xlabel("movement parameter")
            ax7.set_ylabel("annuli parameter")
            ax7.set_title("Spurious Pixels: Weight = "+str(weights[4]))
            divider = make_axes_locatable(ax7)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im7, cax=cax, orientation='vertical')


            im8 = ax8.imshow(metric_cube[5,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
            ax8.set_xlabel("movement parameter")
            ax8.set_ylabel("annuli parameter")
            ax8.set_title("Contrast: Weight = "+str(weights[5]))
            divider = make_axes_locatable(ax8)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im8, cax=cax, orientation='vertical')

            # plot metric
            im9 = ax9.imshow(agg_cube[kl,pl,:,:], origin='lower', vmin=0, vmax=np.sum(weights))
            ax9.set_ylabel("annuli parameter")
            ax9.set_xlabel("movement parameter")
            ax9.set_title("Aggregate Parameter Quality")
            divider = make_axes_locatable(ax9)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im9, cax=cax, orientation='vertical')

            ind = np.where(agg_cube[kl,pl,:,:] == np.nanmax(agg_cube[kl,pl,:,:]))
            label_text = 'a' + str(ann_val[kl][pl]) + 'm' + str(movm_val[kl][pl])
            rect = patches.Rectangle((ind[1][0] - 0.5, ind[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='r', facecolor='none')
            ax9.add_patch(rect)
            ax9.text(ind[1][0] + 0.75, ind[0][0], label_text, color='red')

            plt.suptitle(writestr)
            gs.tight_layout(fig, rect=[0, 0.03, 1, 0.95])

            plt.savefig(outdir+writestr+'_kl'+str(kllist[kl])+'_pl'+str(pl)+'_paramqual.png')
    
    return

def make_klip_snrmaps(d, pedir='./', outdir='proc/', smooth=0.5, snrmeth='stdev', mode='Line', scale=1):
    """
    Generate SNR maps with specified parameters for all KLIP images in a directory

    REQUIRED INPUTS:
    d:      The dictionary storing the KLIPed images

    OPTIONAL INPUTS:
    outdir:     location to store output files
    smooth:     gaussian FWHM with which to smooth the map
    mode:       what type of image to pull and generate SNR maps for
    snrmeth:    one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value),
                and 'all' to average the two
    scale:      for SDI mode, how to scale the images

    RETURNS:
    d:      Ths input dictionary with SNR map images added

    """

    d["smooth"]=smooth

    #figure out how many pes we have
    flist = glob.glob(pedir+'paramexplore*KLmodes-all.fits')
    flist = sorted(flist)
    npes=len(flist)

    #loop over the pes
    for i in np.arange(npes):

        #pull various naming strings
        strklip=d["pe{0}strklip".format(i+1)]
        prefix=d["pe{0}pfx".format(i+1)]
        
        #pull fwhm - needed for Mawet correction
        head=fits.getheader(pedir+d["pe{0}name".format(i+1)])
        fwhm = float(head["SNRFWHM"])
        d["pe{0}fwhm".format(i+1)]=fwhm
        
        #pull the correct image(s) to create the map
        if mode=='Line' or mode=='SDI':
            haim = d["pe{0}haklipim".format(i+1)]
            klipim = haim
            writeprefix = prefix.replace('Cont','Line')
        if mode=='Cont' or mode=='SDI':
            contim = d["pe{0}contklipim".format(i+1)]
            klipim = contim
            writeprefix=prefix

        #make SDI image
        if mode=='SDI':
            klipim = haim - scale*contim
            writeprefix = prefix.replace('Cont','SDI')
            
        #name and make SNR Map
        outname = outdir+writeprefix+strklip + '_sm'+str(smooth)+'_'+mode+'SNRMap_'+snrmeth+'.fits'

        #check whether file already exists
        if os.path.exists(outname):
            print("This file already exists. I am NOT re-running SNRMap, but just reading the existing image in. Check and make sure you weren't intending to change the name")
            Output = fits.getdata(outname)
        else:   
            Output = snr.create_map(klipim, fwhm, smooth=smooth, saveOutput=True, outputName=outname[:-5], checkmask=False, method=snrmeth)
        
        #store maps
        if mode=='Line':
            d["pe{0}hasnrmap".format(i+1)]=Output
        if mode=='Cont':
            d["pe{0}contsnrmap".format(i+1)]=Output
        if mode=='SDI':
            d["pe{0}sdisnrmap".format(i+1)]=Output
            
    return(d)

def compare_pes(d, pedir='./', outdir='klipims/', nrows=False, ncols=False, save=None, title=None):
    
    """
    Generates an image with a grid of collapsed PE aggregate metric score maps for all PEs in a dictionary.

    REQUIRED INPUTS
    d:      The dictionary storing the PE info

    """

    #how many PEs are there?
    flist = glob.glob(pedir+'paramexplore*KLmodes-all.fits')
    npes=len(flist)
    
    #make string lists of relevant dictionary keys
    imlist = ["pe{0}agg".format(i+1) for i in np.arange(npes)]
    annlist = ["pe{0}ann".format(i+1) for i in np.arange(npes)]
    movmlist = ["pe{0}movm".format(i+1) for i in np.arange(npes)]
    namelist = ["pe{0}dset".format(i+1) for i in np.arange(npes)]

    #if grid rows and columns not specified, set 3 columns and enough rows to fit
    if nrows == False:
        ncols=3
        nrows = int(np.ceil(npes/ncols))
    
    #fix for indexing error when don't give it enough spots
    if npes>ncols*nrows:
        print('Not enough spots for number of pes. Increase ncols or nrows.')

    #size figure according to how many rows and columns there are
    figsz=(ncols*4, nrows*5)

    #set up plot
    f, ax = plt.subplots(nrows, ncols, figsize=figsz)
    ax = ax.ravel()

    #add master title for the grid
    if title != None:
        f.suptitle(title, size=20)
    
    #loop over all of the images
    for i in np.arange(len(imlist)):

        #pull the aggregate PE metric map
        im = ax[i].imshow(d[imlist[i]], origin='lower', vmin=1, vmax=np.nanmax(d[imlist[i]]))
        ax[i].set_ylabel("annuli parameter")
        ax[i].set_xlabel("movement parameter")

        #add colorbar
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        plt.colorbar(im, cax=cax, orientation='vertical', label="Parameter Quality Metric")
        
        #find and mark the peak location
        ind = np.where(d[imlist[i]] == np.nanmax(d[imlist[i]]))
        label_text = 'a' + str(d[annlist[i]]) + 'm' + str(d[movmlist[i]])
        rect = patches.Rectangle((ind[1][0] - 0.5, ind[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='r', facecolor='none')
        ax[i].add_patch(rect)
        ax[i].text(ind[1][0] + 0.75, ind[0][0], label_text, color='red')

        #label subplot with dataset title
        ax[i].set_title("\n".join(textwrap.wrap(d[namelist[i]], 30)))

    plt.tight_layout()
    if save != None:
        plt.savefig(outdir+save+d["pestring"])
    return()


def compare_ims(d, pedir='./', kllist=[5,10,20,50], outdir='klipims/', mode='Line', nrows=False, ncols=False, save=None, title=None, boxsz=50):
    """
    Generates an image with a grid of collapsed PE aggregate metric score maps for all PEs in a dictionary.

    REQUIRED INPUTS
    d:      The dictionary storing the PE info

    OPTIONAL INPUTS:
    boxsz:  radius of box around image center to show

    """

    flist = glob.glob(pedir+'paramexplore*KLmodes-all.fits')
    npes=len(flist)

    #same as compare_pes, but generates grid of snrmaps
    if mode=='SDI':
        snmaplist = ["pe{0}sdisnrmap".format(i+1) for i in np.arange(npes)]
    if mode=='Line':
        snmaplist = ["pe{0}hasnrmap".format(i+1) for i in np.arange(npes)]
    if mode=='Cont':
        snmaplist = ["pe{0}contsnrmap".format(i+1) for i in np.arange(npes)]
    
    namelist = ["pe{0}dset".format(i+1) for i in np.arange(npes)]

    #if grid rows and columns not specified, set 3 columns and enough rows to fit
    if nrows == False:
        ncols=3
        nrows = int(np.ceil(npes/ncols))

        #fix for indexing error when don't give it enough spots
    if npes>ncols*nrows:
        print('Not enough spots for number of pes. Increase ncols or nrows.')

    #loop over KL modes
    for k in kllist:

        #size figure according to how many rows and columns there are
        figsz=(ncols*4, nrows*5)
    
        #set up plot
        f, ax = plt.subplots(nrows, ncols, figsize=figsz)
        ax = ax.ravel()

        #add master title for the grid
        if title != None:
            f.suptitle(title + ' - KL ' + str(k), size=20)

        #loop over all of the images   
        for i in np.arange(len(snmaplist)):

            klind = [i for i in range(len(kllist)) if kllist[i] == k]
            print(k, 'kl modes is index', klind)

            #figure out where the image center is
            dims = d[snmaplist[i]][0,0,:,:].shape
            cen = int((dims[0]-1)/2.)

            #pull the region around the image center for plotting
            im = ax[i].imshow(d[snmaplist[i]][0,klind[0],cen-boxsz:cen+boxsz,cen-boxsz:cen+boxsz], cmap='magma',
                              origin='lower', vmin=-2, vmax=5)
            ax[i].set_ylabel("")
            ax[i].set_xlabel("")

            #add colorbar
            divider = make_axes_locatable(ax[i])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im, cax=cax, orientation='vertical', label="SNR")
            
            #label subplot with dataset title
            ax[i].set_title("\n".join(textwrap.wrap(d[namelist[i]], 30)))

        plt.tight_layout()
        if save != None:
            plt.savefig(outdir+save+'_kl'+str(k)+d["pestring"])
    return()


def save_pe_dict(d, dwritename, doutdir='./dicts/'):
    if not os.path.exists(doutdir):
        os.makedirs(outdir)
    pickle.dump(d,open(doutdir+dwritename+".p","wb"))


def pull_dsets(false_dir, real_dir, out_dir, namestr='*',exceptstr = False, dateorder=True):
  """
  Pull all the real planet H-alpha parameter explorers and look for continuum false planet matches.

  Required Input:
  Directories for false planets, real planets, and output

  Optional Inputs:
  namestr = a string (e.g. objname) that all paramexplores you want compiled have

  Outputs:
  list of dataset names (objname_date)
  list of real planet pe filenames
  list of false planet pe filenames
  """
  thisdir = os.getcwd()
  os.chdir(real_dir)
  reallist = glob.glob('paramexplore*'+namestr+'*.fits')
  os.chdir(false_dir)
  falselist = glob.glob('paramexplore*'+namestr+'*.fits')
  os.chdir(thisdir)

  #put reallist in date order, as this will be the processing/display order
  if dateorder==True:
    monthstrs = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    t=[]
    for i in np.arange(len(reallist)):
        #print(reallist[i])
        linefname = reallist[i].split('_')
        line_date = linefname[2]
        monthstr = [s for s in monthstrs if s in line_date]
        monthno = [i for i in np.arange(1,13) if monthstrs[i-1] in line_date]
        year = int('20'+line_date[-2:])
        day = line_date.split(monthstr[0])[0]
        #print(line_date, '= day:', day, ' month:', monthno[0], ' year:', year)
        t.append(date(year=year, month=int(monthno[0]), day=int(day)))

    #now sort
    reallist = [x for _, x in sorted(zip(t, reallist))]

  dset_stringlist=[]  
  rlist = []
  flist = []
  #extract info and find matches 
  for i in np.arange(len(reallist)):
    linefname = reallist[i].split('_')
    #extract object name, date, wl, cut
    line_obj = linefname[1]
    line_date = linefname[2]
    line_wl = linefname[3]
    line_cut = linefname[4]
    line_klip_params = linefname[-2]
    line_hpval = linefname[-1].split('.')[0]

    match=False
    #look for matching continuum
    for j in np.arange(len(falselist)):
      contfname = falselist[j].split('_')
      cont_obj = contfname[1]
      cont_date = contfname[2]
      cont_wl = contfname[3]
      cont_cut = contfname[4]
      cont_klip_params = contfname[-2]
      cont_hpval = contfname[-1].split('.')[0]
      if line_obj == cont_obj:
        if line_date == cont_date:
          if (line_wl == 'Line') and (cont_wl=='Cont'):
            if line_cut == cont_cut:
              if line_klip_params == cont_klip_params:
                if line_hpval == cont_hpval:
                  match=True
                  falseind=j
                  dset_str = line_obj+'_'+line_date
                  if str(exceptstr) not in dset_str:
                      dset_stringlist.append(dset_str)
                      rlist.append(reallist[i])
                      flist.append(falselist[j])

  return(dset_stringlist, rlist, flist)


def proc_one_dset(dset_string, realpe,falsepe, false_dir, real_dir, out_dir, pklstr='*', overwrite=False):
  """
  process one dataset at a time, generating all possible klcombo, parametercombo pairs and save normalized difference
  stats and images in a dictionary

  INPUTS:
  dset_string: the unique OBJNAME_DATE string for the dataset
  realpe:  the filename of the real planet H-alpha pe matching that string
  falsepe: the filename of the false planet Continuum pe matching both 
  """
  outfname = out_dir+dset_string+pklstr+'.p'
  if os.path.exists(outfname):
    print('pickle file' + outfname + 'already exists.')
    if overwrite==True:
        print('but overwrite is true. Rerunning.')
    else:
        return()

  print('processing', dset_string)

  sepkl=True
  snt=False
  
  i=0
  skip=0
  skip2=0
  dtn={}
  dtn2={}
  
  kllist = [1,2,3,4,5,10,20,50,100]
  #how many unique combos for each object?
  klpossible = it.product((np.nan,1), repeat=9)

  bar = progressbar.ProgressBar(maxval=52, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
  barctr=0
  bar.start()

  for klcombo in klpossible:
    if np.nansum(klcombo)>=1:
      j=0
      klchoice = np.multiply(kllist,list(klcombo))
      kls = [int(kl) for kl in klchoice if (np.isnan(kl)==False)]
      
      metricpossible = it.product((0,1), repeat=6)
      
      for mcombo in metricpossible:

        if np.nansum(mcombo)>=1:
          go=True
          
          #will break on false pos only (degenerate)
          if (np.sum(mcombo)==1) and (mcombo[4]==1):
            skip+=1
            go = False

          #stdev across kl modes has no meaning if nkl !>1 skip all metrics
          #elif (np.nansum(klcombo)<3) and (np.sum(mcombo[4:6])>=1):
            #skip+=1
            #go = False

          else:

            if go==True:

              try:
                fake_metric_cube, fake_agg_cube, fake_ann_val, fake_movm_val, fake_metric_scores, fake_metric_fname = find_best_new(falsepe, kls, pedir=false_dir, writestr=False, writefiles=False, weights=list(mcombo), outdir=false_dir+'proc/', 
                                                                                                                                            oldpe=False, debug=False, smt=3, snrmeth='stdev',separate_planets=False, separate_kls=sepkl, snrthresh=snt)
                real_metric_cube, real_agg_cube, real_ann_val, real_movm_val, real_metric_scores, real_metric_fname = find_best_new(realpe, kls, pedir=real_dir, writestr=False, writefiles=False, weights=list(mcombo), outdir=real_dir+'proc/', 
                                                                                                                                            oldpe=False, debug=False, smt=3, snrmeth='stdev',separate_planets=False, separate_kls=sepkl, snrthresh=snt) 
                go=True  

              except:
                go=False
                skip2+=1 

            if go==True:
              fake_cube = fake_metric_cube
              real_cube = real_metric_cube

              ##NORMALIZE TO 90th pctile AND SUBTRACT 10th pctile for each KL mode
              sig_diff=np.zeros((len(kls)))*np.nan
              sig_diff_wt=np.zeros((len(kls)))*np.nan

              for kl in np.arange(len(kls)):
                fake_cube[-1,:,kl,:,:,:]=fake_cube[-1,:,kl,:,:,:]-np.nanpercentile(fake_cube[-1,:,kl,:,:,:],10)
                real_cube[-1,:,kl,:,:,:]=real_cube[-1,:,kl,:,:,:]-np.nanpercentile(real_cube[-1,:,kl,:,:,:],10)
                        
                #normalize to max
                fake_cube[-1,:,kl,:,:,:]=fake_cube[-1,:,kl,:,:,:]/np.nanpercentile(fake_cube[-1,:,kl,:,:,:],90)
                real_cube[-1,:,kl,:,:,:]=real_cube[-1,:,kl,:,:,:]/np.nanpercentile(real_cube[-1,:,kl,:,:,:],90)
                                  
                diff = fake_cube-real_cube
                #weight it by the normalized aggregate metric (downweights low SNR regions)
                diff_wt = diff*fake_cube[-1,:,kl,:,:,:]

                #differences in metric score for real planet at fake planet peak
                sig_diff[kl] = diff[-1,0,kl,0,int(fake_ann_val[kl][0]-1),int(fake_movm_val[kl][0])]
                sig_diff_wt[kl] = diff_wt[-1,0,kl,0,int(fake_ann_val[kl][0]-1),int(fake_movm_val[kl][0])]
                
              #kl x ann x movm grids of aggregate param qual metric 
              diff_wt_agg = diff_wt[-1,0,:,0,:,:]  
              diff_agg = diff[-1,0,:,0,:,:]

              #avg metric score diffs across KL mode
              sig_diff_avg = np.mean(sig_diff)
              sig_diff_wt_avg = np.mean(sig_diff_wt)

              #difference in metric score at fake planet peak
              coll_string = dset_string+'_wts'+''.join([str(x) for x in mcombo])+'_kls'+''.join([str(x) if str(x)=='1' else '0' for x in list(klcombo)])+'_sepkl'+str(sepkl)+'_snthresh'+str(snt)
              dtn[coll_string]=(np.nansum(diff_wt_agg), np.nanstd(diff_wt_agg), np.nanmedian(diff_wt_agg), sig_diff_wt_avg)
              dtn2[coll_string]=(np.nansum(diff_agg), np.nanstd(diff_agg), np.nanmedian(diff_agg), sig_diff_avg) 
              
              j+=1
        else:
          #print('skipping zero sum!', mcombo)
          continue

      if i%10==0:
        barctr+=1
        bar.update(barctr)
        sleep(0.5)
      i+=1
      #print('klcombo', i, '/511 done')
    else:
      print('skipped kl')
      continue

  bar.finish()
  print('skipped', skip)
  print('no unique soln', skip2)
  print("dumping dict for dataset", dset_string, 'with', i, 'klcombos and', j, 'metric combos')
  pickle.dump(dtn,open(out_dir+dset_string+pklstr+'wt.p','wb'))
  pickle.dump(dtn2,open(out_dir+dset_string+pklstr+'.p','wb'))
  
  return(dtn,dtn2)

def split_by_completeness(false_dir, real_dir, out_dir, namestr='*', pklstr='*', exceptstr=False):

    dset_stringlist, rlist, flist = pull_dsets(false_dir, real_dir, out_dir, namestr=namestr, exceptstr=exceptstr)

    ##process only the ones that have not yet been run
    i=0
    new_dsetlist=[]
    new_rlist=[]
    new_flist=[]
    dset_done=[]
    rlist_done=[]
    flist_done=[]
    for dset in dset_stringlist:
        dictlist = glob.glob(out_dir+dset+pklstr+'.p')
        if len(dictlist)>0:
            print('match exists already for', dset, ', namely:', dictlist)
            dset_done.append(dset_stringlist[i])
            rlist_done.append(rlist[i])
            flist_done.append(flist[i])
        else:
            new_dsetlist.append(dset_stringlist[i])
            new_rlist.append(rlist[i])
            new_flist.append(flist[i])
        i+=1

    toproc = (new_dsetlist, new_rlist, new_flist)
    done = (dset_done, rlist_done, flist_done)

    return(toproc, done)

def batch_dset_proc(dset_stringlist, toproc, false_dir, real_dir, out_dir,  pklstr='*', overwrite=False, parallel=False):

    dsetlist, rlist, flist = toproc

    if parallel==True:
        threadlist = []
        for k in np.arange(len(dsetlist)):
            t = threading.Thread(target=proc_one_dset,args=(dsetlist[k],rlist[k],flist[k],false_dir, real_dir, out_dir,  pklstr, overwrite))
            threadlist.append(t)
            t.start()
        
        for tr in threadlist:
            tr.join()
            print("Finished")
    else:
        for i in np.arange(len(dset_stringlist)):
            print(i)
            proc_one_dset(dsetlist[i],rlist[i],flist[i], false_dir, real_dir, out_dir,  pklstr=pklstr, overwrite=overwrite)
    return()

def compile_keys(out_dir, done, pklstr):

    dset_done, rlist_done, flist_done = done

    #sanity check - number of dictionary keys is consistent with the number of unique combos 
    i=0
    
    #how many unique combos for each object?
    klpossible = it.product((np.nan,1), repeat=9)
    skip=0
    k=0
    
    for klcombo in klpossible:
      if np.nansum(klcombo)>=1:
        j=0
        metricpossible = it.product((0,1), repeat=6)
        for mcombo in metricpossible:
          if np.sum(mcombo)>=1:
            test = list(mcombo)
            #will break on false pos only (degenerate)
            if (np.sum(mcombo)==1) and (test[4]==1):
              skip+=1
            #stdev has no meaning if nkl !>2 skip all metrics weighting this
            #elif (np.nansum(klcombo)<2) and (np.sum(test[4:6])>=1):
              #skip+=1         
            else:
              if np.sum(mcombo)==0:
                print('WTF')
              i+=1
              j+=1
          else:
            continue
        k+=1
    print('skipped', skip)
    print('total possible combos:', i)
    print('total metric combos:', j)
    print('total kl combos:', k)

    #now make a list of all keys in all dicts
    dset_nkeys =[]
    current_keys = []
    
    for dictname in dset_done:
        print(dictname)
        d = pickle.load( open( out_dir+dictname+pklstr+'.p', "rb" ) )
      
        keys = d.keys()
        keys_compiled = []

        for key in keys:
            keys_compiled.append(key)
        
        dset_nkeys.append(len(keys_compiled))

        if len(keys_compiled)>=len(current_keys):
            current_keys = keys_compiled
        
        if len(np.unique(current_keys)) != len(current_keys):
            print('why arent the keys unique?')
        else:
            print('dataset', dictname, 'has', len(keys_compiled), 'unique keys')
      
    print('total keys:', np.sum(dset_nkeys))
    print('max keys per dataset:', len(current_keys))
    return(current_keys)

def compute_bulk_diagnostics(current_keys, pklstr, done, out_dir):

    dset_done, rlist_done, flist_done = done

    #compare across combos
    #set up empty arrays to fill
    #dimensions are: 0) keys per dataset, 1) datasets, number of metrics, number of kls
    #for collapsing by kl or metric
    wt_sum_matrix = np.zeros((len(current_keys)*len(dset_done),6,9))*np.nan
    sum_matrix = np.zeros((len(current_keys)*len(dset_done),6,9))*np.nan
    wt_std_matrix = np.zeros((len(current_keys)*len(dset_done),6,9))*np.nan
    std_matrix = np.zeros((len(current_keys)*len(dset_done),6,9))*np.nan
    wt_med_matrix = np.zeros((len(current_keys)*len(dset_done),6,9))*np.nan
    med_matrix = np.zeros((len(current_keys)*len(dset_done),6,9))*np.nan
    wt_mdiff_matrix = np.zeros((len(current_keys)*len(dset_done),6,9))*np.nan
    mdiff_matrix = np.zeros((len(current_keys)*len(dset_done),6,9))*np.nan


    p=0 #dataset counter
    k=0 #key counter
    #loop through 
    refkeys=[]
    for pname in dset_done:
      d = pickle.load( open( out_dir + pname + pklstr + 'wt.p', "rb" ) )
      d2 = pickle.load( open( out_dir + pname + pklstr + '.p', "rb" ) )
      keys = d.keys()

      ##loop over possible klcombo, mcombo pairs
      skip=0
      for key in keys:
        key_split  = key.split('_')
        #pull out parameter weights
        key_wts = key_split[2].split('wts')
        wt_list = [int(x) if int(x)==1 else np.nan for x in key_wts[1]]

        #pull out kllist
        key_kls = key_split[3].split('kls')
        kl_list = [int(x) if int(x)==1 else np.nan for x in key_kls[1]]

        #fill in sum, stdev and mean values for difference maps (weighted and unweighted)       
        #np.outer of the wt_list and kl_list is a nmetric x nkl grid of nans and 1s
        #so sum/std/mean will be repeated for every combo of metric/kl for which it's "on" 
        #add only if the weights and kl modes match
        
        #check same
        wt_sum_matrix[k,:,:]=np.outer(wt_list,kl_list)*(d[key][0])
        sum_matrix[k,:,:]=np.outer(wt_list,kl_list)*(d2[key][0])
        wt_std_matrix[k,:,:]=np.outer(wt_list,kl_list)*(d[key][1])
        std_matrix[k,:,:]=np.outer(wt_list,kl_list)*(d2[key][1])
        wt_med_matrix[k,:,:]=np.outer(wt_list,kl_list)*(d[key][2])
        med_matrix[k,:,:]=np.outer(wt_list,kl_list)*(d2[key][2])
        wt_mdiff_matrix[k,:,:]=np.outer(wt_list,kl_list)*(d[key][-1])
        mdiff_matrix[k,:,:]=np.outer(wt_list,kl_list)*(d2[key][-1])
        k+=1
      p+=1

      wt_diagnostics = (wt_sum_matrix, wt_std_matrix, wt_med_matrix, wt_mdiff_matrix)
      diagnostics = (sum_matrix, std_matrix, med_matrix, mdiff_matrix)
      print("skipped",skip)

    return(wt_diagnostics,diagnostics)


def diaghists_bykl(diagnostics):

    sum_matrix, std_matrix, med_matrix, mdiff_matrix = diagnostics

    ##plot histograms for the different kl modes individually
    kls = [1,2,3,4,5,10,20,50,100]
    for i in np.arange(9):
        f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(16.,3))
        #np.ravel to tally in 1D
        thissum = np.ravel(sum_matrix[:,:,i])
        ax1.hist(thissum,bins=50)
        ax1.set_xlim(-100,500)
        ax1.set_title('Sum of weighted differences \n KL = '+ str(kls[i]) + ' mean:{:.2f} std:{:.2f}'.format(np.nanmean(thissum),np.nanstd(thissum)) )
        thisstd = np.ravel(std_matrix[:,:,i])
        ax2.hist(thisstd,bins=100)
        ax2.set_xlim(0,1)
        ax2.set_title('Stdev of weighted differences \n KL = '+ str(kls[i]) + ' mean:{:.2f} std:{:.2f}'.format(np.nanmean(thisstd),np.nanstd(thisstd)) )
        thismed = np.ravel(med_matrix[:,:,i])
        ax3.hist(thismed,bins=50)
        ax3.set_xlim(-0.1,0.3)
        ax3.set_title('Median of weighted differences \n KL = '+ str(kls[i]) + ' mean:{:.2f} std:{:.2f}'.format(np.nanmean(thismed),np.nanstd(thismed)) )
        
        thismdiff = np.ravel(mdiff_matrix[:,:,i])
        ax4.hist(thismdiff,bins=50)
        ax4.set_xlim(-0.1,1)
        ax4.set_title('Diff. of metrics at cont peak \n KL = '+ str(kls[i]) + ' mean:{:.2f} std:{:.2f}'.format(np.nanmean(thismdiff),np.nanstd(thismdiff)) )
    plt.show()

def diaghists_bymetric(diagnostics):

    sum_matrix, std_matrix, med_matrix, mdiff_matrix = diagnostics

    ##plot histograms for the different metrics individually
    ms = ['Peak SNR','NQ peak SNR','Avg SNR','NQ avg SNR','stdev across KL','NQ stdev', 'spur pix', 'contrast']

    for i in np.arange(6):
        f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(16.,3))
        #np.ravel to tally in 1D
        thissum = np.ravel(sum_matrix[:,i,:])
        ax1.hist(thissum,bins=50)
        ax1.set_xlim(-100,500)
        ax1.set_title('Sum of weighted differences \n '+ str(ms[i]) + ' mean:{:.2f} std:{:.2f}'.format(np.nanmean(thissum),np.nanstd(thissum)) )
        thisstd = np.ravel(std_matrix[:,i,:])
        ax2.hist(thisstd,bins=100)
        ax2.set_xlim(0,1)
        ax2.set_title('Stdev of weighted differences \n '+ str(ms[i]) + ' mean:{:.2f} std:{:.2f}'.format(np.nanmean(thisstd),np.nanstd(thisstd)) )
        thismed = np.ravel(med_matrix[:,i,:])
        ax3.hist(thismed,bins=50)
        ax3.set_xlim(-0.1,0.3)
        ax3.set_title('Median of weighted differences \n '+ str(ms[i]) + ' mean:{:.2f} std:{:.2f}'.format(np.nanmean(thismed),np.nanstd(thismed)) )
        thismdiff = np.ravel(mdiff_matrix[:,i,:])
        ax4.hist(thismdiff,bins=50)
        ax4.set_xlim(-0.1,1)
        ax4.set_title('Diff. of metrics at cont peak \n '+ str(ms[i]) + ' mean:{:.2f} std:{:.2f}'.format(np.nanmean(thismdiff),np.nanstd(thismdiff)) )
        plt.show()

def compute_dset_diagnostics(current_keys, done, pklstr, out_dir):
    
    dset_done, rlist_done, flist_done = done

    #set up empty arrays to fill
    #dimensions are: 0) keys per dataset, 1) datasets

    #for analyzing by dataset
    wt_sum_list = np.zeros((len(current_keys),len(dset_done)))*np.nan
    sum_list = np.zeros((len(current_keys),len(dset_done)))*np.nan
    wt_std_list = np.zeros((len(current_keys),len(dset_done)))*np.nan
    std_list = np.zeros((len(current_keys),len(dset_done)))*np.nan
    wt_med_list = np.zeros((len(current_keys),len(dset_done)))*np.nan
    med_list = np.zeros((len(current_keys),len(dset_done)))*np.nan

    wt_mdiff_list = np.zeros((len(current_keys),len(dset_done)))*np.nan
    mdiff_list = np.zeros((len(current_keys),len(dset_done)))*np.nan

    p=0 #dataset counter
    #loop through 
    refkeys=[]
    for pname in dset_done:
        d = pickle.load( open( out_dir + pname + pklstr + 'wt.p', "rb" ) )
        d2 = pickle.load( open( out_dir + pname + pklstr + '.p', "rb" ) )
        keys = d.keys()
        k=0 #reset key counter

        ##loop over possible klcombo, mcombo pairs
        skip=0
        for key in keys:
            key_split  = key.split('_')
            #pull out parameter weights
            key_wts = key_split[2].split('wts')
            wt_list = [int(x) if int(x)==1 else np.nan for x in key_wts[1]]

            #pull out kllist
            key_kls = key_split[3].split('kls')
            kl_list = [int(x) if int(x)==1 else np.nan for x in key_kls[1]]
            
            uniqname = str(key_wts[1])+str(key_kls[1])
            if p==0:
                refkeys.append(uniqname)

            #fill in sum, stdev and mean values for difference maps (weighted and unweighted)       
            #np.outer of the wt_list and kl_list is a nmetric x nkl grid of nans and 1s
            #so sum/std/mean will be repeated for every combo of metric/kl for which it's "on" 
            #add only if the weights and kl modes match
            
            #check same
            if uniqname==refkeys[k]:
                #also in list form
                wt_sum_list[k,p]=(d[key][0])
                sum_list[k,p]=(d2[key][0])
                wt_std_list[k,p]=(d[key][1])
                std_list[k,p]=(d2[key][1])
                wt_med_list[k,p]=(d[key][2])
                med_list[k,p]=(d2[key][2])
                wt_mdiff_list[k,p]=(d[key][-1])
                mdiff_list[k,p]=(d2[key][-1])
        
            #if not same, find the right place to put it
            else:
                if uniqname in refkeys:
                    refind = [i for i in range(len(refkeys)) if refkeys[i] == uniqname]
                    wt_sum_list[refind[0],p]=(d[key][0])
                    sum_list[refind[0],p]=(d2[key][0])
                    wt_std_list[refind[0],p]=(d[key][1])
                    std_list[refind[0],p]=(d2[key][1])
                    wt_med_list[refind[0],p]=(d[key][2])
                    med_list[refind[0],p]=(d2[key][2])  
                    wt_mdiff_list[refind[0],p]=(d[key][-1])
                    mdiff_list[refind[0],p]=(d2[key][-1])       
                else:
                    skip+=1

            k+=1
        p+=1

    wt_diaglists = (wt_sum_list, wt_std_list, wt_med_list, wt_mdiff_list)
    diaglists = (sum_list, std_list, med_list, mdiff_list)
    print("skipped", skip)
    return(wt_diaglists, diaglists, refkeys)

def avg_diag_across_dsets(diaglists, refkeys):

    sum_list, std_list, med_list, mdiff_list = diaglists
    #loop over all possible combos
    mean_sum = []
    mean_std = []
    mean_med = []
    mean_mdiff = []
    for i in np.arange(len(refkeys)):
      mean_sum.append(np.mean(sum_list[i,:]))
      mean_std.append(np.mean(std_list[i,:]))
      mean_med.append(np.mean(med_list[i,:]))
      mean_mdiff.append(np.mean(mdiff_list[i,:]))

    masterlists = (mean_sum, mean_std, mean_med, mean_mdiff)
    return(masterlists)

def find_best_diag(masterlists, cutoffpct, refkeys):

    mean_sum, mean_std, mean_med, mean_mdiff = masterlists

    sumcut = np.nanpercentile(np.abs(mean_sum),cutoffpct)
    stdcut = np.nanpercentile(mean_std,cutoffpct)
    medcut = np.nanpercentile(np.abs(mean_med),cutoffpct)
    diffcut = np.nanpercentile(np.abs(mean_mdiff), cutoffpct)

    cutoffs = (sumcut, stdcut, medcut, diffcut)
    bestlist=[]
    for i in np.arange(len(refkeys)):
      if np.abs(mean_sum[i]) < sumcut:
        if mean_std[i] < stdcut:
          if np.abs(mean_med[i]) < medcut:
            if np.abs(mean_mdiff[i]) < diffcut:
              print(refkeys[i][:6], refkeys[i][6:], mean_sum[i],mean_std[i],mean_med[i],mean_mdiff[i])
              bestlist.append((refkeys[i][:6], refkeys[i][6:],mean_sum[i],mean_std[i],mean_med[i],mean_mdiff[i]))

    return(bestlist, cutoffs)

def diaghist_wcutoff(current_keys, done, pklstr, out_dir, cutoffpct, wt=True):

    wt_diaglists, diaglists, refkeys = compute_dset_diagnostics(current_keys, done, pklstr, out_dir)
    if wt==True:
        diaglists = wt_diaglists
    masterlists = avg_diag_across_dsets(diaglists, refkeys)
    bestlist, cutoffs = find_best_diag(masterlists, cutoffpct, refkeys)

    mean_sums, mean_stds, mean_meds, mean_mdiffs = masterlists
    sumcut, stdcut, medcut, diffcut = cutoffs

    wtslist = []
    klslist = []

    #loop through all metric combos meeting criteria
    for i in np.arange(len(bestlist)):
        
        wtstr, klstr, mean_sum, mean_std, mean_med, mean_mdiff = bestlist[i]

        wts = [int(i) for i in wtstr]
        wtslist.append(wts)
        kl_binaries = [int(k) for k in klstr]
        kllist = np.array([1,2,3,4,5,10,20,50,100])*kl_binaries
        kls = [int(kl) for kl in kllist if int(kl) > 0] 
        klslist.append(kls)

        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,8))

        #Jea's preferred formatting
        #plt.rcParams["legend.frameon"] = False
        #plt.rcParams["legend.fontsize"] = 12
        #plt.rcParams["legend.borderpad"] = 0.3
        #plt.rcParams["legend.labelspacing"] = 0.3
        #plt.rcParams["legend.handletextpad"] = 0.3
        #plt.rcParams["font.family"] = "serif"
        #plt.rcParams["font.serif"] = "DejaVu Serif"
        #plt.rcParams["font.size"] = 14
        #plt.rcParams["xtick.top"] = True
        #plt.rcParams["ytick.right"] = True
        #plt.rcParams["xtick.direction"] = "in"
        #plt.rcParams["ytick.direction"] = "in"

        n, b1, patches = ax1.hist(mean_sums, range=(0,200), bins=40,density=True, histtype='step', lw=3)
        n, b1, patches = ax1.hist(mean_sums, range=(0,200), bins=40,density=True, alpha=0)
        last = np.where(b1<sumcut)[0][-1]
        plt.setp(patches[:last+1], 'facecolor', 'b', 'alpha', 0.2) 
        ax1.set_xlabel('Avg. Sum of DQ$_{CF}$-DQ$_{HR}$')
        yl = ax1.get_ylim()
        ax1.plot((mean_sum,mean_sum),(0,yl[1]), color='r')
        ax1.set_xlim(0,200)
        ax1.set_ylabel("Density")

        n, b1, patches = ax2.hist(mean_stds, range = (0.1,0.45), bins=40,density=True, lw=3, histtype='step')
        n, b1, patches = ax2.hist(mean_stds, range = (0.1,0.45), bins=40,density=True, alpha=0)
        last = np.where(b1<stdcut)[0][-1]
        plt.setp(patches[:last+1], 'facecolor', 'b', 'alpha', 0.2) 
        ax2.set_xlabel(r'Avg. $\sigma$ of DQ$_{CF}$-DQ$_{HR}$')
        y2 = ax2.get_ylim()
        ax2.plot((mean_std,mean_std),(0,y2[1]), color='r')
        ax2.set_xlim(0.1,0.45)

        n,b1,patches = ax3.hist(mean_meds, range=(-0.05,0.1),bins=40,density=True, lw=3, histtype='step')
        n,b1,patches = ax3.hist(mean_meds, range=(-0.05,0.1),bins=40,density=True, lw=3, alpha=0)
        last = np.where(b1<medcut)[0][-1]
        plt.setp(patches[:last+1], 'facecolor', 'b', 'alpha', 0.2) 
        ax3.set_xlabel(r'Avg. Median of DQ$_{CF}$-DQ$_{HR}$')
        y3 = ax3.get_ylim()
        ax3.plot((mean_med,mean_med),(0,y3[1]), color='r')
        ax3.set_xlim(-0.05,0.1)
        ax3.set_ylabel("Density")

        n, b1, patches = ax4.hist(mean_mdiffs, range =(0,0.9), bins=40,density=True, lw=3, histtype='step')
        n, b1, patches = ax4.hist(mean_mdiffs, range =(0,0.9), bins=40,density=True, lw=3, alpha=0)
        last = np.where(np.abs(b1)<diffcut)[0][-1]
        first = np.where(np.abs(b1)<diffcut)[0][0]
        plt.setp(patches[:last+1], 'facecolor', 'b', 'alpha', 0.2) 
        ax4.set_xlabel(r'Avg. DQ$_{CF}$-DQ$_{HR}$ at Continuum Peak')
        y4 = ax4.get_ylim()
        ax4.plot((mean_mdiff,mean_mdiff),(0,y4[1]), color='r')

        plt.suptitle('Weights: '+str(wts)+', KLs: '+str(kls))
        plt.savefig('wts'+wtstr+'_kls'+klstr+'_metricselect_hists.png')
        plt.show()
    return(wtslist, klslist)

def false_true_compare(false_dir, real_dir, out_dir,  namestr='*', pklstr = '*', exceptstr=False, wts=[1,1,1,1,1,1], kls=[1,2,3,4,5,10,20,50,100], datestr=False, sepkl=False, seppl=False, snt=False, makeplot=True, writefits=False, maxx=25, maxy=25):

    #find all completed pe combos
    dset_stringlist, reallist, falselist = pull_dsets(false_dir, real_dir, out_dir, namestr=namestr, exceptstr=exceptstr)

    f = plt.figure(figsize=(14.9,4*len(dset_stringlist)))
    plt.axis('off')

    #Jea's preferred formatting
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["legend.fontsize"] = 12
    plt.rcParams["legend.borderpad"] = 0.1
    plt.rcParams["legend.labelspacing"] = 0.1
    plt.rcParams["legend.handletextpad"] = 0.1
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "DejaVu Serif"
    plt.rcParams["font.size"] = 12
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    
    #set up gridspec
    nrows = 4*len(dset_stringlist)+1
    outer = gridspec.GridSpec(nrows, 1, wspace=0., hspace=0.)

    #print(reallist,falselist)
    #extract info and find matches 
    for i in np.arange(len(reallist)):
        linefname = reallist[i].split('_')
        #extract object name, date, wl, cut
        line_obj = linefname[1]
        line_date = linefname[2]
        line_wl = linefname[3]
        line_cut = linefname[4]
        line_klip_params = linefname[-2]

        coll_string = line_obj+'_'+line_date+'_wts'+''.join([str(x) for x in wts])+'_kls'+''.join([str(x) for x in kls])+'_sepkl'+str(sepkl)+'_snthresh'+str(snt)
        #if keylist!=False and i>0:
            #if coll_string in keylist:
                #return(dtn)

        match=False
        #look for matching continuum

        for j in np.arange(len(falselist)):
            contfname = falselist[j].split('_')
            cont_obj = contfname[1]
            cont_date = contfname[2]
            cont_wl = contfname[3]
            cont_cut = contfname[4]
            cont_klip_params = contfname[-2]
            if line_obj == cont_obj:
                if line_date == cont_date:
                    if (line_wl == 'Line') and (cont_wl=='Cont'):
                        if line_cut == cont_cut:
                            if line_klip_params == cont_klip_params:
                                match=True
                                falseind=j


        #compute the aggregate parameter quality metric for the real and fake data
        if match==True:
            fake_metric_cube, fake_agg_cube, fake_ann_val, fake_movm_val, fake_metric_scores, fake_metric_fname = find_best_new(falselist[falseind], kls, pedir=false_dir, writestr=False, writefiles=False, weights=wts, outdir=false_dir+'proc/', 
                                                                                                                                oldpe=False, debug=False, smt=3, snrmeth='stdev',separate_planets=seppl, separate_kls=sepkl, snrthresh=snt, maxx=maxx, maxy=maxy)
            real_metric_cube, real_agg_cube, real_ann_val, real_movm_val, real_metric_scores, real_metric_fname = find_best_new(reallist[i], kls, pedir=real_dir, writestr=False, writefiles=False, weights=wts, outdir=real_dir+'proc/', 
                                                                                                                                oldpe=False, debug=False, smt=3, snrmeth='stdev',separate_planets=seppl, separate_kls=sepkl, snrthresh=snt, maxx=maxx, maxy=maxy)
            ##NORMALIZE TO MAX AND SUBTRACT MIN
            fake_cube = fake_metric_cube
            real_cube = real_metric_cube

            #normalize the agg slices for each KL mode
            if sepkl==False:
                kls_forloop=['avg']
            else:
                kls_forloop=kls
            for kl in np.arange(len(kls_forloop)):
            
                fake_cube[-1,:,kl,:,:,:]=fake_cube[-1,:,kl,:,:,:]-np.nanpercentile(fake_cube[-1,:,kl,:,:,:],10)
                real_cube[-1,:,kl,:,:,:]=real_cube[-1,:,kl,:,:,:]-np.nanpercentile(real_cube[-1,:,kl,:,:,:],10)

                #normalize to max
                fake_cube[-1,:,kl,:,:,:]=fake_cube[-1,:,kl,:,:,:]/np.nanpercentile(fake_cube[-1,:,kl,:,:,:],90)
                real_cube[-1,:,kl,:,:,:]=real_cube[-1,:,kl,:,:,:]/np.nanpercentile(real_cube[-1,:,kl,:,:,:],90)
             
            diff = fake_cube-real_cube
            diff_wt = diff*fake_cube[-1,:,:,:,:,:]
            diff_agg = diff[-1,0,0,0,:,:]
            diff_wt_agg = diff_wt[-1,0,0,0,:,:]

            if writefits==True:
                #write out difference cube
                fits.writeto(out_dir+coll_string+'_diff.fits', diff, overwrite=True)

            if makeplot==True:

                #now make weighted plots
                inner = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[i*4:i*4+4], wspace=0, hspace=0)
                ax1 = plt.Subplot(f, inner[0])
                ax2 = plt.Subplot(f, inner[1])
                ax3 = plt.Subplot(f, inner[2])
                ax4 = plt.Subplot(f, inner[3])
                
                f1 = ax1.imshow(fake_cube[-1,0,0,0,:,:], cmap='magma', vmax=0, vmin=1, origin='lower left', alpha=1)
                if i==0:
                    ax1.set_title('Cont. False Planet ADQ (ADQ$_{CF}$)', size=12)
                ax1.set_xlabel('movement')
                ax1.set_ylabel('annuli')
                #ax1.set_facecolor("white")
                f.colorbar(f1, ax=ax1, shrink=0.75, pad=-0.2, alpha=0)#, label='Relative Parameter Quality')
                if datestr!=False:
                    ax1.text(-7,4, line_obj + ' ' + str(datestr[i]), rotation=90, fontsize=14)
                else:
                    ax1.text(-7,4, line_obj + ' ' + line_date, rotation=90, fontsize=14)
                if i<4:
                    ax1.axes.get_xaxis().set_visible(False)

                ind = np.where(fake_cube[-1,0,0,0,:,:] == np.nanmax(fake_cube[-1,0,0,0,:maxy,:maxx]))
                label_text = 'a' + str(int(fake_ann_val[0][0])) + 'm' + str(int(fake_movm_val[0][0]))
                rect = patches.Rectangle((ind[1][0] - 0.5, ind[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='k', facecolor='none')
                ax1.add_patch(rect)
                ax1.text(ind[1][0] + 0.75, ind[0][0]-0.5, label_text, color='black')

                f.add_subplot(ax1)

                f2 = ax2.imshow(real_cube[-1,0,0,0,:,:], cmap='magma', vmax=0, vmin=1, origin='lower left')
                if i==0:
                    ax2.set_title(r'H$\alpha$ Real Object DQ (ADQ$_{HR}$)', size=12)
                if i==4:
                    ax2.set_xlabel('movement')

                f.colorbar(f2, ax=ax2, shrink=0.75, pad=-0.2,  label='Relative Parameter Quality')
                ind2 = np.where(real_cube[-1,0,0,0,:,:] == np.nanmax(real_cube[-1,0,0,0,:maxy,:maxx]))
                label_text2 = 'a' + str(int(real_ann_val[0][0])) + 'm' + str(int(real_movm_val[0][0]))
                rect2 = patches.Rectangle((ind2[1][0] - 0.5, ind2[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='r', facecolor='none')
                ax2.add_patch(rect2)
                if i in [1,3,4]:
                    ax2.text(ind2[1][0] + 0.75, ind2[0][0]-0.5, label_text2, color='red')
                else:
                    ax2.text(ind2[1][0] - 2.5, ind2[0][0]+1, label_text2, color='red')
                rect4 = patches.Rectangle((ind[1][0] - 0.5, ind[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='k', facecolor='none')
                ax2.add_patch(rect4)
                ax2.text(ind[1][0] + 0.75, ind[0][0]-0.5, label_text, color='black')

                if i<4:
                    ax2.axes.xaxis.set_ticks([])
                    ax2.axes.yaxis.set_ticks([])
                else:
                    ax2.axes.get_yaxis().set_visible(False)

                f.add_subplot(ax2)

                f3 = ax3.imshow(diff_wt_agg, cmap='coolwarm', vmin=-0.3, vmax=0.3, origin='lower left')
                if i==0:
                    ax3.set_title('Parameter Quality Difference', size=12)
                if i==4:
                    ax3.set_xlabel('movement')
                #ax3.set_ylabel('annuli')
                f.colorbar(f3, ax=ax3, shrink=0.75, pad=-0.2)#, label = r'(ADQ$_{CF}$ - DQ$_{HR}$)$\times$ADQ$_{CF}$')
                rect3 = patches.Rectangle((ind[1][0] - 0.5, ind[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='k', facecolor='none')
                ax3.add_patch(rect3)
                ax3.text(ind[1][0] + 0.75, ind[0][0]-0.5, label_text, color='black')
                if i<4:
                    ax3.axes.xaxis.set_ticks([])
                    ax3.axes.yaxis.set_ticks([])
                else:
                    ax3.axes.yaxis.set_ticks([])
                    ax3.axes.get_yaxis().set_visible(False)

                f.add_subplot(ax3)

                ax4.hist(diff_wt_agg.flatten(), range=(-0.5,0.5), color='k',  label='all',histtype='step', bins=20, density=True)
                if i==0:
                    ax4.set_title('Norm. Diff. Values', size=12)
                ax4.set_xlabel('ADQ$_{CF}$ - DQ$_{HR}$')
                ax4.set_ylabel('density')
                ax4.yaxis.set_label_position("right")
                ax4.yaxis.tick_right()
                ymax = ax4.get_ylim()
                ax4.yaxis.set_ticks(np.arange(1,ymax[1]))
                if i<4:
                    ax4.axes.xaxis.set_ticks([])
                    ax4.axes.get_xaxis().set_visible(False)
                ax4.set_xlim(-0.5,0.5)
                box = ax4.get_position()
                box.x0 = box.x0 + 0.008
                box.x1 = box.x1 
                ax4.set_position(box)
                ax4.set_aspect(1.05/ax4.get_data_ratio(), adjustable='box')
                plt.legend()

                f.add_subplot(ax4)

    if makeplot==True:
        lastrow = gridspec.GridSpecFromSubplotSpec(1, 100, subplot_spec=outer[-1], wspace=0, hspace=0)
        ax7 = plt.Subplot(f, lastrow[3:50])
        ax8 = plt.Subplot(f, lastrow[55:75])

        #colorbars along the bottom
        f.colorbar(f1, ax=ax7, fraction=0.35, orientation='horizontal', label='Relative Parameter Quality')
        f.colorbar(f3, ax=ax8,  fraction=0.35, orientation='horizontal', label = r'(ADQ$_{CF}$ - DQ$_{HR}$)$\times$ADQ$_{CF}$', aspect=10)

        size = f.get_size_inches()

        plt.savefig(out_dir+coll_string+'_realvfalse.png')

    return()


