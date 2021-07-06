import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
from numpy.core.numeric import _outer_dispatcher
import pandas as pd
from traitlets.traitlets import parse_notifier_name
import gaplanets as gp
from importlib import reload
import re
from scipy import ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import glob
import pyklip.klip as klip
import pyklip.instruments.MagAO as MagAO
import pyklip.parallelized as parallelized
import SNRMap_new as snr
import os
import pdb
import pickle
import textwrap


class CollapsedPE():
    """
    Collapse the parameter explorer
    """

    def __init__(self, pename, pedir='./', outdir='./',  kllist = [5,10,20,50],  calflux = False, datadir = './', writestr=False, snrthresh=False, oldpe=False, separate_planets=False, iwa = 0,  wts = [1,1,1,1,1,1,1,1],
     hpval = 5, header = True, mode = 'Line', xname = '',  collmode = None, separate_kls = False, owa = None, savefig = True, snrmeth = 'stdev'):
        
        """
        Averages over the planet dimension of a parameter explorer file

        Args:
            pename (str):     name of paramexplore file
            pedir (str, optional:      directory holding parameter explorer file
            writestr (str, otpional):   filename prefix for saved file. If not specified, preserves name of parameter explorer
            snrthresh (int, optional):  SNR threshhold for peak pixel under mask. All lower values will be masked.
        
        Returns:
            pecube:     parameter explorer matrix with planet dimension collapsed to 1
            writename:  name of output file

        """

        self.pename = pename
        self.pedir = pedir
        self.outdir = outdir
        self.writestr = writestr
        self.snrthresh = snrthresh
        self.calflux = calflux
        self.oldpe = oldpe
        self.separate_planets = separate_planets
        self.iwa = iwa
        self.wts = wts
        self.hpval = hpval
        self.header = header
        self.mode = mode
        self.savefig = savefig
        self.xname = xname
        self.collmode = collmode
        self.separate_kls = separate_kls
        self.owa = owa
        self.header = header
        self.kllist = kllist
        self.datadir = datadir
        self.snrmeth = snrmeth




    def collapse_planets(self):
        """
        Averages over the planet dimension of a parameter explorer file

        Args:
            pename (str):     name of paramexplore file
            pedir (str, optional:      directory holding parameter explorer file
            writestr (str, otpional):   filename prefix for saved file. If not specified, preserves name of parameter explorer
            snrthresh (int, optional):  SNR threshhold for peak pixel under mask. All lower values will be masked.
        
        Returns:
            pecube:     parameter explorer matrix with planet dimension collapsed to 1
            writename:  name of output file

        """

        if self.writestr is False:
            #use the parameter explorer name for this file as well, minus the '_highpass_klmodes-all.fits'
            self.writestr = self.pename[:-17]

        # read in image and header
        pecube = fits.getdata(self.pedir + self.pename)
        pehead = fits.getheader(self.pedir + self.pename)
        dims = pecube.shape
        nplanets=dims[3]
        pehead["NPLANET"]=nplanets
        pehead["PLSEP"]=self.separate_planets

        #if snrthresh set, find values where peak pixel SNRs (slices 1/2) are below threshhold and replace with nans
        if (self.snrthresh != False):
            #stdev maps
            ind = np.where(pecube[0,:,:,:,:,:] < self.snrthresh)
            lowsnr_mask=np.ones(dims[1:])
            lowsnr_mask[ind]=np.nan
            #apply to peak (sl 0) and avg under mask (sl 2)
            for sl in [0,2]:
                pecube[sl,:,:,:,:]*=lowsnr_mask
            
            #absmed maps
            ind_absmed = np.where(pecube[1,:,:,:,:,:] < self.snrthresh)
            lowsnr_mask_absmed=np.ones(dims[1:])
            lowsnr_mask_absmed[ind_absmed]=np.nan
            #apply to peak (sl 1) and avg under mask (sl 3)
            for sl in [1,3]:
                pecube[sl,:,:,:,:]*=lowsnr_mask_absmed

       
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
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        #collapse in planet dimension or extract 
        #mean preserves nans so all planets have to be recovered
        #separates planets
        if self.separate_planets==True: 
            npldim=nplanets
            writename= self.writestr+'_sepplanets.fits'
        
        #collapses planet dimension
        else: 
            writename = self.writestr+'_planetcollapse.fits'
            #quick hack for old PEs where planet dimension was 6th
            if self.oldpe==True:
                pecube=np.mean(pecube,axis=5, keepdims=True)
            else:
                pecube=np.mean(pecube,axis=3, keepdims=True)
            npldim=1
        
        fits.writeto(self.outdir+writename, pecube, pehead, overwrite=True)

        return pecube, writename, npldim


    def trimkl(self):
        """
        Reads in a parameter explorer file and collapses it in the KLmode dimension (axis 3)

        Args:
            pename (str):     name of paramexplore file
            kllist (list:int):     list of KL modes you want to collapse over
            pedir (str: optional):      directory holding parameter explorer file
            writestr (str: optional):   filename prefix for saved file. If not specified, preserves name of parameter explorer

        Returns:
            Trimmed KL Cube 

        """

        if self.writestr == False:
            self.writestr = self.pename[:-5]

        writename= self.writestr + '_trimmedkls.fits'

        # read in image and header
        klcube_raw = fits.getdata(self.pedir + self.pename)
        head = fits.getheader(self.pedir + self.pename)

        dims=list(klcube_raw.shape)

        #create an array with same dims except for kl modes
        
        dims[2]=len(self.kllist)
        klkeep = np.zeros(dims)
            
        # pull KL modes of parameter explorer from the header
        allkl = list(map(int, head['KLMODES'][1:-1].split(",")))

        # find only the KL slices in the list of modes to be collapsed and make an array
        i = 0;
        j = 0;
        keepind = []
        for kl in allkl:
            if kl in self.kllist:
                keepind.append(i)
                print('keeping kl', kl, "with index", i)
                #only keep KL modes matching kllist
                klkeep[:, :,j,:,:,:] = klcube_raw[:,:,i,:,:,:]
                j += 1
            i += 1

        # replacing -250 (value when planet too close to annulus) with 0 for ease of display
        klkeep[klkeep == -1000] = np.nan

        # update header to reflect KL modes used
        head["KLMODES"] = str(self.kllist)

        fits.writeto(self.outdir + writename, klkeep, head, overwrite=True)

        #return trimmed cubes
        return klkeep, writename


    def __filter_nan_gaussian_conserving(arr, sigma):
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
        loss = ndimage.gaussian_filter(loss, sigma=sigma, mode='constant', cval=1)

        gauss = arr.copy()
        gauss[nan_msk] = 0
        gauss = ndimage.gaussian_filter(gauss, sigma=sigma, mode='constant', cval=0)
        gauss[nan_msk] = np.nan

        gauss += loss * arr

        return gauss


    def find_best_new(self, debug=False, smt=3):
        """
        Collapses parameter explorer file and extracts the optimal parameter value

        Args:
            smt (int: optional):        FWHM of gaussian smoothing kernel for 'neighbor quality' metrics
            snrmeth (str):    one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value),
                        and 'all' to average the two


        Returns:
            metric_cube 
            agg_cube 
            ann_val 
            movm_val 
            metric_scores 
            metric_fname


        """

        #EXTRACT PLANETS OR COLLAPSE

        pecube, plwritename, npldim = self.collapse_planets(self.pename, pedir=self.pedir, outdir=self.outdir, snrthresh=self.snrthresh, oldpe=self.oldpe, writestr=self.writestr, separate_planets=self.separate_planets)

        #EXTRACT KL MODES OR COLLAPSE
        kltrim, writename = self.trimkl(plwritename, self.kllist, pedir=self.outdir, outdir=self.outdir)

        if self.writestr == False:
            writestr = writename[:-5]

        # If collapsing, make mean and stdev arrays
        if self.separate_kls == False:
            stdevkl = np.std(kltrim, axis=2, keepdims=True)

            #overwrite kltrim with average
            kltrim= np.mean(kltrim, axis=2, keepdims=True)

            #grab header
            head = fits.getheader(self.outdir+writename)
            head["KLCOLL"]='True'

            # write arrays
            fits.writeto(self.outdir+ writestr + '_avgkl.fits', kltrim, head, overwrite=True)
            fits.writeto(self.outdir+ writestr + '_stdevkl.fits', stdevkl, head, overwrite=True)
            #fits.writeto(outdir+writestr + '_sumkl.fits', sumkl, head, overwrite=True)

        ##EXTRACT SNR SLICES
        #pull the SNR map slices with the matching snrmeth value
        if self.snrmeth == "absmed":
            slice = 1
        if self.snrmeth == "stdev":
            slice = 0  

        dims = kltrim.shape

        #extract the appropriate slices
        if self.snrmeth is not 'all':
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
        if len(self.kllist)>1 and self.separate_kls==False:
            klloop = 1
            kllist=[self.kllist]
            stdev_valid=True
        else: 
            klloop = len(self.kllist)
            stdev_valid=False

        #set up cubes for planet, kl loop

        #extracts parameter explorer inputs from the filename (assumes was autonamed by parameter explorer)
        namelist = self.pename.split('_')
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
        nq_stdev = np.zeros([nstepy, nstepx])
        agg_cube=np.zeros([klloop,npldim,nstepy,nstepx])
        agg=np.zeros([nstepy,nstepx])
        ann_val = np.zeros([klloop,npldim])
        movm_val = np.zeros([klloop,npldim])
        ##Note - hard coded for 1 subsection. 
        metric_cube = np.zeros([9, 1, klloop, npldim, nstepy, nstepx])
        qual_cube = np.zeros([10, 1, klloop, npldim, nstepy, nstepx])

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

                #normalize the SNR (where high values = good) 
                #note - hard-coded for susections = 1
                snr_norm = kltrim_snr[0,0,k,p,:,:] / np.nanmax(kltrim_snr[0,0,k,p,:,:])
                snr_norm_umask = kltrim_snr[1,0,k,p,:,:] / np.nanmax(kltrim_snr[1,0,k,p,:,:])

                if stdev_valid==True:
                    #normalize standard deviations across KL modes. Low values = good
                    # divide by SNR so is Stdev in SNR as fraction of SNR itself
                    stdev_norm_cube = stdevkl[0:2,0,k,p,:,:] / kltrim_snr[0:2,0,k,p,:,:]
                    #first slice is peak
                    stdev_norm = 1 - (stdev_norm_cube[0,:,:]/np.nanmax(stdev_norm_cube[0,:,:]))
                    #second slice is under mask
                    stdev_norm_umask = 1 - (stdev_norm_cube[1,:,:]/np.nanmax(stdev_norm_cube[1,:,:]))

                #if extracting separately, fill these arrays with nans 
                else:
                    stdev_norm = np.zeros([nstepy,nstepx])*np.nan
                    stdev_norm_umask = np.zeros([nstepy,nstepx])*np.nan

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
                nq_snr = self.__filter_nan_gaussian_conserving(snr_norm,sig)
                nq_stdev = self.__filter_nan_gaussian_conserving(stdev_norm,sig)
                nq_snr_umask = self.__filter_nan_gaussian_conserving(snr_norm_umask,sig)
                nq_stdev_umask = self.__filter_nan_gaussian_conserving(stdev_norm_umask,sig)

                #normalizes neighbor quality
                nq_snr /= np.nanmax(nq_snr)
                nq_stdev /= np.nanmax(nq_stdev)
                nq_snr_umask /= np.nanmax(nq_snr_umask)
                nq_stdev_umask /= np.nanmax(nq_stdev_umask)

                if debug==True:
                    #make a cube of all these metrics for sanity checking
                    qual_cube[0,:,k,p,:,:]=snr_norm
                    qual_cube[1,:,k,p,:,:]=snr_norm_umask
                    qual_cube[2,:,k,p,:,:]=stdev_norm
                    qual_cube[3,:,k,p,:,:]=stdev_norm_umask
                    qual_cube[4,:,k,p,:,:]=spurpix
                    qual_cube[5,:,k,p,:,:]=nq_snr
                    qual_cube[6,:,k,p,:,:]=nq_snr_umask
                    qual_cube[7,:,k,p,:,:]=nq_stdev
                    qual_cube[8,:,k,p,:,:]=nq_stdev_umask
                    qual_cube[9,:,k,p,:,:]=contrast

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
                metric_cube[4,:,k,p,:,:]= stdev_norm_umask
                metric_cube[5,:,k,p,:,:]= nq_stdev_umask

                #spurious pixel metric = 1 if no spurious pixels and 0 if max number for this dataset
                if np.nanmax(spurpix)>0:
                    spurpix_norm = 1-spurpix/np.nanmax(spurpix)
                else: #edge case - no spurious pixels in any image
                    spurpix_norm= 1+spurpix
                
                metric_cube[6,:,k,p,:,:]= spurpix_norm
                metric_cube[7,:,k,p,:,:]= contrast 

                metriclist = (snr_norm, nq_snr, snr_norm_umask, nq_snr_umask, stdev_norm_umask, nq_stdev_umask, spurpix_norm, contrast)
                
                #make sure weights for stdev slices are 0 if only 1 kl mode or extracting separately
                if stdev_valid==False:
                    self.wts[4]=0
                    self.wts[5]=0

                #calculate an aggregate parameter quality metric by summing the individual metrics * their weights
                for metricind in np.arange(len(metriclist)):
                    #only sum if non-zero (otherwise nans will carry over)
                    if self.wts[metricind]>0:
                        agg+= self.wts[metricind]*metriclist[metricind]

                metric_cube[8,:,k,p,:,:]=agg
                agg_cube[k,p,:,:]=agg

                ##find location or peak of parameter quality metric and print info
                ind = np.where(agg == np.nanmax(agg))

                if agg[ind[0],ind[1]].shape[0]>1:
                    print("the optimal solution for this choice of parameters/weights is not unique")
                    return()

                #extract metric scores for this location    
                metric_scores = [snr_norm[ind][0], nq_snr[ind][0], snr_norm_umask[ind][0], nq_stdev_umask[ind][0], \
                stdev_norm_umask[ind][0], nq_stdev_umask[ind][0], spurpix_norm[ind][0], agg[ind][0]]

                #translate to annuli and movement values
                ann_val[k,p]= ymin + ind[0][0] * ystep
                movm_val[k,p] = xmin + ind[1][0] * xstep
                
                if self.separate_planets==False:
                    plno = 'all'
                else:
                    plno = p+1
                print('peak for planet =', plno, 'klmode = ', kllist[k], 'is at', ind[0][0], ind[1][0], 'corresponding to annuli', ann_val[k][p], ' and movement', movm_val[k][p])
                #print('SNR value for fake planets (avg of SNR methods and planets) is', avgSNR)
                #print('metric scores for (snr peak, snr peak neigbors, snr umask, snr umask neighbors, stdev, stdev neighbors, spurious pix, contrast, agg) are:', metric_scores)



        if debug==True:
            fits.writeto(self.outdir+writename[:-5]+'_paramqual_fullcube.fits', qual_cube, overwrite=True)

        metric_fname = writename[:-5]+'_paramqual_metrics.fits'
        fits.writeto(self.outdir+metric_fname, metric_cube, overwrite=True)
        
        return metric_cube, agg_cube, ann_val, movm_val, metric_scores, metric_fname


    def collapse_pes(self, smt = 3):
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
        klliststr="_".join([str(n) for n in self.kllist])
        wtstr = re.sub(r'\W+', '', str(self.wts))
        xstr = '_'+wtstr+'_'+ self.snrmeth+'_'+klliststr+ self.xname

        #find all the parameter explorer files in the specified directory
        rundir = os.getcwd()
        os.chdir(self.pedir)
        flist = glob.glob('paramexplore*KLmodes-all.fits')
        
        os.chdir(rundir)

        #sort alphabetically so the same objects are together
        flist = sorted(flist)

        #set up dictionary
        npes=len(flist)
        d={}

        #store some global values for this dictionary
        d["dict_name"]= self.xname
        d["weights"]= self.wts
        d["snrmeth"]= self.snrmeth
        d["snr_threshhold"]= self.snrthresh
        d["mode"]= self.mode
        d["datadir"]= self.datadir
 

        #create directory to save in if doesn't yet exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        #if highpass, collapse_mode, owa are not lists, repeat the entered value
        if not isinstance(self.hpval, list):
            hpval = np.repeat(self.hpval, npes)
        if not isinstance(self.collmode, list):
            collmode = np.repeat(self.collmode, npes)
        if not isinstance(self.owa, list):
            owa=np.repeat(self.owa, npes)
     
        if not isinstance(self.calflux, list):
            calflux=np.repeat(self.calflux,npes)

        #loop over number of PEs found 
        for i in np.arange(npes):
            # creates a dictionary with parameters needed to run KLIp for all images in a dataset
            #store PE file
            d["pe{0}".format(i+1)] = fits.getdata(self.pedir+flist[i])

            #store PE name
            pename = flist[i]
            d["pe{0}name".format(i+1)] = pename

            #extract IWA and store
            split=flist[i].split('_')
            s = re.findall(r'iwa(.*?)_', flist[i])

            d['pe{0}iwa'.format(i+1)]=int(s[0])

            head = fits.getheader(self.pedir+flist[i])
            file = fits.getdata(self.pedir+flist[i])
            nplanets = file.shape[3]
            print(file.shape, nplanets)

            #if header keyword set, pull other KLIP values from header
            if self.header == True:
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
                d['pe{0}fpath'.format(i+1)]=self.datadir+split[1]+'/'+split[2]+'_'+split[3]+'/'
                d["pe{0}cut".format(i+1)]=split[5]
                #d["pe{0}pfx".format(i+1)]= split[1]+'_'+split[2]+'_'+split[3]+'_'+split[4]+'_'+split[5]
            else:
                d["pe{0}cut".format(i+1)]=split[4]
                d['pe{0}fpath'.format(i+1)]=self.datadir+split[1]+'/'+split[2]+'/'
                #d["pe{0}pfx".format(i+1)]= split[1]+'_'+split[2]+'_'+split[3]+'_'+split[4]
            
            #store dataset name for plot titles later
            #remove first underscore
            pfx=pfx[1:]
            d["pe{0}pfx".format(i+1)]=pfx
            d["pe{0}dset".format(i+1)] = dsetname
            writename = d["pe{0}pfx".format(i+1)]+xstr

            print(pename)

            ## extract best parameters according to weights and methods
            metric_cube, agg_cube, ann_val, movm_val, metric_scores, metric_fname = self.find_best_new(self.pename, self.kllist, pedir=self.pedir, outdir=self.outdir, \
                writestr=writename, weights=self.wts, snrmeth=self.snrmeth, smt=smt, snrthresh=self.snrthresh, separate_planets=self.separate_planets, separate_kls=self.separate_kls)

            d["pe{0}ann".format(i+1)]=ann_val
            d["pe{0}movm".format(i+1)]=movm_val
            d["pe{0}agg".format(i+1)]=agg_cube
            nkldim = ann_val.shape[0]
            npldim = ann_val.shape[1]
            strklip = [['']*npldim]*nkldim

            #save visualization of the metrics for this PE explorer
            if self.savefig==True:
                self.paramexplore_fig(self.pename, self.kllist, pedir=self.pedir, outdir= self.outdir, weights= self.wts, snrmeth= self.snrmeth, smt= self.smt)
        
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
                    if self.mode=='SDI':
                        modes=['Cont', 'Line']
                    else:
                        runmode=self.mode
                        modes=['single']
                    
                    for j in np.arange(len(modes)):
                        #print(haindir, contindir)
                        if self.mode=='SDI':
                            runmode=modes[j]
                        if runmode=='Cont':
                            filelist = glob.glob(contindir + "sliced*.fits")
                        if runmode=='Line':
                            #default prefix has Cont in it, since we ran the PE on fake planets in cont images
                            prefix = prefix.replace('Cont','Line') 
                            filelist = glob.glob(haindir + "sliced*.fits")
                    
                        #run KLIP with the optimum values for each PE

                        #check if file with this name already exists
                        if os.path.exists("{out}/{pre}-KLmodes-all.fits".format(out=self.outdir, pre=prefix+strklip[kl][pl])):
                            print("This file already exists. I am NOT re-running KLIP, but just reading the existing image in. Check and make sure you weren't intending to change the name")
                        #otherwise, run KLIP
                        else:
            
                            dataset = MagAO.MagAOData(filelist) 
                            dataset.IWA=iwa
                            dataset.OWA=float(owa[i])

                            print('running KLIP. highpass = ', hpval[i], ', calflux = ', calflux[i], ', time collapse = ', collmode[i], ', OWA = ', owa[i], 'prefix =', prefix, 'first file =', filelist[0])
                            parallelized.klip_dataset(dataset, outputdir=self.outdir, fileprefix=prefix+strklip[kl][pl], 
                                    algo='klip', annuli=int(ann_val[kl,pl]), subsections=1, movement=movm_val[kl,pl],
                                    numbasis=self.kllist, maxnumbasis=100, calibrate_flux=calflux[i], mode="ADI", highpass=float(hpval[i]), 
                                    save_aligned=False, time_collapse=collmode[i])
                            #fits.writeto('test'+str(i)+'.fits', dataset.output, overwrite=True)

                        #pull output image 
                        klim = fits.getdata("{out}/{pre}-KLmodes-all.fits".format(out=self.outdir, pre=prefix+strklip[kl][pl]))
                        print(klim.shape)
                        if kl==0 and pl==0:
                            klcube = np.zeros([nkldim, npldim, klim.shape[0],klim.shape[1],klim.shape[2]])
                        klcube[kl,pl,:,:,:]=klim

                        #and store in dictionary (overwrites at present until end of loop - needs work)
                        d["pe{0}strklip".format(i+1)]=strklip
                        if runmode=='Cont':
                            d["pe{0}contklipim".format(i+1)]=klcube
                        if runmode=='Line':
                            d["pe{0}haklipim".format(i+1)]=klcube
        return d

    def paramexplore_fig(self):
    
        metric_cube, agg_cube, ann_val, movm_val, metric_scores, metric_fname = self.find_best_new(self.pename, self.kllist, pedir=self.pedir, 
            outdir=self.outdir, writestr=self.writestr, weights=self.weights, snrmeth=self.snrmeth, smt=self.smt)

        if self.writestr == False:
            writestr = self.pename[:-17]

        namelist = self.pename.split('_')
        params = [s for s in namelist if s[0] == 'a']
        params = params[0]
        params = re.split('a|-|x|m|s|iwa', params)
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
        gs = fig.add_gridspec(2, 6)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[1,0])
        ax4 = fig.add_subplot(gs[1,1])
        ax5 = fig.add_subplot(gs[0,2])
        ax6 = fig.add_subplot(gs[0,3])
        ax7 = fig.add_subplot(gs[1,2])
        ax8 = fig.add_subplot(gs[1,3])
        ax9 = fig.add_subplot(gs[:,4:])

        plt.setp((ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9), xticks=np.arange(0, nstepx + 1, 2), xticklabels=np.arange(xmin, xmax + 1, 2),
                    yticks=np.arange(0, nstepy + 1, 2), yticklabels=np.arange(ymin, ymax + 1, 2))

        #if extracted kl or planets separately, make figs separately
        nkldim = ann_val.shape[0]
        npldim = ann_val.shape[1]

        for kl in np.arange(nkldim):
            for pl in np.arange(npldim):

                im1 = ax1.imshow(metric_cube[0,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
                ax1.set_xlabel("movement parameter")
                ax1.set_ylabel("annuli parameter")
                ax1.set_title("Peak SNR: Weight = "+str(self.wts[0]))
                divider = make_axes_locatable(ax1)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(im1, cax=cax, orientation='vertical')

                im2 = ax2.imshow(metric_cube[1,0,kl,pl,:,:], origin='lower',cmap='magma', vmin=0, vmax=1)
                ax2.set_xlabel("movement parameter")
                ax2.set_ylabel("annuli parameter")
                ax2.set_title("Peak SNR Neighbor Quality: Weight = "+str(self.wts[1]))
                divider = make_axes_locatable(ax2)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(im2, cax=cax, orientation='vertical')

                im3 = ax3.imshow(metric_cube[2,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
                ax3.set_xlabel("movement parameter")
                ax3.set_ylabel("annuli parameter")
                ax3.set_title("Avg SNR Under Mask: Weight = "+str(self.wts[2]))
                divider = make_axes_locatable(ax3)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(im3, cax=cax, orientation='vertical')

                im4 = ax4.imshow(metric_cube[3,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
                ax4.set_xlabel("movement parameter")
                ax4.set_ylabel("annuli parameter")
                ax4.set_title("Avg SNR Neighbor Quality: Weight = "+str(self.wts[3]))
                divider = make_axes_locatable(ax4)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(im4, cax=cax, orientation='vertical')

                
                im5 = ax5.imshow(metric_cube[4,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
                ax5.set_xlabel("movement parameter")
                ax5.set_ylabel("annuli parameter")
                ax5.set_title("Stdev Across KL: Weight = "+str(self.wts[4]))
                divider = make_axes_locatable(ax5)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(im5, cax=cax, orientation='vertical')

                im6 = ax6.imshow(metric_cube[5,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
                ax6.set_xlabel("movement parameter")
                ax6.set_ylabel("annuli parameter")
                ax6.set_title("Stdev Neighbor Quality: Weight = "+str(self.wts[5]))
                divider = make_axes_locatable(ax6)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(im6, cax=cax, orientation='vertical')

                im7 = ax7.imshow(metric_cube[6,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
                ax7.set_xlabel("movement parameter")
                ax7.set_ylabel("annuli parameter")
                ax7.set_title("Spurious Pixels: Weight = "+str(self.wts[6]))
                divider = make_axes_locatable(ax7)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(im7, cax=cax, orientation='vertical')


                im8 = ax8.imshow(metric_cube[7,0,kl,pl,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
                ax8.set_xlabel("movement parameter")
                ax8.set_ylabel("annuli parameter")
                ax8.set_title("Contrast: Weight = "+str(self.wts[7]))
                divider = make_axes_locatable(ax8)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                plt.colorbar(im8, cax=cax, orientation='vertical')

                # plot metric
                im9 = ax9.imshow(agg_cube[kl,pl,:,:], origin='lower', vmin=0, vmax=np.sum(self.wts))
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

                plt.savefig(self.outdir+writestr+'_kl'+str(self.kllist[kl])+'_pl'+str(pl)+'_paramqual.png')
        
        return