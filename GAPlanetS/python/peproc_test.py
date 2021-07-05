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

    def __init__(self, pename, pedir='./', outdir='./', writestr=False, snrthresh=False, oldpe=False, separate_planets=False, iwa = 0):
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
        self.oldpe = oldpe
        #self.separate_planets = separate_planets
        self.iwa = iwa


    def collapse_all(self):
        """
        Collapses PE's across all KL modes and Planets
        """

        if self.writestr == False:
        # Use the parameter explorer name for this minus the '.fits'
            writestr = self.pename[:-5]

        # Read in image and header
        pecube = fits.getdata(self.pedir + self.pename)
        pehead = fits.getheader(self.pedir + self.pename)
        dims = pecube.shape
        nplanets=dims[3]
        pehead["NPLANET"] = nplanets
        pehead["PLSEP"]= separate_planets

        # If snrthresh set, find values where peak pixel SNRs (slices 1/2) are below threshhold and replace with nans
        if (self.snrthresh is not False):
            ind = np.where(pecube[0,:,:,:,:,:] < self.nrthresh)
            lowsnr_mask = np.ones(dims[1:])
            lowsnr_mask[ind] = np.nan

            # Apply to peak (slice 0) and avg under mask (slice 2)
            for sl in [0,2]:
                pecube[sl,:,:,:,:]*=lowsnr_mask
            

        
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
        
        fits.writeto(outdir+writename, pecube, pehead, overwrite=True)

        return (pecube, writename, npldim)

    def trimkl(pename, kllist, pedir='./', outdir='proc/', writestr=False):
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

        if writestr == False:
            writestr = pename[:-5]

    writename=writestr + '_trimmedkls.fits'

    # read in image and header
    klcube_raw = fits.getdata(pedir + pename)
    head = fits.getheader(pedir + pename)

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
            print('keeping kl', kl, "with index", i)
            #only keep KL modes matching kllist
            klkeep[:, :,j,:,:,:] = klcube_raw[:,:,i,:,:,:]
            j += 1
        i += 1

    # replacing -250 (value when planet too close to annulus) with 0 for ease of display
    klkeep[klkeep == -1000] = np.nan

    # update header to reflect KL modes used
    head["KLMODES"] = str(kllist)

    fits.writeto(outdir+writename, klkeep, head, overwrite=True)

    #return trimmed cubes
    return (klkeep, writename)
