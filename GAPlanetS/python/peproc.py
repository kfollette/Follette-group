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
import glob
import pyklip.klip as klip
import pyklip.instruments.MagAO as MagAO
import pyklip.parallelized as parallelized
import SNRMap_new as snr
import os
import pdb
import pickle
import textwrap


def collapse_planets(pename, pedir='./', outdir='proc/', writestr=False, snrthresh=False, oldpe=False, separate=False):
	"""
	Averages over the planet dimension of a parameter explorer file

	REQUIRED INPUTS:
	pename: 	name of paramexplore file

	OPTIONAL INPUTS
	pedir:  	directory holding parameter explorer file
	writestr: 	filename prefix for saved file (suffix is _planetcollapse). 
				If not specified, preserves name of parameter explorer
	snrthresh:	SNR threshhold for peak pixel under mask. All lower values will be masked.
 
   RETURNS:
	pecube:		parameter explorer matrix with planet dimension collapsed to 1
	writename:	name of output file

	"""

	if writestr == False:
		#use the parameter explorer name for this file as well, minus the '_highpass_klmodes-all.fits'
		writestr = pename[:-17]

	# read in image and header
	pecube = fits.getdata(pedir + pename)
	pehead = fits.getheader(pedir + pename)
	dims = pecube.shape


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
	if separate==True: 
		writename = writestr+'_pl'
		for i in np.arange(dims[-1]):
			fits.writeto(outdir+writename+str(i+1)+'.fits', pecube[:,:,i,:,:,:], pehead, overwrite=True)
	
	#collapses planet dimension
	else: 
		writename = writestr+'_planetcollapse.fits'
		#quick hack for old PEs where planet dimension was 6th
		if oldpe==True:
			pecube=np.mean(pecube,axis=5)
		else:
			pecube=np.mean(pecube,axis=3)
	
	fits.writeto(outdir+writename, pecube, pehead, overwrite=True)
			
	return (pecube, writename)

def collapsekl(pename, kllist, pedir='./', outdir='proc/', snrmeth='stdev', writestr=False):
	"""
	Reads in a parameter explorer file and collapses it in the KLmode dimension (axis 3)

	REQUIRED INPUTS:
	pename: 	name of paramexplore file
	kllist: 	list of KL modes you want to collapse over

	OPTIONAL INPUTS:
	pedir: 		directory holding parameter explorer file
	snrmeth: 	one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value),
				and 'all' to average the two
	writestr: 	filename prefix for saved file (suffixes are _avgkl and _stdevkl). 
				If not specified, preserves name of parameter explorer

	RETURNS:
	avgkl:  	array containing averages over the specified KL modes
	stdevkl: 	array containing standard deviations over the specified KL modes
	"""

	if writestr == False:
		writestr = pename[:-20]

	# read in image and header
	klcube_raw = fits.getdata(pedir + pename)
	head = fits.getheader(pedir + pename)

	#pull only the SNR map slice with the matching snrmeth value
	if snrmeth == "absmed":
		slice = 1
	if snrmeth == "stdev":
		slice = 0  
	dims = klcube_raw.shape

	#if only 1 metric
	if snrmeth != 'all':
		#for single method cubes, slices are SNR peak, avg SNR, total >thresh pixels, >thresh pixels inside CR
		print('keeping only', snrmeth, 'maps')
		klcube_trimmed=klcube_raw[slice::2,:,:,:,:]

		#if snrmeth = absmed, grab contrast slice too
		if slice==1:
			dmns = np.array(klcube_trimmed.shape)
			dmns[0]+=1
			klcube_plusone = np.zeros(dmns)
			klcube_plusone[0:-1,:,:,:,:]=klcube_trimmed
			#fill last one with contrast cube dimension, which is not dependent on snr method
			klcube_plusone[-1,:,:,:,:]=klcube_raw[-1,:,:,:,:]
			klcube_trimmed = klcube_plusone

		klkeep = np.zeros([5, dims[3], dims[4], len(kllist)])
	else:
		klkeep = np.zeros([9, dims[3], dims[4], len(kllist)])
		
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
			klkeep[:, :, :, j] = klcube_trimmed[:, 0, i, :, :]
			j += 1
		i += 1

	# replacing -250 (value when planet too close to annulus) with 0 for ease of display
	klkeep[klkeep == -1000] = np.nan

	# make mean and stdev arrays
	# preserve nans so non-recoveries are removed
	avgkl = np.mean(klkeep, axis=3)
	stdevkl = np.std(klkeep, axis=3)
	sumkl = np.sum(klkeep, axis=3)

	# update header to reflect KL modes used
	head["KLMODES"] = str(kllist)

	# write arrays
	fits.writeto(outdir+writestr + '_avgkl.fits', avgkl, head, overwrite=True)
	fits.writeto(outdir+writestr + '_stdevkl.fits', stdevkl, head, overwrite=True)
	fits.writeto(outdir+writestr + '_sumkl.fits', sumkl, head, overwrite=True)
	return (avgkl, stdevkl, sumkl)


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

def find_best_new(pename, kllist, pedir='./', writestr=False, weights=[1,1,1,1,1,1,1,1], outdir='proc/', debug=False, smt=3, snrmeth='all'):
	"""
	collapses parameter explorer file and extracts the optimal parameter value

	REQUIRED INPUTS:
	pename: 	name of paramexplore file
	kllist: 	list of KL modes you want to collapse over
   
	OPTIONAL INPUTS:
	pedir: 		directory holding parameter explorer file
	weights: 	weights for parameter quality metric. Metrics are listed in order below. 
				[peak SNR, peak SNR neighbors, avg SNR, avg SNR neighbors, stdev, stdev neighbors, spurious pixels]
	smt: 		FWHM of gaussian smoothing kernel for 'neighbor quality' metrics


	RETURNS:


	"""
	
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
	nq_stdev = np.zeros([nstepy, nstepx])

	#collapse KL mode cube or extract slice if single
	if len(kllist)>1:
		print("collapsing ", snrmeth, "slices for KL modes", kllist)
	else: 
		print("extracting ", snrmeth, "slices for KL mode", kllist[0])

	#absmed map values
	if snrmeth in ["absmed", "all"]:

		#output cubes from collapsekl have dimensions metric (n=5 for single method, n=9 for 'all'), ann, movm
		avgkl_absmedSNR, stdevkl_absmedSNR, sumkl_absmedSNR = collapsekl(pename, kllist, pedir=pedir, outdir=outdir, snrmeth='absmed', writestr=writestr)

		#finds locations of peaks
		maxind_absmedSNR = np.where(avgkl_absmedSNR[0,:,:] == np.nanmax(avgkl_absmedSNR[0,:,:]))
		maxind_absmedSNR_umask = np.where(avgkl_absmedSNR[1,:,:] == np.nanmax(avgkl_absmedSNR[1,:,:]))

		#translates to x and y coordinates
		xcoord_absmedSNR = maxind_absmedSNR[0][0]
		ycoord_absmedSNR = maxind_absmedSNR[1][0]
		print("peak value for median absolute value SNR is at coordinates:", xcoord_absmedSNR, ycoord_absmedSNR)

		#translates to x and y coordinates
		xcoord_absmedSNR_umask = maxind_absmedSNR_umask[0][0]
		ycoord_absmedSNR_umask = maxind_absmedSNR_umask[1][0]
		print("peak value for median absolute value SNR under mask is at coordinates:", xcoord_absmedSNR, ycoord_absmedSNR)

		#normalize the SNR (where high values = good) maps normally
		snr_norm_absmedSNR = avgkl_absmedSNR[0,:,:] / np.nanmax(avgkl_absmedSNR[0,:,:])
		snr_norm_absmedSNR_umask = avgkl_absmedSNR[1,:,:] / np.nanmax(avgkl_absmedSNR[1,:,:])

		#normalize standard deviations across KL modes. Low values = good
		#slice 0 = peak, slice 1 = average SNR under mask, Going with slice 1 for now
		# divide by SNR so is Stdev in SNR as fraction of SNR itself
		stdev_norm_absmedSNR_cube = stdevkl_absmedSNR[0:2,:,:] / avgkl_absmedSNR[0:2,:,:]
		#first slice is peak
		stdev_norm_absmedSNR = 1 - (stdev_norm_absmedSNR_cube[0,:,:]/np.nanmax(stdev_norm_absmedSNR_cube[0,:,:]))
		#second slice is under mask
		stdev_norm_absmedSNR_umask = 1 - (stdev_norm_absmedSNR_cube[1,:,:]/np.nanmax(stdev_norm_absmedSNR_cube[1,:,:]))

		#spurious pixels metrics - pulling slice 3 (in between IWA and CR only)
		spurpix_absmedSNR = sumkl_absmedSNR[3,:,:]
		
		#make a contrast metric
		#returned quantity is -1*contrast. turn back into contrast
		#log so that higher = better
		logcontrast = np.log10(-1*avgkl_absmedSNR[4,:,:])
		#filter out unphysical contrasts
		logcontrast[logcontrast>0]=np.nan
		#now take absolute value - smaller is better
		logcontrast = np.abs(logcontrast)
		#and subtract the minimum so goes min=0 to max
		logcontrast = logcontrast - np.nanmin(logcontrast)
		#now divide by the max so goes 0-->1
		contrast = logcontrast/np.nanmax(logcontrast)

	#stdev map values
	if snrmeth in ["stdev", "all"]:
		avgkl_stdevSNR, stdevkl_stdevSNR, sumkl_stdevSNR = collapsekl(pename, kllist, pedir=pedir, outdir=outdir, snrmeth='stdev', writestr=writestr)
		
		maxind_stdevSNR = np.where(avgkl_stdevSNR[0,:,:] == np.nanmax(avgkl_stdevSNR[0,:,:]))
		maxind_stdevSNR_umask = np.where(avgkl_stdevSNR[1,:,:] == np.nanmax(avgkl_stdevSNR[1,:,:]))

		xcoord_stdevSNR = maxind_stdevSNR[0][0]
		ycoord_stdevSNR = maxind_stdevSNR[1][0]
		print("peak value for standard deviation SNR is at coordinates:", xcoord_stdevSNR, ycoord_stdevSNR)


		xcoord_stdevSNR_umask = maxind_stdevSNR_umask[0][0]
		ycoord_stdevSNR_umask = maxind_stdevSNR_umask[1][0]
		print("peak value for standard deviation SNR under mask is at coordinates:", xcoord_stdevSNR, ycoord_stdevSNR)

		snr_norm_stdevSNR = avgkl_stdevSNR[0,:,:] / np.nanmax(avgkl_stdevSNR[0,:,:])
		snr_norm_stdevSNR_umask = avgkl_stdevSNR[1,:,:] / np.nanmax(avgkl_stdevSNR[1,:,:])

		stdev_norm_stdevSNR_cube = stdevkl_stdevSNR[0:2,:,:] / avgkl_stdevSNR[0:2,:,:]
		stdev_norm_stdevSNR = 1 - (stdev_norm_stdevSNR_cube[0,:,:]/np.nanmax(stdev_norm_stdevSNR_cube[0,:,:]))
		stdev_norm_stdevSNR_umask = 1 - (stdev_norm_stdevSNR_cube[1,:,:]/np.nanmax(stdev_norm_stdevSNR_cube[1,:,:]))

		spurpix_stdevSNR = sumkl_stdevSNR[3,:,:]

		#make a contrast metric
		#returned quantity is -1*contrast. turn back into contrast
		#log so that higher = better
		logcontrast = np.log10(-1*avgkl_stdevSNR[4,:,:])
		#filter out unphysical contrasts
		logcontrast[logcontrast>0]=np.nan
		#now take absolute value - smaller is better
		logcontrast = np.abs(logcontrast)
		#and subtract the minimum so goes min=0 to max
		logcontrast = logcontrast - np.nanmin(logcontrast)
		#now divide by the max so goes 0-->1
		contrast = logcontrast/np.nanmax(logcontrast)

	
	#spurpix_norm_absmedSNR = 1 - (avgkl_absmedSNR[3,:,:]/np.nanmax(avgkl_absmedSNR[3,:,:]))
	#spurpix_norm_stdevSNR = 1 - (avgkl_stdevSNR[3,:,:]/np.nanmax(avgkl_stdevSNR[3,:,:]))

	#average the two SNR computation methods
	if snrmeth=='all':
		snr_norm_avg = (snr_norm_absmedSNR + snr_norm_stdevSNR) / 2.
		snr_norm_avg_umask = (snr_norm_absmedSNR_umask + snr_norm_stdevSNR_umask) / 2.
		stdev_norm_avg = (stdev_norm_absmedSNR + stdev_norm_stdevSNR) / 2.
		stdev_norm_avg_umask = (stdev_norm_absmedSNR_umask + stdev_norm_stdevSNR_umask) / 2.
		spurpix_avg = (spurpix_absmedSNR + spurpix_stdevSNR) / 2.

	#otherwise, pull just the relevant one
	elif snrmeth=='absmed':
		snr_norm_avg = snr_norm_absmedSNR
		snr_norm_avg_umask = snr_norm_absmedSNR_umask
		stdev_norm_avg = stdev_norm_absmedSNR
		stdev_norm_avg_umask = stdev_norm_absmedSNR_umask
		spurpix_avg = spurpix_absmedSNR
	elif snrmeth=='stdev':
		snr_norm_avg = snr_norm_stdevSNR
		snr_norm_avg_umask = snr_norm_stdevSNR_umask
		stdev_norm_avg = stdev_norm_stdevSNR
		stdev_norm_avg_umask = stdev_norm_stdevSNR_umask
		spurpix_avg = spurpix_stdevSNR
	else:
		print('please provide a valid snr computation method')
		return

	# stdevkl[avgkl<3]=np.nan
	# stdev_norm = 1/(stdevkl/np.nanmin(stdevkl))
	

	#kern = conv.Gaussian2DKernel(x_stddev=2)
	#nq_snr = conv.convolve(snr_norm_avg, kern, preserve_nan=True, nan_treatment='interpolate')
	#nq_stdev = conv.convolve(stdev_norm_avg, kern, preserve_nan = True, nan_treatment='interpolate')
	#nq_snr_umask = conv.convolve(snr_norm_avg_umask, kern, preserve_nan=True, nan_treatment='interpolate')
	#nq_stdev_umask = conv.convolve(stdev_norm_avg_umask, kern, preserve_nan = True, nan_treatment='interpolate')

	#computes neighbor quality by smoothing with Gaussian
	sig=smt
	nq_snr = filter_nan_gaussian_conserving(snr_norm_avg,sig)
	nq_stdev = filter_nan_gaussian_conserving(stdev_norm_avg,sig)
	nq_snr_umask = filter_nan_gaussian_conserving(snr_norm_avg_umask,sig)
	nq_stdev_umask = filter_nan_gaussian_conserving(stdev_norm_avg_umask,sig)

	#normalizes neighbor quality
	nq_snr /= np.nanmax(nq_snr)
	nq_stdev /= np.nanmax(nq_stdev)
	nq_snr_umask /= np.nanmax(nq_snr_umask)
	nq_stdev_umask /= np.nanmax(nq_stdev_umask)

	if debug==True:
		#make a cube of all these metrics for sanity checking
		qual_cube = np.zeros([10,nstepy,nstepx])
		#qual_cube[0,:,:]=snr_norm_absmedSNR
		#qual_cube[1,:,:]=snr_norm_stdevSNR
		qual_cube[0,:,:]=snr_norm_avg
		#qual_cube[3,:,:]=snr_norm_absmedSNR_umask
		#qual_cube[4,:,:]=snr_norm_stdevSNR_umask
		qual_cube[1,:,:]=snr_norm_avg_umask
		#qual_cube[6,:,:]=stdev_norm_absmedSNR
		#qual_cube[7,:,:]=stdev_norm_stdevSNR
		qual_cube[2,:,:]=stdev_norm_avg
		#qual_cube[9,:,:]=stdev_norm_absmedSNR_umask
		#qual_cube[10,:,:]=stdev_norm_stdevSNR_umask
		qual_cube[3,:,:]=stdev_norm_avg_umask
		#qual_cube[12,:,:]=spurpix_norm_absmedSNR
		#qual_cube[13,:,:]=spurpix_norm_stdevSNR
		qual_cube[4,:,:]=spurpix_avg
		qual_cube[5,:,:]=nq_snr
		qual_cube[6,:,:]=nq_snr_umask
		qual_cube[7,:,:]=nq_stdev
		qual_cube[8,:,:]=nq_stdev_umask
		qual_cube[9,:,:]=contrast

		fits.writeto(outdir+pename[:-5]+'_paramqual_cube.fits', qual_cube, overwrite=True)

	#average under mask and peak pixel estimates
	#snr_norm_combo = (snr_norm_avg + snr_norm_avg_umask) / 2.
	#nq_snr_combo = (nq_snr + nq_snr_umask) / 2.
	#stdev_norm_combo = (stdev_norm_avg + stdev_norm_avg_umask) / 2.
	#nq_stdev_combo = (nq_stdev + nq_stdev_umask) / 2.


	#write out the cubes being used for the final metric
	metric_cube = np.zeros([9,nstepy,nstepx])
	metric_cube[0,:,:]= snr_norm_avg
	metric_cube[1,:,:]= nq_snr
	metric_cube[2,:,:]= snr_norm_avg_umask
	metric_cube[3,:,:]= nq_snr_umask
	#going to use the under mask values for the final answer - probably more stable
	metric_cube[4,:,:]= stdev_norm_avg_umask
	metric_cube[5,:,:]= nq_stdev_umask

	#spurious pixel metric = 1 if no spurious pixels and 0 if max number for this dataset
	if np.nanmax(spurpix_avg)>0:
		spurpix_norm_avg = 1-spurpix_avg/np.nanmax(spurpix_avg)
	else: #edge case - no spurious pixels in any image
		spurpix_norm_avg= 1+spurpix_avg
	 
	metric_cube[6,:,:]= spurpix_norm_avg
	metric_cube[7,:,:]= contrast 

	metriclist = (snr_norm_avg, nq_snr, snr_norm_avg_umask, nq_snr_umask, stdev_norm_avg_umask,nq_stdev_umask, spurpix_norm_avg, contrast)
	
	#make sure weights for stdev slices are 0 if only 1 kl mode
	if len(kllist) ==1:
		weights[4:6]=0

	#calculate an aggregate parameter quality metric by summing the individual metrics * their weights
	agg=np.zeros((nstepy,nstepx))
	for metricind in np.arange(len(metriclist)):
		#only sum if non-zero (otherwise nans will carry over)
		if weights[metricind]>0:
			agg+=weights[metricind]*metriclist[metricind]

	metric_cube[8,:,:]=agg

	fits.writeto(outdir+pename[:-5]+'_paramqual_metrics.fits', metric_cube, overwrite=True)

	##find location or peak of parameter quality metric and print info
	ind = np.where(agg == np.nanmax(agg))
	if agg[ind].shape[0]>1:
		print("the optimal solution for this choice of parameters/weights is not unique")
		return()

	#extract metric scores for this location    
	metric_scores = [snr_norm_avg[ind][0], nq_snr[ind][0], snr_norm_avg_umask[ind][0], nq_stdev_umask[ind][0], \
	stdev_norm_avg_umask[ind][0], nq_stdev_umask[ind][0], spurpix_norm_avg[ind][0], agg[ind][0]]

	#translate to annuli and movement values
	ann_val = ymin + ind[0][0] * ystep
	movm_val = xmin + ind[1][0] * xstep
	

	#avgSNR = avgkl_absmedSNR[0][ind] + avgkl_absmedSNR[1][ind]  + avgkl_stdevSNR[0][ind] + avgkl_stdevSNR[1][ind] 
	#avgSNR /= 4

	print('peak is at', ind[0][0], ind[1][0], 'corresponding to annuli', ann_val, ' and movement', movm_val)
	#print('SNR value for fake planets (avg of SNR methods and planets) is', avgSNR)
	print('metric scores for (snr peak, snr peak neigbors, snr umask, snr umask neighbors, stdev, stdev neighbors, spurious pix, contrast, agg) are:', metric_scores)
	return (snr_norm_avg, nq_snr, snr_norm_avg_umask, nq_snr_umask, stdev_norm_avg_umask, nq_stdev_umask, spurpix_norm_avg, contrast, agg, ann_val, movm_val, metric_scores)


def collapse_pes(pedir='./', kllist=[5,10,20,50], wts = [1,1,1,1,1,1,1,1], mode='Line', 
				snrmeth='stdev', snrthresh=False, outdir='proc/', xname='', 
				datadir='../', header=True, smt=3, savefig=True,
				hpval=None, collmode=None, owa=None, oldpe=False, calflux=False):
	"""
	Collapses ALL parameter explorer files in a given directory according to the specified combination of KL modes,
	metric weights, and SNR computation method and runs the corresponding KLIP reductions for that set of parameters

	OPTIONAL INPUTS:
	pedir: 		directory holding parameter explorer files
	kllist:		kl modes to keep. will average over them if more than one is specified.
	weights: 	weights for parameter quality metric. Metrics are listed in order below. 
				[peak SNR, peak SNR neighbors, avg SNR, avg SNR neighbors, stdev, stdev neighbors, spurious pixels]
	mode:		what type of images to run the KLIP reductions on. Options are Line, Cont, or SDI (which will generate
				both and subtract them)
	snrmeth: 	one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value),
				and 'all' to average the two
	snrthresh:	SNR threshhold for peak pixel under mask. All lower values will be masked.
	outdir: 	location to store output files
	xname:		additional string to add to all filenames 
	datadir: 	path to the group's dropbox on your machine (note you can selective sync only relevant dq_cuts subdirectories)
	hpval:		value of the highpass keyword for KLIP. can be True/False or value. 
	collmode:	value for time_collapse keyword for KLIP. Options are 'median', 'mean', 'weighted-mean'
	oldpe:		set to True for older parameter explorers where planet dimension is axis=5 rather than axis=3
	header: 	if True, will ignore specified values for owa, collmode, hpval, and calflux and get them from PE header

	RETURNS:
	d:			a dictionary with all relevant parameter explorer info, KLIP parameters, and KLIPed images


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

		#if header keyword set, pull other KLIP values from header
		if header == True:
			head = fits.getheader(pedir+flist[i])
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

		## runs planet collapse on each PE in the dictionary
		pecube, pcolname = collapse_planets(pename, pedir=pedir, outdir=outdir, snrthresh=snrthresh, oldpe=oldpe)

		#runs PE collapse on planet collapsed cube with specified weights and snrmeth
		#stores optimal annuli and movement values and the aggregate parameter quality map for each PE
		#pedir is outdir here since that's where it was saved by collapse_planets
		j1, j2, j3, j4, j5, j6, j7, j8, agg, ann_val, movm_val, j9 = find_best_new(pcolname, kllist, pedir=outdir, outdir=outdir, writestr=writename, weights=wts, snrmeth=snrmeth, smt=smt)
		d["pe{0}ann".format(i+1)]=ann_val
		d["pe{0}movm".format(i+1)]=movm_val
		d["pe{0}agg".format(i+1)]=agg

		#save visualization of the metrics for this PE explorer
		if savefig==True:
			paramexplore_fig(pcolname, kllist, pedir=outdir, outdir=outdir, writestr=writename, weights=wts, snrmeth=snrmeth, smt=smt)
	
		#define image input direcotries for KLIP based on PE filename
		haindir=d["pe{0}fpath".format(i+1)]+'dq_cuts/'+'Line_'+d["pe{0}cut".format(i+1)]+'_sliced/'
		contindir=d["pe{0}fpath".format(i+1)]+'dq_cuts/'+'Cont_'+d["pe{0}cut".format(i+1)]+'_sliced/'
		
		#set up some unique naming info
		prefix=d["pe{0}pfx".format(i+1)]+xstr
		iwa =d["pe{0}iwa".format(i+1)]
		strklip = '_a'+str(ann_val)+'m'+str(movm_val)+'iwa'+str(iwa)
		d["pe{0}strklip".format(i+1)]=strklip
		
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
			if os.path.exists("{out}/{pre}-KLmodes-all.fits".format(out=outdir, pre=prefix+strklip)):
				print("This file already exists. I am NOT re-running KLIP, but just reading the existing image in. Check and make sure you weren't intending to change the name")
			#otherwise, run KLIP
			else:
				
				dataset = MagAO.MagAOData(filelist) 
				dataset.IWA=iwa
				dataset.OWA=float(owa[i])
				print('running KLIP. highpass = ', hpval[i], ', calflux = ', calflux[i], ', time collapse = ', collmode[i], ', OWA = ', owa[i], 'prefix =', prefix, 'first file =', filelist[0])
				parallelized.klip_dataset(dataset, outputdir=outdir, fileprefix=prefix+strklip, 
                        algo='klip', annuli=ann_val, subsections=1, movement=movm_val,
                        numbasis=kllist, maxnumbasis=100, calibrate_flux=calflux[i], mode="ADI", highpass=float(hpval[i]), 
                        save_aligned=False, time_collapse=collmode[i])
				#fits.writeto('test'+str(i)+'.fits', dataset.output, overwrite=True)

			#pull output image 
			klim = fits.getdata("{out}/{pre}-KLmodes-all.fits".format(out=outdir, pre=prefix+strklip))

			#and store in dictionary
			if runmode=='Cont':
				d["pe{0}contklipim".format(i+1)]=klim
			if runmode=='Line':
				d["pe{0}haklipim".format(i+1)]=klim
	return(d)

def paramexplore_fig(pename, kllist, pedir='proc/', outdir='proc/', writestr=False, weights=[1,1,1,1,1,1,1,1], snrmeth='all', smt=3):
    
    snr_norm_avg, nq_snr, snr_norm_avg_umask, nq_snr_umask, stdev_norm_avg_umask, nq_stdev_umask, spurpix_norm_avg, contrast, agg, ann_val, movm_val, metric_scores = \
        find_best_new(pename, kllist, pedir=pedir, outdir=outdir, writestr=writestr, weights=weights, snrmeth=snrmeth, smt=smt)

    if writestr == False:
        writestr = pename[:-17]

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

    #plt.setp((ax1, ax2, ax3, ax4, ax5), xticks=np.arange(nstepx + 1), yticks=np.arange(nstepy + 1), xticklabels=[], yticklabels=[])

    im1 = ax1.imshow(snr_norm_avg, origin='lower', cmap='magma', vmin=0, vmax=1)
    ax1.set_xlabel("movement parameter")
    ax1.set_ylabel("annuli parameter")
    ax1.set_title("Peak SNR: Weight = "+str(weights[0]))
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im1, cax=cax, orientation='vertical')

    im2 = ax2.imshow(nq_snr, origin='lower',cmap='magma', vmin=0, vmax=1)
    ax2.set_xlabel("movement parameter")
    ax2.set_ylabel("annuli parameter")
    ax2.set_title("Peak SNR Neighbor Quality: Weight = "+str(weights[1]))
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im2, cax=cax, orientation='vertical')

    im3 = ax3.imshow(snr_norm_avg_umask, origin='lower', cmap='magma', vmin=0, vmax=1)
    ax3.set_xlabel("movement parameter")
    ax3.set_ylabel("annuli parameter")
    ax3.set_title("Avg SNR Under Mask: Weight = "+str(weights[2]))
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im3, cax=cax, orientation='vertical')

    im4 = ax4.imshow(nq_snr, origin='lower', cmap='magma', vmin=0, vmax=1)
    ax4.set_xlabel("movement parameter")
    ax4.set_ylabel("annuli parameter")
    ax4.set_title("Avg SNR Neighbor Quality: Weight = "+str(weights[3]))
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im4, cax=cax, orientation='vertical')

    
    im5 = ax5.imshow(stdev_norm_avg_umask, origin='lower', cmap='magma', vmin=0, vmax=1)
    ax5.set_xlabel("movement parameter")
    ax5.set_ylabel("annuli parameter")
    ax5.set_title("Stdev Across KL: Weight = "+str(weights[4]))
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im5, cax=cax, orientation='vertical')

    im6 = ax6.imshow(nq_stdev_umask, origin='lower', cmap='magma', vmin=0, vmax=1)
    ax6.set_xlabel("movement parameter")
    ax6.set_ylabel("annuli parameter")
    ax6.set_title("Stdev Neighbor Quality: Weight = "+str(weights[5]))
    divider = make_axes_locatable(ax6)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im6, cax=cax, orientation='vertical')

    im7 = ax7.imshow(spurpix_norm_avg, origin='lower', cmap='magma', vmin=0, vmax=1)
    ax7.set_xlabel("movement parameter")
    ax7.set_ylabel("annuli parameter")
    ax7.set_title("Spurious Pixels: Weight = "+str(weights[6]))
    divider = make_axes_locatable(ax7)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im7, cax=cax, orientation='vertical')


    im8 = ax8.imshow(contrast, origin='lower', cmap='magma', vmin=0, vmax=1)
    ax8.set_xlabel("movement parameter")
    ax8.set_ylabel("annuli parameter")
    ax8.set_title("Contrast: Weight = "+str(weights[7]))
    divider = make_axes_locatable(ax8)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im8, cax=cax, orientation='vertical')

    # plot metric
    im9 = ax9.imshow(agg, origin='lower', vmin=0, vmax=np.sum(weights))
    ax9.set_ylabel("annuli parameter")
    ax9.set_xlabel("movement parameter")
    ax9.set_title("Aggregate Parameter Quality")
    divider = make_axes_locatable(ax9)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im9, cax=cax, orientation='vertical')

    ind = np.where(agg == np.nanmax(agg))
    label_text = 'a' + str(ann_val) + 'm' + str(movm_val)
    rect = patches.Rectangle((ind[1][0] - 0.5, ind[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='r', facecolor='none')
    ax9.add_patch(rect)
    ax9.text(ind[1][0] + 0.75, ind[0][0], label_text, color='red')

    plt.suptitle(writestr)
    gs.tight_layout(fig, rect=[0, 0.03, 1, 0.95])

    plt.savefig(outdir+writestr+'_paramqual.png')
    
    return(ann_val, movm_val, agg)

def make_klip_snrmaps(d, pedir='./', outdir='proc/', smooth=0.5, snrmeth='stdev', mode='Line', scale=1):
	"""
	Generate SNR maps with specified parameters for all KLIP images in a directory

	REQUIRED INPUTS:
	d: 		The dictionary storing the KLIPed images

	OPTIONAL INPUTS:
	outdir: 	location to store output files
	smooth:		gaussian FWHM with which to smooth the map
	mode:		what type of image to pull and generate SNR maps for
	snrmeth: 	one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value),
				and 'all' to average the two
	scale:		for SDI mode, how to scale the images

	RETURNS:
	d:		Ths input dictionary with SNR map images added

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
	d: 		The dictionary storing the PE info

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
	d: 		The dictionary storing the PE info

	OPTIONAL INPUTS:
	boxsz:	radius of box around image center to show

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
