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


def collapse_planets(pename, pedir='./', outdir='proc/', writestr=False, snrthresh=False, oldpe=False, separate_planets=False):
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
		writestr = pename[:-16]
		print(writestr)

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
	
	fits.writeto(outdir+writename, pecube, pehead, overwrite=True)

	return (pecube, writename, npldim)

def trimkl(pename, kllist, pedir='./', outdir='proc/', writestr=False):
	"""
	Reads in a parameter explorer file and collapses it in the KLmode dimension (axis 3)

	REQUIRED INPUTS:
	pename: 	name of paramexplore file
	kllist: 	list of KL modes you want to collapse over

	OPTIONAL INPUTS:
	pedir: 		directory holding parameter explorer file
	writestr: 	filename prefix for saved file 
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

def find_best_new(pename, kllist, pedir='./', writestr=False, weights=[1,1,1,1,1,1,1,1], outdir='proc/', snrthresh=False,
	oldpe=False, debug=False, smt=3, snrmeth='all',separate_planets=False, separate_kls=False):
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
	snrmeth: 	one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value),
				and 'all' to average the two


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

	#EXTRACT PLANETS OR COLLAPSE
	if separate_planets==False:
		print("COLLAPSING IN PLANET DIMENSION")
	else:
		print("EXTRACTING PLANETS SEPARATELY")
	pecube, writename, npldim = collapse_planets(pename, pedir=pedir, outdir=outdir, snrthresh=snrthresh, oldpe=oldpe, writestr=writestr, separate_planets=separate_planets)

	#EXTRACT KL MODES OR COLLAPSE
	print("EXTRACTING ONLY KL MODES SPECIFIED")
	kltrim, writename = trimkl(writename, kllist, pedir=pedir, outdir=outdir,writestr=writestr)

	if writestr==False:
		writestr=writename[:-5]

	# if collapsing, make mean and stdev arrays
	if separate_kls==False:
		print("COLLAPSING IN KL DIMENSION")
		stdevkl = np.std(kltrim, axis=2, keepdims=True)
		#sumkl = np.sum(kltrim, axis=2, keepdims=True)
		#overwrite kltrim with average
		kltrim= np.mean(kltrim, axis=2, keepdims=True)
		print(stdevkl.shape)

		#grab header
		head = fits.getheader(outdir+writename)
		head["KLCOLL"]='True'

		# write arrays
		fits.writeto(outdir+ writestr + '_avgkl.fits', kltrim, head, overwrite=True)
		fits.writeto(outdir+ writestr + '_stdevkl.fits', stdevkl, head, overwrite=True)
		#fits.writeto(outdir+writestr + '_sumkl.fits', sumkl, head, overwrite=True)

	##EXTRACT SNR SLICES
	#pull the SNR map slices with the matching snrmeth value
	if snrmeth == "absmed":
		slice = 1
	if snrmeth == "stdev":
		slice = 0  

	dims = kltrim.shape

	#extract the appropriate slices
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
		stdev_valid=True
	else: 
		print("EXTRACTING SLICES FOR KL MODES", kllist, 'SEPARATELY')
		klloop = len(kllist)
		stdev_valid=False

	#set up cubes for planet, kl loop
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
			nq_snr = filter_nan_gaussian_conserving(snr_norm,sig)
			nq_stdev = filter_nan_gaussian_conserving(stdev_norm,sig)
			nq_snr_umask = filter_nan_gaussian_conserving(snr_norm_umask,sig)
			nq_stdev_umask = filter_nan_gaussian_conserving(stdev_norm_umask,sig)

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
				weights[4]=0
				weights[5]=0

			#calculate an aggregate parameter quality metric by summing the individual metrics * their weights
			for metricind in np.arange(len(metriclist)):
				#only sum if non-zero (otherwise nans will carry over)
				if weights[metricind]>0:
					agg+=weights[metricind]*metriclist[metricind]

			metric_cube[8,:,k,p,:,:]=agg
			agg_cube[k,p,:,:]=agg

			fits.writeto(outdir+writename[:-5]+'_paramqual_metrics.fits', metric_cube, overwrite=True)

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
			
			if separate_planets==False:
				plno = 'all'
			else:
				plno = p+1
			print('peak for planet =', plno, 'klmode = ', kllist[k], 'is at', ind[0][0], ind[1][0], 'corresponding to annuli', ann_val[k][p], ' and movement', movm_val[k][p])
			#print('SNR value for fake planets (avg of SNR methods and planets) is', avgSNR)
			#print('metric scores for (snr peak, snr peak neigbors, snr umask, snr umask neighbors, stdev, stdev neighbors, spurious pix, contrast, agg) are:', metric_scores)


	if debug==True:
		fits.writeto(outdir+writename[:-5]+'_paramqual_cube.fits', qual_cube, overwrite=True)
	
	return metric_cube, agg_cube, ann_val, movm_val, metric_scores


def collapse_pes(pedir='./', kllist=[5,10,20,50], wts = [1,1,1,1,1,1,1,1], mode='Line', 
				snrmeth='stdev', snrthresh=False, outdir='proc/', xname='', 
				datadir='../', header=True, smt=3, savefig=True,
				hpval=None, collmode=None, owa=None, oldpe=False, calflux=False, 
				separate_planets=False, separate_kls=False):
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

		head = fits.getheader(pedir+flist[i])
		file = fits.gethdata(pedir+flist[i])
		nplanets = file.shape[3]
		print(file, nplanets)

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

		## extract best parameters according to weights and methods

		#colapse in planets dimension (or extract separately)
		pecube, pcolname = collapse_planets(pename, pedir=pedir, outdir=outdir, snrthresh=snrthresh, oldpe=oldpe, separate_planets=separate_planets)
			
		#collapse in kl dimension

		for i in np.arange(nplanets):
			#runs PE collapse for each planet with specified weights and snrmeth
			metric_cube, ann_val_single, movm_val_single, metric_scores= find_best_new(pcolname+str(i)+'.fits', kllist, pedir=outdir, outdir=outdir, writestr=writename, weights=wts, snrmeth=snrmeth, smt=smt, separate_kls=separate_kls)
				
			#define some arrays for storing planet and/or kl slices
			if i == 0:
				shape = [sh for sh in agg_val_single.shape]
				shape.append(nplanets)
				if separate_kls==True:
					shape.append(nkls)
					ann_val = np.zeros(nplanets,nkls)
					movm_val = np.zeros(nplanets,nkls)
				else:
					ann_val = np.zeros(nplanets)
					movm_val = np.zeros(nplanets)
				agg=np.zeros(shape)

			agg[:,:,i]=agg_single
			ann_val.append(ann_val_single)
			movm_val.append(movm_val_single)

		#option 2 - collapse planets and extract signle optimum value
		else:
			pecube, pcolname = collapse_planets(pename, pedir=pedir, outdir=outdir, snrthresh=snrthresh, oldpe=oldpe)

			#runs PE collapse on planet collapsed cube with specified weights and snrmeth
			#stores optimal annuli and movement values and the aggregate parameter quality map for each PE
			#pedir is outdir here since that's where it was saved by collapse_planets
			metric_cube, agg, ann_val, movm_val, metric_scores = find_best_new(pcolname, kllist, pedir=outdir, outdir=outdir, writestr=writename, weights=wts, snrmeth=snrmeth, smt=smt)
			

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
    
    metric_cube, agg, ann_val, movm_val, metric_scores = find_best_new(pename, kllist, pedir=pedir, 
    	outdir=outdir, writestr=writestr, weights=weights, snrmeth=snrmeth, smt=smt)

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

    im1 = ax1.imshow(metric_cube[0,:,k,p,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
    ax1.set_xlabel("movement parameter")
    ax1.set_ylabel("annuli parameter")
    ax1.set_title("Peak SNR: Weight = "+str(weights[0]))
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im1, cax=cax, orientation='vertical')

    im2 = ax2.imshow(metric_cube[1,:,k,p,:,:], origin='lower',cmap='magma', vmin=0, vmax=1)
    ax2.set_xlabel("movement parameter")
    ax2.set_ylabel("annuli parameter")
    ax2.set_title("Peak SNR Neighbor Quality: Weight = "+str(weights[1]))
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im2, cax=cax, orientation='vertical')

    im3 = ax3.imshow(metric_cube[2,:,k,p,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
    ax3.set_xlabel("movement parameter")
    ax3.set_ylabel("annuli parameter")
    ax3.set_title("Avg SNR Under Mask: Weight = "+str(weights[2]))
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im3, cax=cax, orientation='vertical')

    im4 = ax4.imshow(metric_cube[3,:,k,p,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
    ax4.set_xlabel("movement parameter")
    ax4.set_ylabel("annuli parameter")
    ax4.set_title("Avg SNR Neighbor Quality: Weight = "+str(weights[3]))
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im4, cax=cax, orientation='vertical')

    
    im5 = ax5.imshow(metric_cube[4,:,k,p,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
    ax5.set_xlabel("movement parameter")
    ax5.set_ylabel("annuli parameter")
    ax5.set_title("Stdev Across KL: Weight = "+str(weights[4]))
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im5, cax=cax, orientation='vertical')

    im6 = ax6.imshow(metric_cube[5,:,k,p,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
    ax6.set_xlabel("movement parameter")
    ax6.set_ylabel("annuli parameter")
    ax6.set_title("Stdev Neighbor Quality: Weight = "+str(weights[5]))
    divider = make_axes_locatable(ax6)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im6, cax=cax, orientation='vertical')

    im7 = ax7.imshow(metric_cube[6,:,k,p,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
    ax7.set_xlabel("movement parameter")
    ax7.set_ylabel("annuli parameter")
    ax7.set_title("Spurious Pixels: Weight = "+str(weights[6]))
    divider = make_axes_locatable(ax7)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im7, cax=cax, orientation='vertical')


    im8 = ax8.imshow(metric_cube[7,:,k,p,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
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
