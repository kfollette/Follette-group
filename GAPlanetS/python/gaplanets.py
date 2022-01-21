from astropy.io import fits
import os
import pandas as pd
import numpy as np
import astropy.convolution as conv
from astropy.modeling import models, fitting, functional_models
import pyklip_headers as pykh
import glob
import pyklip.klip as klip
import pyklip.instruments.MagAO as MagAO
import pyklip.parallelized as parallelized
import pyklip.fakes as fakes
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import shutil
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import matplotlib.patches as patches
import SNRMap_new as snr
from importlib import reload
from datetime import datetime
from scipy import ndimage
import peproc as pe
import matplotlib.gridspec as gridspec
import datetime as dt
import matplotlib.pylab as pl
import pickle
from matplotlib import rc
import accmodels as am

rc("text", usetex=True)
plt.rcParams.update(plt.rcParamsDefault)

#import accmodels as am

def SliceCube(imfile, rotfile, indir='./', slicedir='sliced/'):
    """
    Slices a 3D cube into individual 2D images numbered sequentially
    dir: directory containing the 3D image
    imfile: name of the 3D image
    rotfile: name of the cube containing the rotoff values needed to rotate images
    slicedir: name of the directory to write the sliced images into
    """

    #read in 3D image, header, and rotoffs
    cube = fits.getdata(imfile)
    head = fits.getheader(imfile)
    rotoffs = fits.getdata(rotfile)

    #record image dimensions
    nims, ydim, xdim = cube.shape

    #check that rotoff and 3rd dimension are same length
    if not nims == len(rotoffs):
        print("number of images and number of rotoffs are not equal")
        return

    #create sliced directory if doesn't already exist
    if not os.path.exists(indir+slicedir):
        os.makedirs(indir+slicedir)

    #loop to save 2D images into the sliced directory
    for z in range(nims):
        single = cube[z,:,:]
        head["rotoff"]=rotoffs[z]
        fits.writeto(indir+slicedir+"sliced_"+str(z+1)+".fits", single, head, overwrite=True)

    print("Done creating individual images for KLIP. These live in the", indir+slicedir, 'directory')


    return

def get_df(dfname):
    # define dataframe if doesn't already exist
    if os.path.exists(dfname):
        df=pd.read_csv(dfname)
        #print("found existing data frame. Reading in.") #with the following contents: \n", df)
    else:
        print("creating new data frame")
        #if none found, make one
        df_cols = ['Dataset', 'wl', 'pctcut', 'nims', 'minpeak', 'medfwhm','rotrange']
        # make a blank pandas dataframe
        df = pd.DataFrame({})
        for name in df_cols:
            df[name] = []
    return(df)

def peak_cut(data_str, wl, rdx_params_dfname='rdx_params.csv', rerun=False,
             debug=False, ghost=False, pctcuts=[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90], electrons=False):
    """
    PURPOSE
    Makes new images and rotoff cubes containing only images that have peaks greater
    than or equal to the percentile threshhold(s) specified with the pctcuts keyword.
    Writes files to disk with appropriate suffixes and slices them for KLIP. Records
    peaks, nimages, etc. in dataframe.

    REQUIRED INPUTS
    data_str = unique identifier for this dataset for naming and table entry
    wl = "Line" or "Cont" depending on which you are reducing
    df = pandas dataframe to store output

    OPTIONAL INPUTS
    ghost = if set to True, will do cuts on ghost peak and write ghost and estimate of star
            peak to header
    pctcuts = a list of percentile threshholds
    electrons = keyword to indicate that units are electrons (just used for sanity check of peak values)

    OUTPUTS
    df = populated pandas dataframes

    written by Kate Follette June 2019
    modified October 2020 by KBF- gaussfitter obsolete. Moving to standard astropy models. Moffatt better than 
    Gaussian. Debug now displays the fit and the image every 10 images in the cube
    """

    #check whether has already been run

    df = get_df(rdx_params_dfname)

    #check whether this dataset has values for 90pct cut (last)
    try: 
        if int(df.loc[(df["Dataset"]==data_str) &  (df["pctcut"]==90),["nims"]].values[0][0]) > 0:
            if rerun==False:
                print("dq cuts have already been run. not rerunning")
                return
    except:
        print('running data quality cuts', pctcuts)

    if wl == 'Line':
        rotstring = 'rotoff_noLinecosmics'
    if wl == 'Cont':
        rotstring = 'rotoff_noContcosmics'

    #read in rotoffs
    rotoffs = fits.getdata('preprocessed/'+rotstring + '.fits')
    
    #figure out what the fname is
    preims = glob.glob('preprocessed/'+wl+'*nocosmics.fits')
    preims = preims[0].split(wl)
    imstring = preims[1][:-5]

    imcube = fits.getdata('preprocessed/'+wl + imstring + '.fits')
    head = fits.getheader('preprocessed/'+wl + imstring + '.fits')

    #create dq_cuts directory if doesn't already exist)
    if not os.path.exists('dq_cuts'):
        os.makedirs('dq_cuts')
    if not os.path.exists('dq_cuts/cubes'):
        os.makedirs('dq_cuts/cubes')


    # Basic parameters needed for fits
    dim = imcube.shape[1]  # dimension
    nims = imcube.shape[0]
    print('found', nims, 'images')
    cen = int((dim - 1) / 2.)
    xcen = cen
    ycen = cen

    #define ghost as center for saturated images
    if ghost == True:
        # will be used as indices so must be integers
        xcen = int(cen + 157.5)
        ycen = int(cen - 7)

    #empty arrays to store fits
    fullfit = []
    peaklist = []
    peaklist2 = []
    fwhmlist = []

    #size of image stamp
    width = 31
    stmpsz = int((width - 1) / 2)

    # loop through images fitting peaks and storing them
    for i in np.arange(0, nims):
        im = imcube[i, ycen - stmpsz-1:ycen + stmpsz, xcen - stmpsz-1:xcen + stmpsz]
        y,x = np.mgrid[:width,:width]
        #g_init=models.Gaussian2D(np.nanmax(im), stmpsz, stmpsz, 2, 2)
        #Gaussian PSF model
        g_init=models.Moffat2D(np.nanmax(im), stmpsz, stmpsz, 6, 1)
        fit_g=fitting.LevMarLSQFitter()
        p=fit_g(g_init,x,y,im)
        fwhm=p.fwhm

        fullfit.append(p)
        peaklist.append(p.amplitude.value)
        fwhmlist.append(fwhm)
        if debug==True:
            if i%10 == 0:
                cmax=np.nanmax(p(x,y))*0.8
                cmin=0

                f, (ax1, ax2, ax3) = plt.subplots(1, 3)
                f1 = ax1.imshow(p(x,y), vmax=cmax, vmin=cmin)
                f.colorbar(f1, ax=ax1, shrink=0.5)

                f2 = ax2.imshow(im, vmax=cmax, vmin=cmin)
                f.colorbar(f2, ax=ax2, shrink=0.5)

                f3 = ax3.imshow(im-p(x,y))
                f.colorbar(f3, ax=ax3, shrink=0.5)
                #f.colorbar(plt3,ax=ax3,)
                plt.show()

        #print('check amplitude:', p.amplitude, '~', p2[0])
        #print('check fwhm:', fwhm, '~', p2[3])

        if electrons == True:
            # if units are electrons, raise the peak threshhold
            if (p.amplitude.value < 0) or (p.amplitude.value > 50000):
                print("warning: unphysical peak value of", p.amplitude.value, 'for image', i + 1)
        if electrons == False:
            if (p.amplitude.value < 0) or (p.amplitude.value > 17000):
                print("warning: unphysical peak value of", p.amplitude.value, 'for image', i + 1)


#    print('median fwhm', np.median(fwhmlist))
#    plt.plot(peaklist, peaklist2, 'gx')
#    plt.xlim(1000,5000)
#    plt.xlabel('new')
#    plt.ylabel('old')
#    plt.ylim(1000,5000)
#    plt.show()

    # empty array to store cuts
    cuts = []

    # find peak values that correspond to those percentiles
    for pct in pctcuts:
        cut = int(np.percentile(peaklist, (pct)))
        print(pct, ' cut = ', cut)
        cuts.append(cut)

    # make the cuts and write out new image and rotoff cubes to directories
    j = 0  # initialize loop counter
    for cut in cuts:
        goodims = []
        goodpeak = []
        goodfwhm = []
        for i in np.arange(0, nims):
            if peaklist[i] > cut:
                goodims.append(i)
                goodpeak.append(peaklist[i])
                goodfwhm.append(fwhmlist[i])
        head["PCTCUT"] = pctcuts[j]
        head["MINPEAK"] = cuts[j]
        if ghost == True:
            head["GHSTPK"] = "True"
        #need fwhm of 0pct cut in all headers
        if j == 0:
            reffwhm = np.nanmedian(goodfwhm)

        #record both fwhm values in headers
        head["0PCTFWHM"]=reffwhm
        head["MEDFWHM"]=np.nanmedian(goodfwhm)

        newim = imcube[goodims, :, :]
        new_nims = newim.shape[0]
        hdu = fits.PrimaryHDU(newim, head)
        imname = 'dq_cuts/cubes/'+ wl + imstring + '_' + str(pctcuts[j]) + 'pctcut.fits'
        hdu.writeto(imname, overwrite=True)
        rotnew = rotoffs[goodims]

        #stupid catch for rot angles that pass through zero
        if rotnew[0] - rotnew[-1] >= 180:
            rotnew[-1]+=360
            rotrange = (rotnew[0] - rotnew[-1]) / ((rotoffs[0]) - (rotoffs[-1]+360)) * 100
        else:
            rotrange = (rotnew[0] - rotnew[-1]) / ((rotoffs[0]) - (rotoffs[-1])) * 100
        hdu2 = fits.PrimaryHDU(rotnew)
        rotname = 'dq_cuts/cubes/' + rotstring + '_' + str(pctcuts[j]) + 'pctcut.fits'
        hdu2.writeto(rotname, overwrite=True)

        if debug == True:
            print(pctcuts[j], "% cut")
            print("Peak > ", cut)
            print("check:", len(goodims) / nims * 100, "pct of images")
            print("shape of new data cube is", newim.shape)
            print(rotrange, "pct of rotation")
            print("median FWHM = ", np.nanmedian(goodfwhm))
            print("slicing cube", imname)

        # slice the newly created file in preparation for KLIP
        outdir = 'dq_cuts/'+ wl + '_' + str(pctcuts[j]) + 'pctcut_sliced/'
        SliceCube(imname, rotname, slicedir=outdir)

        # and add star peak to headers of sliced files
        if ghost == True:
            pykh.addstarpeak(outdir, debug=debug, mask=True, ghost=True, wl=wl)
        else:
            pykh.addstarpeak(outdir, debug=debug, mask=True)

        # record relevant quantities for this data cut to dataframe
        datalist = pd.DataFrame([[data_str, wl, pctcuts[j], new_nims, cuts[j], np.nanmedian(goodfwhm),
                                  rotrange]], columns=['Dataset', 'wl', 'pctcut', 'nims', 'minpeak', 'medfwhm','rotrange'])
        df = df.append(datalist)
        j += 1

    df.to_csv(rdx_params_dfname, index=False)

    return


def get_control_rad():
    filelist = glob.glob('raw/'+'V47*.fits')
    dummyhead = fits.getheader(filelist[10])
    aobin = int(dummyhead["AOBIN"])
    ctrlrads = [30,15] #bin 1, bin 
    ctrlrad = ctrlrads[aobin-1]
    return(ctrlrad)

def make_prefix(data_str, wl, cut):
    prefix = data_str + '_' + wl + '_' + str(cut) + 'pctcut'
    return(prefix)

def compute_thrpt(data_str, wl, cut, outputdir = 'dq_cuts/contrastcurves/', numann=3, movm=4, KLlist=[10], IWA=0,
                  contrast=1e-2, theta=0., clockang=85, debug=False, record_seps=[0.1, 0.25, 0.5, 0.75, 1.0],
                  highpass=True, ghost=False, savefig=False, overwrite=False, iterations=3,rdx_params_dfname='rdx_params.csv',
                  timecoll='median', usecutfwhm=False):
    """
    PURPOSE
    Injects false planets separated by the measured fwhm of the dataset in radius and the clockang parameter
    in azimuth and calculates the  throughput.

    REQUIRED INPUTS
    data_str = unique identifier for this dataset for naming and table entry
    wl = "Line" or "Cont" depending on which you are reducing
    cut = what percentage of images do you want to retain?

    OPTIONAL INPUTS
    IWA = interior mask for KLIP. first planet will be injected at IWA+FWHM
    theta = beginning clock angle
    clockang = azimuthal separation of injected planets
    numann, movm, KLlist = the usual klip annuli, movement and KL mode parameters
    debug = if set to True, prints some sanity checks about injected planet locations, etc.
    savefig = if set to true, saves a throughput figure in the contrastcurve directory
    ghost = when set to True, will inject scaled copy of ghost as fake planet
    record_seps = locations of throughputs to be interpolated and recorded
    df = if set, will read in and add throughputs at record_seps locations


    OUTPUT
    thrpt_out = computed separations (row 1), throughput averages (row 2) and individual throughputs (remaining rows)
    zone_boundaries = locations of boundaries between annular zones, in pixels
    thrpt_table_vals OR dframe = either list of throughput values at recod_seps locations OR a dataframe where those
    throughputs have been recorded

    written by Kate Follette June 2019
    """

    #read in data frame
    df = get_df(rdx_params_dfname)

    # if directory doesn't already exist, create it
    if os.path.exists(outputdir) == False:
        print(outputdir, 'doesn\'t exist yet')
        os.makedirs(outputdir)
        
    prefix = make_prefix(data_str, wl, cut)

    #set up some naming stuff
    if KLlist==[1,2,3,4,5,10,20,50,100]:
      klstr='all'
    else:
      if isinstance(KLlist,int):
        klstr=str(KLlist)
      else:
        klstrlist = [str(kl) for kl in KLlist]
        klstr='_'.join(klstrlist)

    dataset_prefix = prefix +  '_a' + str(numann) + 'm' + str(movm) + 'iwa'+str(IWA) + 'hp'+str(highpass) + '_kl'+klstr

    prefix_fakes = '_initPA'+str(theta)+'_CA'+str(clockang)+'_ctrst'+str(contrast)+'_'+str(iterations)+'FAKES'

    dataset_prefix+=prefix_fakes

    tpt_fname = outputdir + dataset_prefix + '_throughputs.fits'

    # to be added to baseline info
    cols = ["ctrst_fkpl", "tpt ann", "tpt movm", "tpt KL", "tpt IWA", "tpt hp", "tpt theta", "tpt clockang", "tpt iters","ghost", "uniq rdx str"]

    #check whether thisDQ cut has been run
    if True in (df.loc[(df["Dataset"]==data_str) &  (df["pctcut"]==cut),["nims"]].values > 0):
        #find locations where this is true
        idx = df.index[(df["Dataset"]==data_str) &  (df["pctcut"]==cut)].tolist()
        print('this dataset and cut has ', len(idx), 'entries already')

        #if only peakcut columns exist, add the relevant new ones
        if "uniq rdx str" not in df.columns:
            for c in cols:
                df[c]=np.nan
                
        ##are any of them filled in? 
        thesedata_strs = df.loc[(df["Dataset"]==data_str) &  (df["pctcut"]==cut),["uniq rdx str"]].values
        uniq_str_list=[]
        for s in thesedata_strs:
            if isinstance(s[0],str):
                uniq_str_list.append(s[0])

        #if just blank entries (just peak cut data exists)
        if len(uniq_str_list)<1:
            print("But no thrpt/contrast data filled in yet. Filling in.")
            #add the uniq rdx str to the existing line
            df.loc[(df["Dataset"]==data_str) &  (df["pctcut"]==cut),["uniq rdx str"]]=dataset_prefix

        #if at least one string entry exists
        else:
            #check whether this reduction is the same   
            if dataset_prefix in uniq_str_list:
                print('and the uniq rdx str is the same.')
                if overwrite==True:
                    print('Overwriting')
                else:
                    print('will read in existing file')
            #if not the same
            else:
                print('but these are new parameters. copying to a new line')
                print('check', dataset_prefix, uniq_str_list)
                #if more than one entry, just copy the first occurrence 
                if len(idx)>1:
                    idx=idx[0]
                dfnewrow=df.loc[idx]
                #clear old values in copied row
                for c in cols:
                    dfnewrow[c]=np.nan
                #fill in this new str
                dfnewrow["uniq rdx str"]=dataset_prefix
                #add to df
                df=df.append(dfnewrow,ignore_index=True)

    #if this cut has not yet been run for this wl
    else:
        print('dataset does not have a basic entry from peak_cut yet')    
        #if hasn't been run, run this peak cut
        peak_cut(data_str, wl, rdx_params_dfname=rdx_params_dfname, rerun=False,
             debug=False, ghost=ghost, pctcuts=[cut])
        #needs to pull the new df now
        df = get_df(rdx_params_dfname)
        df.loc[(df["Dataset"]==data_str) &  (df["pctcut"]==cut),["uniq rdx str"]]=dataset_prefix
    
    #new values and column names to store in df
    vals = (contrast, numann, movm,  ','.join(map("'{0}'".format, KLlist)), IWA, highpass, theta, clockang, iterations, ghost, dataset_prefix)        
    #add contrast to df
    for i in np.arange((len(cols))):
        df.loc[(df["Dataset"]==data_str) &  (df["pctcut"]==cut) & (df["uniq rdx str"]==dataset_prefix),[cols[i]]]=vals[i]
    
    #df["ctrst_fkpl"]=contrast

    calcthrpt=True
    
    if (os.path.exists(tpt_fname)) and (overwrite == False):
        print("throughput file", tpt_fname, "already exists.")
        calcthrpt=False
        thrpt_cube = fits.getdata(tpt_fname)
        #check whether KL modes same
        thrpt_header = fits.getheader(tpt_fname)
        file_kl = list(map(int, thrpt_header['KLMODES'].split(","))) 
        IWA=thrpt_header['IWA']
        if numann !=1:
            zone_boundaries_h = thrpt_header['ZONEBDRY'].split(',')
            zone_boundaries = [int(bdry) for bdry in zone_boundaries_h]
        #for ann=1 zone boundary is OWA
        else:
            sep_spacing = thrpt_cube[0,0,-1] - thrpt_cube[0,0,-2]
            zone_boundaries = [thrpt_cube[0,0,-1]+sep_spacing]
            
        if file_kl == KLlist:
            print("KLmodes match. reading in values.")
            #seps are same for all KL modes. Just pull from first.
            thrpt_seps = thrpt_cube[0,0,:]
            thrpt_avgs = thrpt_cube[:,1,:]
            thrpts = thrpt_cube[:,2:,:]
        else: 
            print("but KL mdoes don't match.")
            print("KL modes are", file_kl, "for file and", KLlist, "in function call")
            calcthrpt = True

    if calcthrpt==True:

        ###load data
        filelist = glob.glob('dq_cuts/'+wl + '_' + str(cut) + 'pctcut_sliced' + "/sliced*.fits")

        # sorts file list so it's sequential
        filelist.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

        ### pull the values of the star peak from the headers
        starpeak = []
        for i in np.arange(len(filelist)):
            head = fits.getheader(filelist[i])
            starpeak.append(head["STARPEAK"])

        if usecutfwhm==True:
            fwhm=head["MEDFWHM"]
        else:
            #pull fwhm of 0pctcut for consistent spacings
            fwhm=head["0PCTFWHM"]
        # import the dataset. Wait to highpass filter until the KLIP call in the next cell
        dataset = MagAO.MagAOData(filelist, highpass=False)
        # set IWA
        dataset.IWA = IWA

        # divide by starpeak values again to get input images in contrast units
        for i in np.arange(0, len(starpeak)):
            dataset.input[i, :, :] /= starpeak[i]

        # compute number of planets to inject
        first = IWA + fwhm*0.5
        n_planets = (dataset.input.shape[1] / 2. - IWA) / fwhm * 2
        #stay away from very outer part of image
        n_planets = int(n_planets)-2
        if debug == True:
            print('I will inject ', n_planets, 'planets')

        # calculates separations (in pixels) where planets will be injected
        thrpt_seps = []
        sep = int(fwhm)/2
        thissep=first
        for i in np.arange(n_planets):
            thrpt_seps.append(thissep)
            thissep += sep

        if debug == True:
            print("injecting planets at the following radii:", thrpt_seps)

        size = dataset.input.shape
        nims = size[0]
        xcen = int((size[1] - 1) / 2)
        ycen = int((size[2] - 1) / 2)

        #save a version of the input dataset without any fakes
        dataset_copy = np.copy(dataset.input)
        imsz = dataset.input.shape[1]

        thrpts = np.zeros((len(KLlist), iterations, len(thrpt_seps)))

        # full sequence of planet injection, recovery and throughput calculation as many times as the iterations keyword, with
        # starting planet locations clocked by 75 degrees each cycle
        for iter in np.arange(iterations):
            #clock starting planet location by 75deg each cycle
            theta=75.*iter
            
            #make prefixes for output files
            pfx=dataset_prefix + '_set' + str(iter+1) 

            #check whether fake files or throughputs with these parameters have already been calculated
            runfakes=True

            if (os.path.exists(outputdir+pfx+'-KLmodes-all.fits')) and (overwrite==False):
                print("file with prefix", pfx, "already exists")
                #check whether KL modes same
                fakes_header = fits.getheader(outputdir+pfx+'-KLmodes-all.fits')
                #pyklip writes out each KL mode value to separate keyword
                kls = fakes_header["KLMODE*"]
                #make a list and populate it with KL mode values
                fake_kls = []
                for i in np.arange(len(kls)):
                    fake_kls.append(fakes_header["KLMODE"+str(i)])
                if fake_kls == KLlist:
                    print('and KL modes match')
                    runfakes=False
                else:
                    print('but KL modes don\'t match. Regenerating.')
                    runfakes = True

            #special case where pipeline broke in between the two steps - needs whole cube in memory to do throughput calc.
            if runfakes == False and calcthrpt == True:
                runfakes = True
                print("but I can't find a throughput file so I am regenerating it.")

            if runfakes==True:

                #pull the clean input data every time
                dataset.input = np.copy(dataset_copy)

                # replace any nans in image with zeros (edges, padded) for fake injection
                dataset.input[np.isnan(dataset.input) == True] = 0.

                #create copy of dataset.input * contrast for injection
                if ghost==False:
                    to_inject = contrast * dataset.input
                
                #inject planets in raw images
                for sep in thrpt_seps:
                    if ghost == True:
                        fakes.inject_planet(dataset.input, dataset.centers, np.repeat(contrast, nims)[0], dataset.wcs, sep, theta, fwhm=fwhm)
                    else:
                        fakes.inject_planet(dataset.input, dataset.centers, to_inject, dataset.wcs, sep, theta, fwhm=fwhm,
                                        stampsize=imsz)  # , thetas=thetas)
                    theta += clockang

                # put NaNs back
                dataset.input[np.isnan(dataset_copy) == True] = np.nan

                # KLIP dataset with fake planets. Highpass filter here.
                parallelized.klip_dataset(dataset, outputdir=outputdir, fileprefix=pfx, algo='klip', annuli=numann,
                                          subsections=1, movement=movm, numbasis=KLlist, calibrate_flux=False, mode="ADI",
                                          highpass=highpass, save_aligned=False, time_collapse=timecoll, maxnumbasis=100)

            # reset initial theta for recovery loop
            inittheta = 75.*iter

            # create a list that will store throughputs for this iteration
            thrpt_list = np.zeros((len(KLlist),len(thrpt_seps)))


            #calculate a few fixed quantities
            annspacing = (imsz / 2. - IWA) / numann
            zone_boundaries_p = np.arange(1, numann) * annspacing + IWA
            zone_boundaries = [int(bdry) for bdry in zone_boundaries_p]

            nimages = dataset.output.shape[1]
            fake_fluxes = np.zeros((n_planets, nimages))

            #klloop
            klctr=0
            for KL in KLlist:

            #planet recovery loop

                plctr=0
                theta=inittheta
                for sep in thrpt_seps:

                    # thetas for retrieve planet call are measured from +x axis, so should always add 90 to pa for thetas keyword
                    # dataset shape is [KL modes, n images, 1(wl dim), xdim, ydim]
                    fake_flux = fakes.retrieve_planet_flux(dataset.output[klctr, :, 0, :, :], dataset.centers, dataset.wcs, sep, theta,
                                                           searchrad=8, guessfwhm=fwhm,
                                                           guesspeak=contrast * np.median(dataset.star_flux),
                                                           thetas=np.repeat(theta + 90, dataset.output.shape[1]), refinefit=True)
                    #record fake flux for this planet for each image in cube
                    fake_fluxes[plctr,:] = fake_flux
                    #translate to throughput (images come out in contrast units)
                    newthpt = np.nanmedian(fake_flux) / (contrast)
                    if (newthpt > 1) or (newthpt < 0):
                        print('invalid throughput value of', newthpt, 'for planet at separation', sep, '. Replacing with NaN')
                        newthpt = np.nan
                   
                    #check for bad value - throughput should increase with distance. decrease suggests inner point may not have been a true recovery
                    if len(thrpt_list) > 1:
                        if newthpt < thrpt_list[klctr,-1]:
                            if abs(newthpt - thrpt_list[klctr,-1]) > 0.1:
                                print('suspicious trend in throughput. Value at', sep, 'is', newthpt, 'but was', thrpt_list[klctr,-1], 'in previous iteration.')
                                #if planet is within 10 pixels of zone boundary, relax the requirement that the trend has to be upward
                                if len(zone_boundaries) > 0:
                                    if min(abs(zone_boundaries - sep)) < 10:
                                        print(min(abs(zone_boundaries - sep)), 'pixels from zone boundary. Keeping.')
                                    else:
                                        print('planet is', min(abs(zone_boundaries - sep)), 'pixels from nearest zone boundary. Setting to NaN.')
                                        newthpt=np.nan

                    thrpt_list[klctr,plctr]=newthpt
                    
                    #if debug == True:
                        #print("fake planet at ", sep, " pixels has a mean throughput of", np.nanmean(fake_flux / (contrast)),
                          #"median of", np.nanmedian(fake_flux / (contrast)), " and stdev of ",
                          #np.nanstd(fake_flux / (contrast)))
                    theta += clockang
                    plctr += 1

                thrpts[klctr,iter,:] = thrpt_list[klctr,:]
                
                if debug==True:
                    print('throughputs for iteration', iter, 'KL mode', KL, 'are', thrpt_list[klctr,:])
                klctr += 1
                #print('check: dataset output dims are', dataset.output.shape)

    #reset KL mode counter
    klctr=0
    thrpt_avgs = np.zeros((len(KLlist),len(thrpt_seps)))

    for KL in KLlist:
        #average iterations
        thrpt_avgs[klctr,:]=np.nanmean(thrpts[klctr,:,:],axis=0)

        #up to 10 iterations. will break if do more. 
        cx = ['rx','gx','yx','cx','mx','rs','gs','ys','cs','ms']
        #if on last iteration, make plot and save
        ctrl_rad = get_control_rad() 

        if (savefig == True):
            #plot the individual points
            plt.figure(figsize=(7,4), dpi=750)
            for iter in np.arange(iterations):
                plt.plot(thrpt_seps, thrpts[klctr,iter,:], cx[iter], label="set"+str(iter+1))
            # plot the throughput averages (should all be <1 and should increase outward until they hit a zone boundary)
            plt.plot(thrpt_seps, thrpt_avgs[klctr,:], 'bo', label="average")
            plt.plot((IWA, IWA),(0,1),'k-', label="IWA")
            plt.plot((ctrl_rad, ctrl_rad),(0,1),'m--', label="control radius")

            for bd in zone_boundaries:
                if bd == zone_boundaries[0]:
                    plt.plot((bd, bd), (0, 1), '--',color='grey',label='zone boundary')
                else:
                    plt.plot((bd, bd), (0, 1), '--', color='grey',)
            plt.xlabel("separation in pixels")
            plt.ylabel("throughput")
            plt.title(dataset_prefix.replace('_',' ')+ ' KL '+str(KL)+ ' Throughput')
            plt.legend()
            plt.savefig(outputdir + dataset_prefix + '_KL'+str(KL)+ '_throughput.jpg')
            plt.show()
            plt.clf()

            # locations to record throughputs in table
            platescale = 0.00795

            # blank array to store interpolated throughputs
            thrpt_table_vals = []
            # loop through locations, find throughput and record
            for loc in record_seps:
                # interpolate over throughputs to find value at loc
                if calcthrpt==True:
                    imsz=dataset.input.shape[1]
                else:
                    imsz=int(thrpt_header["IMSZ"])
                if loc < imsz*platescale:
                    tpt = np.interp(loc / platescale, thrpt_seps, thrpt_avgs[klctr,:])
                else:
                    tpt = np.nan
                thrpt_table_vals.append(tpt)
                # if df keyword is set, records in dataframe.
                if len(df) > 0:
                    colname = "tpt_" + str(loc)+'_KL'+str(KL)
                    # if column doesn't already exist, create and fill with nans
                    if colname not in df:
                        df[colname] = np.nan
                        print("creating column", colname)
                    df.loc[(df["Dataset"]==data_str) &  (df["pctcut"]==cut) & (df["uniq rdx str"]==dataset_prefix),[colname]]=tpt
                    #df[colname].loc[(df.Dataset == data_str) & (df.pctcut == cut)] = tpt
            klctr+=1

    thrpt_out=np.zeros((len(KLlist), 2+iterations, len(thrpt_seps)))
    thrpt_out[:,0,:]=thrpt_seps
    thrpt_out[:,1,:]=thrpt_avgs
    thrpt_out[:,2:,:]=thrpts

    #if this is the first time generating it, write out
    if calcthrpt==True:

        head["NPLANET"]=n_planets
        head["CTRST"]=contrast
        head["ITERS"]=iterations
        head["PASTART"]=0
        head["CLOCKANG"]=clockang
        head["KLMODES"]=str(KLlist)[1:-1]
        head["NUMANN"]=numann
        head["MOVM"]=movm
        head["IWA"]=IWA
        head["ZONEBDRY"]=str(zone_boundaries)[1:-1]
        head["ROW1"]="separations (pix)"
        head["ROW2"]="average throughput"
        head["OTHROWS"]="individual throughputs"
        head["IMSZ"]=dataset.input.shape[1]
                    
        fits.writeto(tpt_fname, thrpt_out, header=head, overwrite=True)

        df.to_csv(rdx_params_dfname, index=False)

    return (thrpt_out, zone_boundaries, df, dataset_prefix)


def make_contrast_curve(data_str, wl, cut, thrpt_out, dataset_prefix, outputdir = 'dq_cuts/contrastcurves/', numann=3,
                        movm=4, KLlist=[10], IWA=0, rdx_params_dfname='rdx_params.csv', record_seps=[0.1, 0.25, 0.5, 0.75, 1.0], 
                        savefig=False, debug=False, overwrite=False, highpass=True, maxsep=1.5, timecoll='median'):
    """
    PURPOSE
    calculates raw contrast by running KLIP on images and then corrects for throughput
    to generate a calibrated contrast curve.

    REQUIRED INPUTS
    data_str = unique identifier for this dataset for naming and table entry
    wl = "Line" or "Cont" depending on which you are reducing
    cut = what percentage of images do you want to retain?
    thrpt_out = output of compute_thrpt.

    OPTIONAL INPUTS
    IWA = interior mask for KLIP. first planet will be injected at IWA+FWHM
    theta = beginning clock angle
    clockang = azimuthal separation of injected planets
    numann, movm, KLlist = the usual klip annuli, movement and KL mode parameters
    debug = if set to True, prints some sanity checks about injected planet locations, etc.
    savefig = if set to true, saves a throughput figure in the contrastcurve directory

    OUTPUT
    contrast_seps = separations in pixels where contrasts were computed
    corrected_contrast_curve = computed (throughput corrected) contrasts at those searations
    IWA = inner working angle

    written by Kate Follette June 2019
    """

    #read in or make data frame
    df = get_df(rdx_params_dfname)

    iterations = len(thrpt_out[0,:,0])-2
    platescale = 0.00795

    # set up directories and naming
    
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    #pull throughput info 
    thrpt_seps = thrpt_out[0,0,:]
    thrpt_avgs = thrpt_out[:,1,:]
    thrpts = thrpt_out[:,2:,:]

    if (os.path.exists(outputdir+dataset_prefix+'-KLmodes-all.fits')) and (overwrite == False):
        klcube = fits.getdata(outputdir+dataset_prefix+'-KLmodes-all.fits')
        klheader = fits.getheader(outputdir+dataset_prefix+'-KLmodes-all.fits')
        print('KLIPed image without fakes already exists for these parameters')

    else:
        # specify directory and read in data again. This time we will NOT inject fakes
        filelist = glob.glob('dq_cuts/'+wl + '_' + str(cut) + 'pctcut_sliced' + "/sliced*.fits")
        filelist.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
        ## create dataset object for processing
        dataset = MagAO.MagAOData(filelist, highpass=False)
        ## set inner working angle to match previous
        dataset.IWA = IWA

        ### pull the values of the star peak from the headers
        starpeak = []
        for i in np.arange(len(filelist)):
            head = fits.getheader(filelist[i])
            starpeak.append(head["STARPEAK"])

        # get image units in terms of contrast rather than counts by dividing each image by the star peak
        for i in np.arange(0, len(starpeak)):
            dataset.input[i, :, :] /= starpeak[i]

        # klip the dataset with same set of KLIP parameters as fake planets
        parallelized.klip_dataset(dataset, outputdir=outputdir, fileprefix=dataset_prefix, annuli=numann, subsections=1,
                                  algo='klip', movement=movm, numbasis=KLlist, calibrate_flux=False,
                                  mode="ADI", highpass=highpass, time_collapse=timecoll, maxnumbasis=100)

        ##pull some needed info from headers
        kl_hdulist = fits.open("{out}/{pre}-KLmodes-all.fits".format(out=outputdir, pre=dataset_prefix))
        klcube = kl_hdulist[1].data
        klheader = kl_hdulist[0].header
    
    #read in fwhm
    dataset_fwhm = klheader["0PCTFWHM"]

    #set up some naming stuff
    if KLlist==[1,2,3,4,5,10,20,50,100]:
      klstr='all'
    else:
      if isinstance(KLlist,int):
        klstr=str(KLlist)
      else:
        klstrlist = [str(kl) for kl in KLlist]
        klstr='_'.join(klstrlist)

    #raw contrast is independent of fake injection so simpler name string
    rawc_prefix = make_prefix(data_str, wl, cut)
    rawc_prefix += '_a' + str(numann) + 'm' + str(movm) + 'iwa'+str(IWA) + 'hp' + str(highpass)+'_kl'+klstr

    #set OWA
    klim = klcube[0, :, :]

    #thrpt doesn't work for planets near the edge. truncating to match those seps
    OWA = klim.shape[1]/2-2*int(dataset_fwhm)

    #check whether has already been computed
    if (os.path.exists(outputdir+rawc_prefix+'_rawcontrast.fits')) and (overwrite==False):
        ctrst = fits.getdata(outputdir+rawc_prefix+'_rawcontrast.fits')
        contrast_seps = ctrst[0,:]
        contrast = ctrst[1,:]
    
    #if not, generate raw contrast and write out as .jpg and .fits
    else:
        dataset_center = [klheader['PSFCENTX'], klheader['PSFCENTY']]
        
        klctr=0
        
        #loop over KL modes
        for KL in KLlist:

            contrast_seps, contrast = klip.meas_contrast(klcube[klctr,:,:], IWA, OWA, dataset_fwhm*2,
                                                     center=dataset_center, low_pass_filter=False)
            if klctr == 0:
                contrast_out=np.zeros((len(KLlist), 3+iterations, len(contrast_seps)))
            
            contrast_out[klctr,0,:]=contrast_seps
            contrast_out[klctr,1,:]=contrast

            ctrl_rad = get_control_rad()
            ctrl_rad = ctrl_rad * platescale

            imsz = klim.shape[1]
            annspacing = (imsz / 2. - IWA) / numann
            zone_boundaries = np.arange(1, numann) * annspacing + IWA

            if savefig == True and debug==True:

                plt.figure(figsize=(7,4), dpi=750)
                plt.plot(contrast_seps * platescale, contrast)
                plt.plot(contrast_seps * platescale, contrast, 'bo')
                plt.yscale("log")
                plt.ylim(np.nanmin(contrast), 1e-1)
                plt.xlim(0, maxsep)
                plt.xlabel("distance in arcseconds")
                plt.ylabel("contrast")
                if IWA > 0:
                    plt.plot((IWA * platescale, IWA * platescale), (1e-5, 1e-1), 'k-', label='IWA')
                plt.plot((ctrl_rad, ctrl_rad),(0,1),'m--', label="control radius")
                
                for bd in zone_boundaries * platescale:
                    if bd < OWA * platescale:
                        if bd == zone_boundaries[0]*platescale:
                            plt.plot((bd, bd), (0, 1), '--', color='grey', label='zone boundary')
                        else:
                            plt.plot((bd, bd), (0, 1), '--', color='grey')
                plt.legend()
                plt.title(rawc_prefix.replace('_',' ')+ "KL"+str(KL)+ " Raw Contrast")
                plt.savefig(outputdir + rawc_prefix + "_KL"+str(KL)+ '_rawcontrast.jpg')
                plt.show()
                plt.clf()  # clear figure


    ## corrected contrast figure
            corrected_contrast_curve = np.copy(contrast)
            interp_thrpts = []

            #interpolate average throughputs at the same separations as the contrasts
            for i, sep in enumerate(contrast_seps):
                thrpt_interp = np.interp(sep, thrpt_seps, thrpt_avgs[klctr,:])
                if debug == True:
                    closest_throughput_index = np.argmin(np.abs(thrpt_seps - sep))
                    print('for separation', sep, " closest throughput is at separation ", thrpt_seps[closest_throughput_index])
                    print('interpolated throughput is', thrpt_interp, 'for separation', sep)
                corrected_contrast_curve[i] /= thrpt_interp
                interp_thrpts.append(thrpt_interp)

            ctrsts = np.zeros((iterations, len(contrast_seps)))
            #compute contrast curves independently for each set of throughputs (use for range)
            for iter in np.arange(iterations):
                thrpts_thisset = np.interp(contrast_seps, thrpt_seps, thrpts[klctr,iter,:])
                ctrst_thisset = contrast/thrpts_thisset
                ctrsts[iter,:]=ctrst_thisset

            if debug==True:
                #check throughput interpolations
                plt.plot(contrast_seps, interp_thrpts, 'k-', label="interpolated")
                plt.plot(thrpt_seps, thrpt_avgs[klctr,:], 'r--', label = "measured")
                plt.title("checking interpolated throughputs")
                plt.legend()
                plt.show()
                plt.clf()

            if savefig == True:
                plt.figure(figsize=(7,4), dpi=750)
                cx = ['r-','g-','y-','c-','m-','r:','g:','y:','c:','m:']
                #plot the individual points
                for iter in np.arange(iterations):
                    thrpts_thisset = np.interp(contrast_seps, thrpt_seps, thrpts[klctr,iter,:])
                    plt.plot(contrast_seps * platescale, contrast/thrpts_thisset, cx[iter], label="set"+str(iter+1))
                plt.plot(contrast_seps * platescale, corrected_contrast_curve, label=r'average 5$\sigma$ contrast')
                plt.plot(contrast_seps * platescale, contrast, label=r'raw 5$\sigma$ contrast', color='gray')
                plt.yscale("log")
                plt.ylim(np.nanmin(contrast), 1e-1)
                plt.xlabel("distance in arcseconds")
                plt.ylabel("contrast")
                plt.xlim(0,maxsep)
                if IWA > 0:
                    plt.plot((IWA * platescale, IWA * platescale), (1e-5, 1e-1), 'k-', label='IWA')
                plt.plot((ctrl_rad, ctrl_rad),(0,1),'m--', label="control radius")
                if numann!=1:
                    for bd in zone_boundaries * platescale:
                        if bd < OWA * platescale:
                            if bd == zone_boundaries[0]*platescale:
                                plt.plot((bd, bd), (0, 1), '--', color='grey', label='zone boundary')
                            else:
                                plt.plot((bd, bd), (0, 1), '--', color='grey')
                plt.legend()
                plt.title(rawc_prefix.replace('_',' ') + " KL"+str(KL)+" Corrected Contrast")
                plt.savefig(outputdir + dataset_prefix + "_KL"+str(KL)+ '_contrastcurve.jpg')
                plt.show()
                plt.clf()  # clear figure

                contrast_out[klctr,2,:]=corrected_contrast_curve
                contrast_out[klctr,3:,:]=ctrsts
            klctr+=1

    fits.writeto(outputdir + dataset_prefix +  '_contrasts.fits', contrast_out, overwrite=True)
    

    # locations to record throughputs in table
    # ctrst_table_loc = np.array(record_seps)/platescale
    # blank array to store interpolated throughputs
    ctrst_table_vals = []
    # loop through locations, find throughput and record
    for loc in record_seps:
        # interpolate over throughputs to find value at loc
        contrast = np.interp(loc / platescale, contrast_seps, corrected_contrast_curve)
        ctrst_table_vals.append(contrast)
        # if df keyword is set, recors in dataframe.
        if len(df) > 0:
            colname = "ctrst_" + str(loc) +'_KL'+str(KL)
            # if column doesn't already exist, create and fill with nans
            if colname not in df:
                df[colname] = np.nan
                print("creating column", colname)

            df.loc[(df["Dataset"]==data_str) &  (df["pctcut"]==cut) & (df["uniq rdx str"]==dataset_prefix),[colname]]=contrast       
            #df[colname].loc[(df.Dataset == data_str) & (df.pctcut == cut)] = contrast
    
    df.to_csv(rdx_params_dfname, index=False)

    return (contrast_out, df, OWA)

def cut_comparison(data_str, wl, outputdir='dq_cuts/contrastcurves/',pctcuts=[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90], record_seps=[0.1, 0.25, 0.5, 0.75, 1.0],
                   contrast=1e-2, numann=3, movm=4, KLlist=[10], IWA=0, savefig=False, ghost=False, rdx_params_dfname='rdx_params.csv', debug=False,
                   iterations=3, theta=0., overwrite=False, highpass=True, clockang=85, usecutfwhm=False):
    """
    PURPOSE
    loop through data quality cuts and compile corrected contrast curves into single array

    REQUIRED INPUTS
    data_str = unique identifier for this dataset for naming and table entry
    wl = "Line" or "Cont" depending on which you are reducing
    fwhm = fwhm of the dataset (we use fwhm of median determined with idl atv aperture photometry)

    OPTIONAL INPUTS
    IWA = interior mask for KLIP. first planet will be injected at IWA+FWHM
    numann, movm, KLlist = the usual klip annuli, movement and KL mode parameters
    savefig = if set to true, saves a throughput figure in the contrastcurve directory

    OUTPUT
    contrast_seps = separations in pixels where contrasts were computed
    contrasts = ncuts x contrast_seps array storing all corrected curves
    zone_boundaries = locations of boundaries between annular zones, in pixels
    IWA = inner working angle

    written by Kate Follette June 2019
    2020 - modified to accommodate multiple KL modes
    """
    #read in or make data frame
    df = get_df(rdx_params_dfname)

    fakestr = '_initPA'+str(theta)+'_CA'+str(clockang)+'_ctrst'+str(contrast)+'_'+str(iterations)+'FAKES'
    klipstr =  '_a' + str(numann) + 'm' + str(movm) + 'iwa'+str(IWA) +'hp'+str(highpass)#+'_KL'+str(KLlist[0])
    outstr = klipstr + fakestr

    i=0

    ##loop through data quality cuts
    for cut in pctcuts:

        ##add logic to find existing files if they've already been computed!
        prefix = make_prefix(data_str, wl, cut)
        namestr = prefix+outstr

        if (os.path.exists(outputdir + namestr + '_throughputs.fits')) and (overwrite==False):
            print ('found existing throughput file', outputdir + namestr + '_throughputs.fits')
            thrpt_out = fits.getdata(outputdir + namestr + '_throughputs.fits')
            head = fits.getheader(outputdir + namestr + '_throughputs.fits')
            zone_boundaries = list(map(int, head['ZONEBDRY'].split(",")))

        else:
            print('computing throughputs for', cut, 'pct cut')
            thrpt_out, zone_boundaries, df, dataset_prefix = compute_thrpt(data_str, wl, cut,
                                                                    savefig=savefig, ghost=ghost, contrast=contrast,
                                                                    record_seps=record_seps, theta=theta,
                                                                    outputdir=outputdir, clockang=clockang,
                                                                    #KLIP parameters
                                                                    numann=numann, movm=movm, highpass=highpass,
                                                                    KLlist=KLlist, IWA=IWA, rdx_params_dfname=rdx_params_dfname, 
                                                                    debug=debug, iterations=iterations, overwrite=overwrite,
                                                                    usecutfwhm=usecutfwhm)

        if (os.path.exists(outputdir + namestr + '_contrasts.fits')) and (overwrite==False):
            print ('found existing contrast curve', outputdir + namestr + '_contrasts.fits')
            ctrsts = fits.getdata(outputdir + namestr + '_contrasts.fits')

        else:
            print('computing contrasts for', cut, 'pct cut')
            ctrsts, df, OWA = make_contrast_curve(data_str, wl, cut,
                                                                 thrpt_out, namestr, record_seps=record_seps,
                                                                 savefig=savefig, outputdir=outputdir, overwrite=overwrite,
                                                                 ##KLIP parameters
                                                                 numann=numann, movm=movm, KLlist=KLlist, IWA=IWA, highpass=highpass,
                                                                 rdx_params_dfname=rdx_params_dfname, debug=debug)

        # compile contrasts for all cuts into array
        contrast_seps = ctrsts[0,0,:]

        #note assumes here just one kl mode in slice with index 0
        if cut == 0:
            # create a blank array to store contrast curves
            contrasts = np.zeros([len(pctcuts), len(KLlist), len(ctrsts[0,2,:])])
        
        klctr=0
        for kl in KLlist:
            #pull only throughput corrected average curve (slice 2 in ctrsts)
            contrasts[i,klctr,:] = ctrsts[klctr,2,:]
            klctr+=1

        # loop counter
        i += 1

    return (contrast_seps, contrasts, zone_boundaries, IWA, df, KLlist, outstr)

def split_objstr(objstr, fakestr=False):

    #prefix = data_str + '_' + wl + '_' + str(cut) + 'pctcut'

    #dataset_prefix = prefix +  '_a' + str(numann) + 'm' + str(movm) + 'iwa'+str(IWA) + 'hp'+str(highpass) + klstr

    #prefix_fakes = '_initPA'+str(theta)+'_CA'+str(clockang)+'_ctrst'+str(contrast)+'_'+str(iterations)+'FAKES'

    #dataset_prefix+=prefix_fakes

    print(objstr)
    #klsplit = objstr.split('kl')
    #klmodelist = klsplit[-1].split('_')
    #kllist = []
    #for k in kllist:
        #try:
            #kllist.append(int(k))
        #except:
            #print('end of kl mode string')

    klipstr = re.split('_a|m|iwa|hp|_', objstr)

    ann = klipstr[1]
    movm = klipstr[2]
    iwa = klipstr[3]
    hp = klipstr[4]

    if fakestr==True:
        fakestr = objstr.split('_initPA')
        fakestr_list = re.split('_CA|_ctrst|_|FAKES',fakestr[1])
        pa = fakestr_list[0]
        ca = fakestr_list[1]
        ctrst = fakestr_list[2]
        numit = fakestr_list[3]
        fakelist=(pa,ca,ctrst,numit)
    else:
        fakelist=()
    return(ann, movm, iwa, hp, fakelist)


def contrastcut_fig(data_str, wl, contrast_seps, contrasts, zone_boundaries, KLlist, outstr, outputdir = 'dq_cuts/contrastcurves/', outer=50, IWA=0,
                    pctcuts=[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90], chosen=None):
    """
    PURPOSE
    Generate figure for comparison of contrasts with different image quality cuts

    REQUIRED INPUT
    data_str = unique identifier for this dataset for naming and table entry
    wl = "Line" or "Cont" depending on which you are reducing
    contrast_seps = separations in pixels where contrasts were computed
    contrasts = ncuts x contrast_seps array storing all corrected curves
    zone_boundaries = locations of boundaries between annular zones, in pixels

    OPTIONAL INPUT
    IWA = inner working angle
    pct cuts = pct cut labels

    OUTPUT
    figure showing contrast curve comparison for different image quality cuts

    written by Kate Follette June 2019
    """

    platescale = 0.00795
    ctrl_rad = get_control_rad()
    ctrl_rad*=platescale

    outer_asec = outer*platescale
    #what is index of last contrast within specified zone
    last_idx=-1
    for dist in contrast_seps:
        if dist<outer:
            last_idx+=1

    if 'FAKES' in outstr:
        fakestr=True
    else:
        fakestr=False

    ann, movm, iwa, hpstr, fakelist = split_objstr(outstr, fakestr=fakestr)

    text_kwargs = dict(ha='center', va='center', fontsize=28, color='C1')

    klctr=0
    for kl in KLlist:
        # plt.style.use('seaborn-colorblind')
        f, ax2 = plt.subplots(figsize=(7,4), dpi=750)
        ax1 = ax2.twinx()  

        extrap_ctrst=[] 
        for i in np.arange(len(pctcuts)):
            cut = pctcuts[i]
            j = i / len(pctcuts)
            if i % 2 == 0:
                linesty = '-'
            else:
                linesty = '--'

            if chosen!=None and cut==int(chosen):
                clr = 'k'
                linewid=2
                linesty = '-'
                for_oplot = contrasts[i,klctr, :]
            else:
                clr = cm.plasma(j)
                linewid=1
            
            ax1.plot(contrast_seps * platescale, contrasts[i,klctr, :],
                     label=str(cut) + r'\% cut', linestyle=linesty, linewidth = linewid, 
                     color=clr)

            #extrapolate to iwa
            #slope of innermost two points
            m = (np.log10(contrasts[i,klctr,0])-np.log10(contrasts[i,klctr,1]))/((contrast_seps[1]-contrast_seps[0])*platescale)
            IWA_y = np.log10(contrasts[i,klctr,0]) + m*(contrast_seps[0] - IWA)*platescale
            extrap_ctrst.append(10**IWA_y)
            if chosen!=None and cut==int(chosen):
                IWA_y_oplot=IWA_y
            ax1.plot([IWA * platescale, contrast_seps[0]*platescale],[10**IWA_y, contrasts[i,klctr,0]],  linesty, linewidth=linewid, color=clr, alpha=0.3)

        if chosen!=None:
            #overplot chosen cut so on top
            ax1.plot(contrast_seps * platescale, for_oplot, linestyle='-', linewidth = 2, color='k')
            ax1.plot([IWA * platescale, contrast_seps[0]*platescale],[10**IWA_y_oplot, for_oplot[0]],  '-', linewidth=2, color='k', alpha=0.3)
 
        ax1.yaxis.tick_left()
        ax1.yaxis.set_label_position("left")
 
        cfloor = np.nanmin(contrasts[:,klctr,:])
        cmax=np.nanmax(extrap_ctrst)

        if IWA > 0:
            ax2.axvline(x=(IWA * platescale), linestyle='--', color='k', label='IWA')

        ax2.axvline(x=ctrl_rad,linestyle='--', color='m', label="control radius")
        ax2.get_yaxis().set_visible(False)
        
        for bd in np.multiply(zone_boundaries,platescale):
            if bd < outer_asec:
                if bd == zone_boundaries[0]*platescale:
                    ax2.axvline(x=bd, linestyle='--', color='grey', label='zone boundary')
                else:
                    ax2.axvline(x=bd, linestyle='--', color='grey')

        plt.yscale("log")

        dset_label = data_str.split('_')
        #translate to mm/dd/yy format
        monthstrs = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
        objname = dset_label[0]
        datestr = dset_label[1]
        print(objname, datestr)
        monthstr = [s for s in monthstrs if s in datestr]
        monthno = [i for i in np.arange(1,13) if monthstrs[i-1] in datestr]
        year = int('20'+datestr[-2:])
        day = datestr.split(monthstr[0])[0]
        #print(line_date, '= day:', day, ' month:', monthno[0], ' year:', year)
        date_label = dt.date(year=year, month=int(monthno[0]), day=int(day))


        plt.title(objname + ' ' + date_label.strftime("%m/%d/%Y"))
        plt.xlim(0, outer_asec)
        plt.autoscale(enable=True, axis='y', tight=True)
        #plt.ylim(10 ** cfloor, cmax)
        ax2.set_xlabel("distance in arcseconds")
        ax1.set_ylabel("contrast")

        #add some text with reduction parameters
        
        text_block = r'\textbf{KLIP Parameters:}'+ '\n annuli = ' + ann + '\n movement = ' + movm
        text_block+= '\n IWA = ' + iwa + '\nhighpass = ' + hpstr + ' pix'
        if len(fakelist)>1:
            pa, ca, ctrst, numit = fakelist
            text_block+= '\n\n'+r'\textbf{False Planet Parameters}:' +'\n'
            text_block+= 'initial PA = ' + pa + '\n clock angle =' + ca + '\n n iterations =' + numit 
            text_block+= '\ncontrast = '+ ctrst

        ax1.text(outer_asec*0.8, (10**(np.log10(cmax)-(np.log10(cmax)-np.log10(cfloor))/2)), text_block, multialignment="left",size=6)

        #separate legend for cut lines and IWA/CR/ZB
        ax2.legend(loc='upper right', fontsize='x-small', facecolor='white', framealpha=1)
        ax1.legend(loc='lower left', fontsize='x-small', facecolor='white', framealpha=1)

        # write out 
        plt.savefig(outputdir + data_str + outstr + '_KL'+str(kl)+'_contrastsbycut.jpg')
        plt.show()
        plt.clf()
        klctr+=1
    return


def clean_cuts(wl, keeplist, subdir='./dq_cuts/', pctcuts=[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]):
    """
    Deletes superfluous intermediate files for image quality cuts that are not going to be used for analysis
    wl: 'Line' or 'Cont' or other file prefix
    subdir: path to directory where files live
    keeplist: list of image quality cuts to keep
    pctcuts: full list of image quality cuts
    """

    for pct in pctcuts:
        if pct not in keeplist:
            list1=glob.glob(subdir+'*'+str(pct)+'pctcut*.fits')
            for file in list1:
                os.remove(file)
            list2 = glob.glob(subdir + 'contrastcurves/*' + str(pct) + 'pctcut*')
            for file in list2:
                os.remove(file)
            shutil.rmtree(subdir+wl+'_'+str(pct)+'pctcut_sliced')

    return

def clean_fakes(keepstr, fakesdir):
    """
    """
    full_list = glob.glob(fakesdir+'/*')
    keep_list = [keep for keep in full_list if keepstr in keep]
    rm_list = [junk for junk in full_list if junk not in keep_list]
    print(full_list, keep_list, rm_list)
    return



def inject_fakes(data_str, cut, IWA, wl='Line', outputdir='fakes/', numann=6, movm=1, KLlist=[1,2,3,4,5,10,20,50,100],
                 contrasts=[1e-2,1e-2,1e-2], seps=[10, 10, 10], thetas=[0, 120, 240], debug=False,
                 ghost=False, mask=[3, 15], slicefakes=True,ctrlrad=30, highpass=True, meansnr=False, timecoll='median', smooth=False):
    """
    PURPOSE
    Injects false planets at specified locations, runs through KLIP, calculates throughput and
    displays a test image. This is largely the same as compute_thrpt, but for specific injected planet locations

    REQUIRED INPUTS
    data_str = unique identifier for this dataset for naming and table entry
    cut = what percentage of images do you want to cyt?
    fwhm = fwhm of the dataset (we use fwhm of median determined with idl atv aperture photometry)

    OPTIONAL INPUTS
    wl = wavelength you are reducing. set to "Line" by default, but can also specify "Cont"
    IWA = interior mask for KLIP
    seps = separations (in pixels) to inject planets
    thetas = PAs to inject planets
    numann, movm, KLlist = the usual klip annuli, movement and KL mode parameters
    debug = if set to True, prints some sanity checks about injected planet locations, etc.
    ghost = when set to True, will inject scaled copy of ghost as fake planet
    meansnr = return average SNR of positive pixels under planet mask instead of peak pixel

    OUTPUT
    separations where planets were injected = separations in pixels where throughputs were computed
    thrpt_list = computed throughputs at those separations (injected/recovered fake planet brightness)

    written by Kate Follette June 2019
    """
    #figure out what the fname is

    #if doesn't exist yet, make it
    if not os.path.exists('dq_cuts/' + wl + '_' + str(cut) + 'pctcut_sliced'):
        print("this wavelength and cut has not yet been generated. making now.")
        peak_cut(data_str, wl, pctcuts=[cut], ghost=ghost, rerun=True)
    
    # if contrast curve directory doesn't already exist, create it
    if os.path.exists(outputdir) == False:
        os.mkdir(outputdir)
    prefix=make_prefix(data_str, wl, cut)
    klstr=str(KLlist)
    klstr=klstr[1:-1]
    klstr=klstr.replace(', ','_')
    strklip = '_a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA) + '_KL' + klstr
    
    sepstr=''
    thetastr = ''
    ctrststr = ''
    for sep in seps:
        sepstr+=str(sep)+'_'
    for theta in thetas:
        thetastr+=str(theta)+'_'
    for contrast in contrasts:
        ctrststr+=str(contrast)+'_'
    prefix_fakes = prefix +'_thetas_'+thetastr+'seps_'+sepstr+'ctrst_'+ctrststr+'FAKES'

    ###load data
    filelist = glob.glob('dq_cuts/' + wl + '_' + str(cut) + 'pctcut_sliced' + "/sliced*.fits")
    # sorts file list so it's sequential
    filelist.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    ### pull the values of the star peak from the headers
    starpeak = []
    for i in np.arange(len(filelist)):
        head = fits.getheader(filelist[i])
        starpeak.append(head["STARPEAK"])
    head = fits.getheader(filelist[0])
    fwhm=head["MEDFWHM"]

    # import the dataset. Wait to highpass filter until the KLIP call in the next cell
    dataset = MagAO.MagAOData(filelist, highpass=False)
    # set IWA
    dataset.IWA = IWA

    # divide by starpeak values again to get input images in contrast units
    for i in np.arange(0, len(starpeak)):
        dataset.input[i, :, :] /= starpeak[i]

    # compute number of planets to inject
    n_planets = len(seps)
    if debug == True:
        print('I will inject ', n_planets, ' planets at separations', seps)

    size = dataset.input.shape
    nims = size[0]
    xcen = int((size[1] - 1) / 2)
    ycen = int((size[2] - 1) / 2)

    # replace any nans in image with zeros (edges, padded) for fake injection
    dataset_copy = np.copy(dataset.input)
    dataset.input[np.isnan(dataset.input) == True] = 0.
    imsz = dataset.input.shape[1]

    # now inject planets
    i = 0
    stmpsz = imsz
    starstamp = dataset.input[:, int(ycen - stmpsz / 2):int(ycen + stmpsz / 2), int(xcen - stmpsz / 2):int(xcen + stmpsz / 2)]

    for sep in seps:
        if ghost == True:
            if i==0:
                print("injecting planets into saturated data with fwhm=", fwhm)
            fakes.inject_planet(dataset.input, dataset.centers, np.repeat(contrasts[i], nims)[0], dataset.wcs, sep,
                                thetas[i], fwhm=fwhm)
        else:
                    # where can, use actual image as "stamp" to avoid edge effects
            # multiply starstamp image by contrast to create your false planets
            to_inject = contrasts[i] * starstamp
            fakes.inject_planet(dataset.input, dataset.centers, to_inject, dataset.wcs, sep, thetas[i], fwhm=fwhm,
                                stampsize=imsz)  # , thetas=thetas)
        i += 1

    # put NaNs back
    dataset.input[np.isnan(dataset_copy) == True] = np.nan

    #write out fakes cube
    fits.writeto(outputdir + prefix_fakes + '.fits', dataset.input, overwrite=True)

    #slice the final fake dataset in preparation for parameter exploration
    if slicefakes==True:
        outdir=outputdir+prefix_fakes+'_sliced/'
        if os.path.exists(outdir)== False:
            os.mkdir(outdir)
        for i in np.arange(0,len(dataset.prihdrs)):
            fits.writeto(outdir+'sliced_'+str(i)+'.fits', dataset.input[i,:,:], dataset.prihdrs[i], overwrite=True)

    # KLIP dataset with fake planets. Highpass filter here.
    parallelized.klip_dataset(dataset, outputdir=outputdir, fileprefix=prefix_fakes+strklip, algo='klip', annuli=numann,
                              subsections=1, movement=movm, numbasis=KLlist, calibrate_flux=False, mode="ADI",
                              highpass=highpass, save_aligned=False, time_collapse=timecoll, verbose = False, maxnumbasis=100)

    # read in the KLIP cube that was just created
    klcube = fits.getdata("{out}/{pre}-KLmodes-all.fits".format(out=outputdir, pre=prefix_fakes+strklip))
    # compile specs for the injected planets
    planetspecs = (seps, thetas, mask)
    # make a SNRMap cube
    #print(fwhm, outputdir, prefix_fakes,strklip, planetspecs, klcube.shape)
    Output, snrs, snr_sums, snr_spurious = snr.create_map(klcube, fwhm, saveOutput=True, outputName=outputdir+prefix_fakes+strklip + '_SNRMap.fits',
                                       planets=planetspecs, checkmask=False, method='all',ctrlrad=ctrlrad, smooth=smooth)
    # create a list that will store throughputs
    thrpt_list = []

    # loop through planet locations and retrieve flux
    i = 0
    nimages = dataset.output.shape[1]
    fake_fluxes = np.zeros((n_planets, nimages))

    for sep in seps:
        theta = thetas[i]
        contrast = contrasts[i]
        # thetas for retrieve planet call are measured from +x axis, so should always add 90 to pa for thetas keyword
        # dataset shape is [KL modes, n images, 1(wl dim), xdim, ydim]
        guesspeak = contrast * np.median(dataset.star_flux)
        
        fake_flux = fakes.retrieve_planet_flux(dataset.output[0, :, 0, :, :], dataset.centers, dataset.wcs, sep, theta,
                                               searchrad=8, guessfwhm=fwhm,
                                               guesspeak=guesspeak,
                                               thetas=np.repeat(theta + 90, dataset.output.shape[1]), refinefit=True)
        print('replacing', len(fake_flux[fake_flux < 0]), 'negative fluxes with NaN')
        fake_flux[fake_flux < 0]=np.nan

        fake_fluxes[i, :] = fake_flux
        newthpt = np.nanmedian(fake_flux / (contrast))
        
        if (newthpt > 1) or (newthpt < 0):
            print('invalid throughput value of', newthpt, 'for planet at separation', sep, '. Replacing with NaN')
            newthpt = np.nan
        thrpt_list.append(newthpt)
        if debug == True:
            print("fake planet at ", sep, " pixels has a mean throughput of", np.nanmean(fake_flux / (contrast)),
                  "median of", np.nanmedian(fake_flux / (contrast)), " and stdev of ",
                  np.nanstd(fake_flux / (contrast)))
        i += 1
    print("throughputs are", thrpt_list)

    #print("snrs (max pixel) for median SNR map are", snrs[0,:,:])
    print("stdev SNR map results")
    for pl in np.arange(n_planets):
        print("average peak SNR across KL modes for planet", str(pl+1), 'is', np.nanmean(snrs[0,:,pl]), 'and median is', np.nanmedian(snrs[0,:,pl]))

    #print("snrs (avg under mask) for median SNR map is", snr_sums)
    if debug==True:
        for pl in np.arange(n_planets):
            print("average average under mask SNR across KL modes for planet", str(pl+1), 'is', np.nanmean(snr_sums[0,:,pl]), 'and median is', np.nanmedian(snr_sums[0,:,pl]))

        print("median SNR map results")
        #print("snrs (max pixel) for stdev SNR map are", snrs[0,:,:])
        for pl in np.arange(n_planets):
            print("average peak SNR across KL modes for planet", str(pl+1), 'is', np.nanmean(snrs[1,:,pl]), 'and median is', np.nanmedian(snrs[1,:,pl]))

        #print("snrs (avg under mask) for stdev SNR map is", snr_sums)
        for pl in np.arange(n_planets):
            print("average average under mask SNR across KL modes for planet", str(pl+1), 'is', np.nanmean(snr_sums[1,:,pl]), 'and median is', np.nanmedian(snr_sums[1,:,pl]))

    if meansnr==True:
        snrs = snr_sums
    return (dataset.input, dataset.prihdrs, prefix_fakes, snrs)


def radarr(xdim,ydim):
    '''
    make an image of size xdim x ydim where every pixel has a value equal to its distance from the center of the array
    '''
    im = np.zeros((xdim, ydim))
    xcen = (xdim - 1) / 2.
    ycen = (ydim - 1) / 2.
    for index1 in range(xdim):
        for index2 in range(ydim):
            im[index1, index2]=np.sqrt((xcen-index1)**2. + (ycen-index2)**2.)
    return(im)

def ctrlmask(xdim, ydim, rin, rout):
    '''
    makes a xdim x ydim mask where pixels with distances from center between rin and
    rout have been replaced by NaNs
    '''
    arr = radarr(xdim, ydim)
    new = np.ones((xdim, ydim))
    new[(arr >= rin) & (arr <= rout)] = np.nan
    return(new)

def paramexplore_fig(pedir, pename, kllist, writestr=False, weights=[1,1,1,1,1,1], snrmeth='all', smt=3, separate_planets=False, separate_kls=False):
    
    """
    Right now only works with either separate_kls=False or a single KL mode
    Likewise, will only work with planets collapsed or will show innermost (index 0) planet
    """

    metric_cube, agg_cube, ann_val, movm_val, metric_scores, metric_fname = \
        pe.find_best_new(pename, kllist, pedir=pedir, writestr=writestr, weights=weights, snrmeth=snrmeth, smt=smt, separate_planets=separate_planets, separate_kls=separate_kls)

    #reduce dimensions of metric cube (see note above)
    metric_cube = metric_cube[:,0,0,0,:,:]

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

    # set up tick labels according to parameter ranges
    fig_xdim = nstepx*0.35
    fig_ydim = nstepy

    fig = plt.figure(tight_layout=True, figsize=(fig_ydim,fig_xdim))
    gs = fig.add_gridspec(2, 5, hspace=0.2, wspace=0.5)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    ax5 = fig.add_subplot(gs[0,2])
    ax6 = fig.add_subplot(gs[1,2])
    ax7 = fig.add_subplot(gs[:,3:])
    #ax5 = fig.add_subplot(gs[:,2:])

    plt.setp((ax1, ax2, ax3, ax4, ax5, ax6, ax7), xticks=np.arange(0,nstepx + 1,2), xticklabels=np.arange(xmin, xmax + 1,2),
             yticks=np.arange(0, nstepy + 1,2), yticklabels=np.arange(ymin, ymax + 1,2))

    labelstyle=dict(size=16)
    titlestyle=dict(size=18)
    bigtitlestyle=dict(size=24)

    #plt.setp((ax1, ax2, ax3, ax4, ax5), xticks=np.arange(nstepx + 1), yticks=np.arange(nstepy + 1), xticklabels=[], yticklabels=[])
    im1 = ax1.imshow(metric_cube[0,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)

    ax1.set_xlabel("movement parameter", **labelstyle)
    ax1.set_ylabel("annuli parameter", **labelstyle)
    ax1.set_title("Peak SNR", **titlestyle)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = plt.colorbar(im1, cax=cax, orientation='vertical')
    cb.set_label(label='Normalized Score', size=16)
    cb.ax.tick_params(labelsize='large')

    im2 = ax2.imshow(metric_cube[1,:,:], origin='lower',
                     cmap='magma', vmin=0, vmax=1)
    ax2.set_xlabel("movement parameter", **labelstyle)
    ax2.set_ylabel("annuli parameter", **labelstyle)
    ax2.set_title("Peak SNR Neighbor Quality", **titlestyle)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = plt.colorbar(im2, cax=cax, orientation='vertical')
    cb.set_label(label='Normalized Score', size=16)
    cb.ax.tick_params(labelsize='large')

    im3 = ax3.imshow(metric_cube[2,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
    ax3.set_xlabel("movement parameter", **labelstyle)
    ax3.set_ylabel("annuli parameter", **labelstyle)
    ax3.set_title("Avg. SNR Under Planet Mask", **titlestyle)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = plt.colorbar(im3, cax=cax, orientation='vertical')
    cb.set_label(label='Normalized Score', size=16)
    cb.ax.tick_params(labelsize='large')
    
    im4 = ax4.imshow(metric_cube[3,:,:], origin='lower', cmap='magma', vmin=0, vmax=1)
    ax4.set_xlabel("movement parameter", **labelstyle)
    ax4.set_ylabel("annuli parameter", **labelstyle)
    ax4.set_title("Avg. SNR Neighbor Quality", **titlestyle)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = plt.colorbar(im4, cax=cax, orientation='vertical')
    cb.set_label(label='Normalized Score', size=16)
    cb.ax.tick_params(labelsize='large')

    
    im5 = ax5.imshow(metric_cube[4,:,:], origin='lower',
                     cmap='magma', vmin=0, vmax=1)
    ax5.set_xlabel("movement parameter", **labelstyle)
    ax5.set_ylabel("annuli parameter", **labelstyle)
    ax5.set_title(r"5$\sigma$ Pixels IWA to CR", **titlestyle)
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = plt.colorbar(im5, cax=cax, orientation='vertical')
    cb.set_label(label='Normalized Score', size=16)
    cb.ax.tick_params(labelsize='large')

    im6 = ax6.imshow(metric_cube[5,:,:], origin='lower',
                     cmap='magma', vmin=0, vmax=1)
    ax6.set_xlabel("movement parameter")
    ax6.set_ylabel("annuli parameter")
    ax6.set_title("Contrast", **titlestyle)
    divider = make_axes_locatable(ax6)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = plt.colorbar(im6, cax=cax, orientation='vertical')
    cb.set_label(label='Normalized Score', size=16)
    cb.ax.tick_params(labelsize='large')


    # plot metric
    agg=metric_cube[6,:,:]
    im7 = ax7.imshow(agg, origin='lower', vmin=1, vmax=np.nanmax(agg))
    ax7.set_ylabel("annuli parameter", **titlestyle)
    ax7.set_xlabel("movement parameter", **titlestyle)
    ax7.set_title("Aggregate Parameter Quality Metric", **bigtitlestyle)
    divider = make_axes_locatable(ax7)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im7, cax=cax, orientation='vertical', label="")
    cb = plt.colorbar(im7, cax=cax, orientation='vertical')
    cb.set_label(label='Weighted Aggregate Score', size=20)
    cb.ax.tick_params(labelsize=16)

    ind = np.where(agg == np.nanmax(agg))
    label_text = 'a' + str(int(ann_val[0][0])) + 'm' + str(int(movm_val[0][0]))
    rect = patches.Rectangle((ind[1][0] - 0.5, ind[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='r', facecolor='none')
    ax7.add_patch(rect)
    ax7.text(ind[1][0] + 0.75, ind[0][0]-0.5, label_text, color='red', **titlestyle, weight='bold')

    plt.savefig(pedir+writestr+'_paramqual.png')
    
    return(ann_val, movm_val, agg)


def get_pe_df(dfname):
    # define dataframe if doesn't already exist
    try:
        df=pd.read_csv(dfname)
        print("found existing data frame") # with the following contents: \n", df)
        df_cols = df.columns.tolist()
    except:
        print("creating new data frame")
        #if none found, make one
        df_cols =['Dataset','Object', 'Date','subset','Wavelength','pctcut', 'optimal movement', 'optimal annuli', 'IWA', 'KL combo used',
        'SNR value', 'SNR metric', 'SNR neighbors metric', 'stdev metric', 'stdev neighbors metric', 'spurious pixels metric', 'total metric', 'date added']
        # make a blank pandas dataframe
        df = pd.DataFrame({})
        for name in df_cols:
            df[name] = []
    return(df, df_cols)

def add_to_pe_df(data_str, pedir, pename, kllist, weights=[1,1,1,1,1,1], pe_dfname='../../optimal_params.csv'):
    df, df_cols =get_pe_df(pe_dfname)
    avgkl, stdevkl = collapsekl(pedir, pename, kllist)
    snr_norm, nq_snr, stdev_norm, nq_stdev, spurpix, agg, ann_val, movm_val, metric_scores, avgSNR = pe.find_best_new(pedir, pename, kllist, weights=weights)
    ind=np.where(agg == np.nanmax(agg))
    data_split=re.split('_',data_str)
    objname = data_split[0]
    date = data_split[1]
    try:
        subset = data_split[2]
    except:
        subset = ''

    if "Line" in pename:
        wl = "Line"
    elif "Cont" in pename:
        wl = "Cont"
    else:
        "Could not parse wavelength from filename. please enter"
        wl = input("please enter wavelength:")

    try:
        head = fits.getheader(pedir+pename)
        print("found header")
        IWA = head["IWA"]
    except:
        IWA = int(input("You are using an old Parameter Explorer file. Please enter IWA:"))

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    pedir_info = re.split('/|_', pedir)

    cut=[s for s in pedir_info if 'cut' in s]
    cut=(cut[0].split('pctcut'))[0]
    datalist=pd.DataFrame([[data_str,objname,date,subset,wl,cut,movm_val[0],ann_val[0],IWA,str(kllist),avgSNR[0],snr_norm[ind][0], nq_snr[ind][0],stdev_norm[ind][0], nq_stdev[ind][0], spurpix[ind][0], agg[ind][0], dt_string]], columns=df_cols)
    df = df.append(datalist)
    df.to_csv(pe_dfname, index=False)
    return(df)

def clean_pe(pedir, keepparams):
    """
    Deletes superfluous parameter explorer that are not going to be used for analysis
    prefix: file prefix
    pedir: path to directory where parameter explorer files live
    keeplist: list of files to keep in [annuli, movem] form
    """
    fulllist = glob.glob(pedir+'*.fits')
    keeplist=[]
    for keep in keepparams:
        lst = glob.glob(pedir+'*a'+str(keep[0])+'m'+str(keep[1])+'.*.fits')
        keeplist+=lst
    pefiles = glob.glob(pedir+'paramexplore*.fits')
    keeplist+=pefiles
    for file in fulllist:
        if file not in keeplist:
            os.remove(file)

    return


def get_klip_inputs(data_str_uniq, pe_dfname='../../optimal_params.csv', rdx_params_dfname='rdx_params.csv'):
    """
    pulls info for naming and KLIP parameters from two data frames storing this info
    pe_dfname: dataframe containing info about optimal parameters
    data_str_uniq = unique dataset identifier 
    :return:
    """
    df, df_cols =get_pe_df(pe_dfname)
    df2=get_df(rdx_params_dfname)
    name_split=re.split('_',data_str_uniq)
    objname = name_split[0]
    date = name_split[1]

    match = [line for line in df["Dataset"] if data_str_uniq == line]

    if len(match) > 1:
        print('there is more than one entry in the dataframe for this dataset. please delete the duplicates.')
        return
    else:
        cut = int(df[df["Dataset"].values == match]["pctcut"].values[0])
        movm = int(df[df["Dataset"].values == match]["optimal movement"].values[0])
        numann = int(df[df["Dataset"].values == match]["optimal annuli"].values[0])
        fwhm = df2[df2["Dataset"] == data_str_uniq]["medfwhm"].values[0]
        IWA = int(df[df["Dataset"] == data_str_uniq]["IWA"].values[0])
        kllist = df[df["Dataset"] == data_str_uniq]["KL combo used"].values[0]
        kllist = eval(kllist)
    return (objname, date, cut, movm, numann, fwhm, IWA, kllist)


def klip_data(data_str, wl, params=False, fakes=False, planets=False, highpass=True, overwrite=False, klinput = False, indir='dq_cuts/', outputdir='final_ims/', ctrlrad=30, timecoll='median', smooth=False):

    if os.path.exists(outputdir) == False:
        os.mkdir(outputdir)

    if params == False:
        objname, date, cut, movm, numann, fwhm, IWA, kllist = get_klip_inputs(data_str)
    else:
        [objname, date, cut, movm, numann, fwhm, IWA, kllist] = params

    if klinput != False:
        kllist = klinput
    
    if fakes==False:
        namestr = wl + '_' + str(cut) +'pctcut_sliced'
    else:
        namestr = fakes + '_sliced'
    
    #if hasn't already been sliced, slice it
    slicedir = indir + namestr + '/'

    if os.path.exists(slicedir) == False:
        preims = glob.glob('preprocessed/'+wl+'*nocosmics.fits')
        preims = preims[0].split(wl)
        imstring = preims[1][:-5]
        print(wl, " image has not yet been sliced. Slicing now.")
        imname = indir+wl+imstring+'_'+str(cut)+'pctcut.fits'
        rotoff_name = indir+'rotoff_no'+wl+'cosmics_'+str(cut)+'pctcut.fits'
        SliceCube(imname, rotoff_name, slicedir=slicedir)
        pykh.addstarpeak(slicedir, debug=True, mask=True)

    prefix =namestr[:-6] + 'a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA) + 'hp' + str(highpass)
    if kllist==[1,2,3,4,5,10,20,50,100]:
      klstr='all'
    else:
      if isinstance(kllist,int):
        klstr='_'+str(kllist)
      else:
        klstrlist = [str(kl) for kl in kllist]
        klstr='_'.join(klstrlist)
    prefix+='_kl'+klstr

    #check whether this KLIP image has already been generated
    runrdx=True
    if (os.path.exists(outputdir+prefix+'-KLmodes-all.fits')) and (overwrite==False):
        print('a file with the name', outputdir+prefix+'-KLmodes-all.fits', 'already exists')
        klcube_header = fits.getheader(outputdir+prefix+'-KLmodes-all.fits')
        #pyklip writes out each KL mode value to separate keyword
        kls = klcube_header["KLMODE*"]
        #make a list and populate it with KL mode values
        im_kls = []
        for i in np.arange(len(kls)):
            im_kls.append(klcube_header["KLMODE"+str(i)])

        #catch for single Kl mode cubes
        if not isinstance(kllist, list):
          kllist=[kllist]

        if im_kls == kllist:
            print('and KL modes match')
            klcube=fits.getdata(outputdir+prefix+'-KLmodes-all.fits')
            runrdx=False
        else:
            print('but KL modes don\'t match. Regenerating.')
            runrdx = True

    if runrdx==True:   
        filelist = glob.glob(slicedir + "/sliced*.fits")
        dataset = MagAO.MagAOData(filelist, highpass=False) 
        parallelized.klip_dataset(dataset, outputdir=outputdir, fileprefix=prefix, algo='klip', annuli=numann, subsections=1, movement=movm,
                              numbasis=kllist, calibrate_flux=False, mode="ADI", highpass=highpass, save_aligned=False, time_collapse=timecoll,
                              maxnumbasis=100)
        klcube = fits.getdata("{out}{pre}-KLmodes-all.fits".format(out=outputdir, pre=prefix))
    
    #check whether this SNRMap has already been generated
    #if os.path.exists(outputdir+prefix + '_SNRMap.fits'):
        #snmap = fits.getdata(outputdir+prefix + '_SNRMap.fits')
    #else:
    if planets != False:
        snmap, snrs, snr_sums, snr_spurious = snr.create_map("{out}{pre}-KLmodes-all.fits".format(out=outputdir, pre=prefix), fwhm, planets=planets, saveOutput=True, outputName=outputdir+prefix + '_SNRMap.fits', ctrlrad=ctrlrad, method='stdev', smooth=smooth)
    else: 
        snmap = snr.create_map("{out}{pre}-KLmodes-all.fits".format(out=outputdir, pre=prefix), fwhm, planets=planets, saveOutput=True, outputName=outputdir+prefix + '_SNRMap.fits', method="stdev", smooth=smooth)
    #snmap, snrs, snr_sums, snr_spurious = snr.create_map(klcube, fwhm, saveOutput=True, outputName=prefix + '_SNRMap.fits')
    return (klcube, snmap, fwhm)


def get_scale_factor(data_str, scalefile = '../../GAPlanetS_Dataset_Table.csv'):
    #scalefile=scalefile.decode('utf-8')
    df3 = pd.read_csv(scalefile)
    data_path = data_str.replace('_','/',1)
    scale = df3[df3["Folder Name"] == data_path]["Scale Factor"].values[0]
    return (scale)


def run_redx(data_str, scale = False, indir='dq_cuts/', highpass=True, params=False, outputdir = 'final_ims/', klinput=False, scalefile = '../../GAPlanetS_Dataset_Table.csv', overwrite=False, timecoll='median', smooth=False, planets=False):
    wls = ['Line', 'Cont']
    if params == False:
        objname, date, cut, movm, numann, fwhm, IWA, kllist = get_klip_inputs(data_str)
    else:
        objname, date, cut, movm, numann, fwhm, IWA, kllist = params

    #print(data_str, imstring, indir, outputdir)
    linecube, linesnr, linefwhm = klip_data(data_str, wls[0], indir=indir, outputdir=outputdir, klinput=klinput, params=params, highpass=highpass, overwrite=overwrite, timecoll=timecoll, smooth=smooth, planets=planets)
    contcube, contsnr, contfwhm = klip_data(data_str, wls[1], indir=indir, outputdir=outputdir, klinput=klinput, params=params, highpass=highpass, overwrite=overwrite, timecoll=timecoll, smooth=smooth, planets=planets)
    
    if scale == False:
        print('pulling scale from file')
        scale = get_scale_factor(data_str, scalefile=scalefile)
    else:
        print('setting scale to', scale)
        scale = float(scale)
    
    sdicube = linecube - scale * contcube
    prefix = data_str + '_' + str(cut) + 'pctcut_' + 'a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA)+ 'hp'+str(highpass)

    if kllist==[1,2,3,4,5,10,20,50,100]:
      klstr='all'
    else:
      if isinstance(kllist,int):
        klstr='_'+str(kllist)
      else:
        klstrlist = [str(kl) for kl in kllist]
        klstr='_'.join(klstrlist)
    prefix+='_kl'+klstr

    sdisnr = snr.create_map(sdicube, (linefwhm + contfwhm) / 2., saveOutput=True, outputName=outputdir+prefix +'_SDI_scl'+'{:.2f}'.format(scale) + '_SNRMap.fits', method='stdev', smooth=smooth, planets=planets, checkmask=False)
    if planets!=False:
        sdisnr = sdisnr[0]

    return (linecube, linesnr, contcube, contsnr, sdicube, sdisnr, prefix, scale)

def grab_planet_specs(df,dset_path):
    """
    Given (properly formatted) dataframe, extract locations and designations for planet(s).
    """
    pllabel = df[df["Path"]==dset_path]["Designation"].values
    plsep = df[df["Path"]==dset_path]["Separation (pix)"].values
    seperr = df[df["Path"]==dset_path]["Sep Error"].values
    plpa = df[df["Path"]==dset_path]["PA"].values
    #for later. this is how to make a wedge-shaped patch. pull pa uncert as well
    #Wedge((.8, .3), .2, 45, 90, width=0.10),  # Ring sector
    return(tuple(pllabel),tuple(plsep),tuple(plpa), tuple(seperr))

def bulk_rdx(sorted_objs, wl, outdir, scalefile, df, base_fpath='/content/drive/Shareddrives/',hpmult=0.5,
    klopt=False, ctrstopt=False, weights=[1,1,1,1,1,1], kllist_coll = [10,100], overwrite=False, maxx=24, maxy=25, 
    minx=0, miny=1, skipdates=False, pldf=False, seppl=False, sepkl=False, smt=1, timecoll='median',combochoice=None,
    smooth=False, haopt=False, maskspec=None):
    
    thisdir=os.getcwd()

    ##turn keywords into unique string
    wts = weights
    wtstr = ''.join(str(val) for val in wts)
    klstr = '_kl'+'_'.join(str(val2) for val2 in kllist_coll)   
    hpstr = "_"+str(hpmult)+"fwhm"
    theseparams ='wts'+wtstr+klstr+'_seppl'+str(seppl)+'_sepkl'+str(sepkl)+hpstr
    if minx!=0:
        theseparams+='_minx'+str(minx)
    if miny!=1:
        theseparams+='_miny'+str(miny)
    if maxx!=24:
        theseparams+='_maxx'+str(maxx)
    if maxy!=25:
        theseparams+='_maxy'+str(maxy)
    outname = outdir+'Optimal_Vals_Coll_'+theseparams+'.csv'
    
    n_dsets = []
    #count how many total datasets
    for obj in sorted_objs:
        filt = df[df["Object Name"]==obj][df["Done"]=='Y']
        if skipdates!=False:
            n_dsets.append(len(filt[~filt.Date.isin(skipdates)]))
        else:    
            n_dsets.append(len(filt))

    totaldsets = np.sum(n_dsets)

    dofds={}
    #dset counter
    i=0 #plot counter

    dupobj=[]
    lastind=[]
    dupyr=[]
    dup_dstrings=[]

    for obj in sorted_objs:
        if skipdates==False:
            subdf = df[df["Object Name"]==obj][df["Done"]=='Y']
        else:
            subdf = df[df["Object Name"]==obj][~df.Date.isin(skipdates)][df["Done"]=='Y']
        
        n_dset=(len(subdf))
        #list of paths for theis object
        dpaths=[]
        for dset in np.arange(n_dset):
            dpaths.append(subdf["Path"].values[dset])
        
        #put in date order
        monthstrs = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
        t=[]
        line_dates=[]
        for j in np.arange(n_dset):
            linefname = dpaths[j].split('/')
            line_date = linefname[1]
            line_dates.append(line_date)
            if '_' in line_date:
                line_date = line_date.split('_')[0]
            monthstr = [s for s in monthstrs if s in line_date]
            monthno = [m for m in np.arange(1,13) if monthstrs[m-1] in line_date]
            year = int('20'+line_date[-2:])
            day = line_date.split(monthstr[0])[0]
            #print(line_date, '= day:', day, ' month:', monthno[0], ' year:', year)
            t.append(dt.date(year=year, month=int(monthno[0]), day=int(day)))

        yrlist=[]
        gencombo=False

        for date in t:
            yrlist.append(date.year)

        unique, counts = np.unique(yrlist, return_counts=True)
        for yr in np.arange(len(unique)):
            if counts[yr]>1:
                print("more than one", obj, "entry for year", unique[yr])
                if combochoice!=True:
                    if combochoice==False:
                        gencombo=False
                    else:
                        answ = input('Should I combine (Y/N)?')
                        if answ=='Y' or 'y':
                            combochoice=True
                if combochoice==True:
                    gencombo = True
                    #find last occurrence
                    dupobj.append(obj)
                    li = len(yrlist)-1-yrlist[::-1].index(unique[yr])
                    lastind.append(li)
                    dupyr.append(yrlist[li])
                    dlist = [line_dates[i] for i in np.arange(len(yrlist)) if yrlist[i]==unique[yr]]
                    dup_dstrings.append(dlist)

        #now sort
        sorted_dpaths = [x for _, x in sorted(zip(t, dpaths))]

        #pull from correct Drive directory
        whichdir = df[df["Object Name"]==obj]['Which AWS'].values[0]
        if whichdir == 1:
            fpath_st = base_fpath+'Follette-Lab-AWS/GAPlanetS_Database/'
        elif whichdir ==2:
            fpath_st = base_fpath+'Follette-Lab-AWS-2/'

        #pull info for each dataset
        for j in np.arange(n_dset):
            #by dataset dictionary
            d = {}
            dset_path = sorted_dpaths[j]
            date=dset_path.split('/')[1]
            #check whether pes are complete
            go=True

            if df[df["Path"]==dset_path]["Done"].values[0]!='Y':
                print('datset not fully processed')
                go=False

            if go==True:
                #optimal values from df
                cut = int(df[df["Path"]==dset_path]["cut chosen"].values[0])
                dq_str = str(cut)+'pctcut'
                #pull various other inputs
                IWA = int(df[df["Path"]==dset_path]["IWA"].values[0])
                fwhm = int(df[df["Path"]==dset_path]["cut fwhm"].values[0])
                contrast = float(df[df["Path"]==dset_path]["injected planet contrast"].values[0])
                satrad = df[df["Path"]==dset_path]["saturation radius"].values[0]
                whichbin = df[df["Path"]==dset_path]["bin"].values[0]
 
                if np.isfinite(whichbin):
                    if int(whichbin)==1:
                        ctrlrad=30
                    elif int(whichbin)==2:
                        ctrlrad=15
                    else:
                        print('bin must be 1 or 2')
                        ctrlrad=None
                else:
                    "no bin specified - can't plot ctrlrad"
                    ctrlrad=None
                
                full_fpath = fpath_st+dset_path

                d["IWA"]=IWA
        
                #set ghost to true if saturated
                if np.isnan(satrad):
                    ghost=False
                    satrad=np.nan        
                else:
                    ghost=True
                    satrad=int(satrad)

                #path to the relevant line directory
                fpath=full_fpath+'/dq_cuts/'+wl+'_'+dq_str+'_sliced/'
                #dataset label = string of filepath
                data_str = dset_path.replace('/','_')

                print('STARTING', data_str)

                #find optimal params according to the collapse strategy
                #find the cont fakes klip dir with paramexplore file
                if haopt==True:
                    fakedir=full_fpath+'/dq_cuts/'+wl+'_'+dq_str+'_sliced_klip/'
                    print(fakedir)
                else:
                    fakepath=data_str+'fakes_redo/'
                    fakedirfiles = glob.glob(full_fpath+'/'+fakepath+'*')
                    dirnames = [d for d in fakedirfiles if os.path.isdir(d)]
                    fakedir = [dir for dir in dirnames if dir[-4:] == 'klip']
                    if len(fakedir) > 1:
                        print('more than one klipped fakes dir.')
                        for opt in np.arange(len(fakedir)):
                            print('option', opt, ': ', fakedir[opt])
                        choiceind = input('Which do you want?')
                        fakedir=fakedir[int(choiceind)]
                    else:
                        fakedir=fakedir[0]
                    
                pefile = glob.glob(fakedir+'/param*.fits')
                if len(pefile) > 1:
                    print('more than one pefile')
                else:
                    pefile=pefile[0]
                    pefname=pefile.split('_klip/')
                    pedir=pefname[0]+'_klip/'
                    pefname = pefname[1]

                #set hpval based on this dset's fwhm
                hpval = hpmult*fwhm    
                 
                #find peak for this set according to kl collapse mode in order to select optimal ann, movm
                cube, agg_cube, anns, movms, scores, mfname = pe.find_best_new(pefname,kllist_coll,pedir=pedir,weights=wts,snrmeth='stdev', outdir=outdir+data_str+'/', separate_planets=seppl, separate_kls=sepkl, maxx=maxx, maxy=maxy, minx=minx, miny=miny)

                ann = int(anns[0][0])
                movm = int(movms[0][0])
                print('peak is at a', ann, 'm', movm)

                #run metric calculation for ALL kl modes so can select optimal kl
                kllist = [1,2,3,4,5,10,20,50,100]     
                cube1, agg_cube1, anns1, movms1, scores1, mfname1 = pe.find_best_new(pefname,kllist,pedir=pedir,weights=wts,snrmeth='stdev', outdir=outdir+data_str+'/', separate_planets=seppl,separate_kls=True, maxx=maxx, maxy=maxy, minx=minx, miny=miny)

                #find kl mode with highest score at ann, movm
                if klopt==True:
                    pk=np.zeros((len(kllist)))
                    for k in np.arange(len(kllist)):
                        #subtract 1 from ann to index from 0. note assumes step size is 1
                        pk[k] = np.nanmax(agg_cube1[k,0,ann-1,movm])
                
                    pk_ind = np.where(pk == np.nanmax(pk))
                    kl = kllist[int(pk_ind[0])]
                
                    print('peak metric for this ann, movm combo among', pk, 'is', pk[int(pk_ind[0])], 'for kl', kl)
                    df.loc[df["Path"]==dset_path,"opt kl"]=kl
                elif ctrstopt==True:
                    peuse=fits.getdata(pedir+pefname)
                    pk=np.zeros((len(kllist)))
                    for k in np.arange(len(kllist)):
                        #subtract 1 from ann to index from 0. note assumes step size is 1
                        #last index in pe is contrast
                        pk[k] = np.nanmin(np.log10(-1*peuse[-1,0,k,0,ann-1,movm]))
                
                    pk_ind = np.where(pk == np.nanmin(pk))
                    kl = kllist[int(pk_ind[0])]
                
                    print('best log contrast for this ann, movm combo among', pk, 'is', pk[int(pk_ind[0])], 'for kl', kl)
                    df.loc[df["Path"]==dset_path,"opt kl"]=kl

                else:
                    #otherwise, do 10 kl modes
                    kl = 10

                #fill in in the table
                df.loc[df["Path"]==dset_path,"opt ann"]=ann
                df.loc[df["Path"]==dset_path,"opt movm"]=movm

                if os.path.exists(outname):
                    dfin = pd.read_csv(outname)
                else:
                    dfin = pd.DataFrame(columns=df.columns)
                
                idx = df.index[df["Path"]==dset_path].tolist()
                dfin = dfin.append(df.loc[idx], ignore_index=True)
                dfin.to_csv(outname, index=False)     

                #compute contrast curve
                os.chdir(full_fpath)
                thrpt_out, zb, df2, dataset_prefix = compute_thrpt(data_str, wl, cut, outputdir = outdir+data_str+'/', numann=ann, movm=movm, KLlist=[kl], IWA=IWA, 
                                                                    contrast=contrast, theta=0., clockang=135, debug=False, record_seps=[0.1, 0.25, 0.5, 0.75, 1.0],
                                                                    ghost=ghost, savefig=True, iterations=3,rdx_params_dfname=outdir+'rdx_params'+hpstr+'.csv', 
                                                                    highpass=hpval, overwrite=overwrite, timecoll=timecoll, usecutfwhm=True)
            
                contrast_out, df2, OWA = make_contrast_curve(data_str, wl, cut, thrpt_out, dataset_prefix,  outputdir=outdir+data_str+'/', 
                                                                numann=ann, movm=movm, KLlist=[kl], IWA=IWA, rdx_params_dfname=outdir+'rdx_params'+hpstr+'.csv', 
                                                                record_seps=[0.1, 0.25, 0.5, 0.75, 1.0], savefig=True, debug=False, highpass=hpval, 
                                                                overwrite=overwrite, timecoll=timecoll)
                
            
                klipstr=' a'+str(ann)+'m'+str(movm)+'kl'+str(kl)
                date_obj = t[j]
                ccurve_lbl = obj + ' ' + date_obj.strftime("%m/%d/%Y") + klipstr
                ccurve_seps = contrast_out[0,0,:]
                ccurve_ctrst = contrast_out[0,2,:]

                #make a dictionary with all these things
                params = [obj,date,cut,movm,ann,fwhm,IWA,kl]
    

                #set planet marking keywords to False by default
                plspecs = False
                plcand = False

                ##compile a bunch of stuff for plotting, labeling, etc.

                planets=False
                #pull planet info, if given a planet specs frame
                if isinstance(pldf, pd.DataFrame):
                    #check whether df contains note to mark planets
                    if df[df["Path"]==dset_path]["mark known planets"].values[0]=="Y":
                        plspecs = grab_planet_specs(pldf,dset_path)
                        if maskspec==None:
                            maskspec=(fwhm/2,15)
                        planets=(plspecs[1],plspecs[2],maskspec)
                        #check whether marked as candidate
                        if df[df["Path"]==dset_path]["candidate"].values[0]=="Y":
                            plcand=True

                ##generate images
                linecube, linesnr, contcube, contsnr, sdicube, sdisnr, prefix, scale = run_redx(data_str,params=params, scalefile=scalefile, highpass=hpval, scale=False, overwrite=overwrite, timecoll=timecoll, smooth=smooth, planets=planets)
                linecube, linesnr, contcube, contsnr, sdicube2, sdisnr2, prefix, scale2 = run_redx(data_str,params=params, highpass=hpval, scale=1, overwrite=overwrite, timecoll=timecoll, smooth=smooth, planets=planets)


                dict_keys = ['prefix', 'full_fpath', 'dset_path', 'satrad', 'ctrlrad', 'scale', 
                'ccurve_seps', 'ccurve_ctrst', 'ccurve_lbl',
                'linecube', 'contcube', 'sdicube', 'sdicube2', 'linesnr', 'contsnr', 'sdisnr', 'sdisnr2', 
                'obj', 'date', 'date_obj', 'plspecs', 'plcand',
                'cut', 'movm', 'ann', 'fwhm', 'IWA', 'kl', 'hpval']

                for k in dict_keys:
                    d[k]=locals()[k]

                os.chdir(thisdir)
                #add dict for this dataset to dict of dicts
                dofds[data_str]=d
                i+=1

    if gencombo==True:
        dofds = add_dsetcombos(dofds,dupobj,dup_dstrings, dupyr)
        totaldsets+=(len(dupobj))

    objlist = '_'.join(sorted_objs)
    master_keys = ["theseparams","totaldsets","objlist","outdir"]
    for ky in master_keys:
        dofds[ky]=locals()[ky]

    dictname = outdir+objlist+'_'+theseparams+'.p'
    pickle.dump(dofds,open(dictname,"wb"))

    return(dictname, dofds)

def add_dsetcombos(dofds,dupobj,dup_dstrings, dupyr):
    
    #add a new combo dictionary for each combined dataset
    for ds in np.arange(len(dupobj)):
        combod = {}
        dicts = []

        #make a list of the dictionaries corresponding to this combo
        for i in np.arange(len(dup_dstrings[ds])):
            print(dupobj[ds],dup_dstrings[ds])
            thiskey =  dupobj[ds]+'_'+dup_dstrings[ds][i]
            print("key", thiskey)
            dicts.append(dofds[thiskey])
    
        #average line/cont snrmaps and ims
        keylist = ['linecube', 'contcube', 'sdicube', 'sdicube2', 'linesnr', 'contsnr', 'sdisnr', 'sdisnr2']
        for key in keylist:
            if 'snr' in key:
                foravg = np.zeros((len(dicts),1,1,451,451))
            else:
                foravg = np.zeros((len(dicts),1,451,451))
            for d in np.arange(len(dicts)):
                foravg[d,:,:]=dicts[d][key]
            combod[key]=np.nanmean(foravg,axis=0)

        #average scales
        quants = ["scale","IWA", "fwhm"]
        for q in quants:
            parlist = []
            for d in np.arange(len(dicts)):
                parlist.append(dicts[d][q])
            print(parlist)
            if q=="IWA":
                combod[q]=np.max(parlist)
            else:
                combod[q]=np.mean(parlist)

        #grab planet specs from first dict
        combod["plspecs"]=dicts[0]["plspecs"]
        combod["plcand"]=dicts[0]["plcand"]
        combod["ctrlrad"]=dicts[0]["ctrlrad"]
        combod["obj"]=dupobj[ds]

        newkey = dupobj[ds]+'_'+str(dupyr[ds])+'_combo'
        print("new key", newkey)
        combod["prefix"]=newkey

        dofds[newkey]=combod
        print(combod.keys())
    return(dofds)  

def plotdict_ims(d, snr=False, smt=False, lims=[-1,4], stampsz=75):
    """

    """
    #if d is string, open it
    if isinstance(d, str):
        d = pickle.load( open(d, "rb" ) )

    totaldsets = d["totaldsets"]
    f1 = plt.figure(figsize=(16.2,totaldsets*4.1))
    plt.axis('off')
    outer = gridspec.GridSpec(totaldsets, 1, wspace=0, hspace=0.01)

    dkeys = list(d.keys())
    dset_strs = []
    datestrs = []
    for key in dkeys:
        if '_' in key:
            key_split =key.split('_')
            #pull out year
            if "combo" in key:    
                #pull out year and add 0.2 so it appears at the end of that year
                datestr = key_split[1]+'.2'
            else:
                datestr = '20'+key_split[1][-2:]
            
            #add to datestring and dataset lists
            datestrs.append(datestr)
            dset_strs.append(key)
    #sort datasets according to date strings (just puts combos in right place)
    sorted_dset_strs = [x for _, x in sorted(zip(datestrs, dset_strs))]

    for i in np.arange((totaldsets)):
        thisd = d[sorted_dset_strs[i]]
        inner = outer[i:i+1].subgridspec(1, 4, wspace=0.01, hspace=0)
        if "combo" in sorted_dset_strs[i]:
            dset_pieces =sorted_dset_strs[i].split('_')
            datestr = dset_pieces[1]+' combo'
        else:
            datestr = thisd["date_obj"].strftime("%b %d %Y")
        
        
        if snr==True:
            if smt!=False:
                smtstr='sm'+str(smt)
            else:
                smtstr=''
            indivobj_fig(thisd["linesnr"][0,0,:,:],thisd["contsnr"][0,0,:,:],thisd["sdisnr"][0,0,:,:], thisd["scale"],
                thisd["prefix"]+'_SNR_'+smtstr,IWA=thisd["IWA"], smooth=smt, secondscale=1,secondscaleim=thisd["sdisnr2"][0,0,:,:],
                snr=True, title = thisd["obj"] + '\n'+ datestr, lims=lims, stampsz=stampsz, plspecs=thisd["plspecs"], plcand=thisd["plcand"], 
                returnfig=True, ax=(f1, inner, i, totaldsets), markctrl=int(thisd["ctrlrad"]), fwhm=int(thisd["fwhm"]))
        else:
            indivobj_fig(thisd["linecube"][0,:,:],thisd["contcube"][0,:,:],thisd["sdicube"][0,:,:], thisd["scale"],
                thisd["prefix"],IWA=thisd["IWA"], secondscale=1,secondscaleim=thisd["sdicube2"][0,:,:], 
                title = thisd["obj"] + '\n'+ datestr, lims=lims, stampsz=stampsz, plspecs=thisd["plspecs"], 
                plcand=thisd["plcand"], returnfig=True, ax=(f1,inner, i, totaldsets), 
                markctrl=int(thisd["ctrlrad"]), fwhm=int(thisd["fwhm"]))

    plt.savefig(str(d["outdir"])+str(d["objlist"])+str(d["theseparams"])+'.png')

    return ()

def plotdict_sdigrid(d, snr=False, lims=[-1,4], secondscale=False, stampsz=75, plcen=None, nperrow=2):
    """
    OPTIONAL INPUTS:
    plcen = a number equal to the index of the planet candidate you'd like to center the image stamp on
            default = None
    """
    #if d is string, open it
    if isinstance(d, str):
        d = pickle.load( open(d, "rb" ) )

    totaldsets = d["totaldsets"]

    nrows = int(np.ceil(totaldsets/nperrow))
    fig = plt.figure(figsize=(5*nperrow,nrows*5))
    gs = fig.add_gridspec(nrows, nperrow, hspace=0.25, wspace=0.4) # ,6)

    dkeys = list(d.keys())
    dset_strs = []

    platescale = 0.00795
    titlestyle=dict(size=16)
    fontprops = fm.FontProperties(size=16)

    for key in dkeys:
        if '_' in key:
            dset_strs.append(key)

    if stampsz=='adjust':
        adjstamp=True
    else:
        adjstamp=False

    for i in np.arange((totaldsets)):
        thisd = d[dset_strs[i]]

        #define left or righthand plot
        if i%nperrow==0:
            ax = fig.add_subplot(gs[int(i/nperrow),0])
            pltctr=1
        else:
            ax = fig.add_subplot(gs[int((i-1)/nperrow),pltctr])
            pltctr+=1

        if secondscale!=False:
            if snr==True:
                sdiim = thisd["sdisnr2"][0,0,:,:]
            else:
                sdiim = thisd["sdicube2"][0,:,:]
        else:
            if snr==True:
                sdiim = thisd["sdisnr"][0,0,:,:]
            else:
                sdiim = thisd["sdicube"][0,:,:]
        
        imsz = sdiim.shape[1]

        #mask IWA
        IWA=thisd["IWA"]
        IWAmask = ctrlmask(imsz, imsz, 0, IWA)
        sdiim*=IWAmask

        plspecs=thisd["plspecs"]
        plcand=thisd["plcand"]
        fwhm=thisd["fwhm"]
        ##if planets need to be marked, find their coordinates and define patches
        if plspecs!=False:
            plsep_y=[]
            plsep_x=[]
            pllabels=[]
            plseperr=[]

            #apply PA offset for VisAO (Balmer+2022) = 0.497+/-0.192 CCW
            paoff=-0.0497

            for j in np.arange(len(plspecs[0])):
                lb = plspecs[0][j]
                plsep = plspecs[1][j]
                plpa = plspecs[2][j]+paoff
                plsep_y.append(float(plsep)*np.sin((float(plpa)+90)*np.pi/180))
                plsep_x.append(float(plsep)*np.cos((float(plpa)+90)*np.pi/180))
                seperr = plspecs[3][j]
                if fwhm/2>seperr:
                    seperr=fwhm/2
                plseperr.append(seperr)
                pllabels.append(str(lb))

            #adjust stamp size according to planet separation if specified
            if adjstamp==True:
                xmax = np.nanmax(np.abs(plsep_x))
                ymax = np.nanmax(np.abs(plsep_y))

                if xmax>ymax:
                    stampsz=xmax*2+100
                else:
                    stampsz=ymax*2+100
        
        #this is dumb and I'm missing something silly in the algebra, but only centered
        #perfectly if odd numbered
        if stampsz%2==0:
            stampsz+=1

        ##COORDS RELATIVE TO STAMP CENTER
        stampcen = (stampsz - 1)/2.
        stampsz_asec = stampsz*platescale
        nticks = np.floor(stampsz_asec/2/0.25)
        ticklabels = np.arange(-1*nticks, nticks+1)*0.25
        ticklabels_str = [str(lab)+'\"' for lab in ticklabels]
        ticks = ticklabels/platescale + stampcen

        cen = (imsz - 1) / 2
        if plcen!=None:
            print('here', plsep_y, plsep_x, plcen)
            ylow = int(cen+plsep_y[plcen] - stampsz / 2 + 1)
            yhigh = int(cen+plsep_y[plcen] + stampsz / 2 )
            xlow = int(cen+plsep_x[plcen] - stampsz / 2 + 1)
            xhigh = int(cen+plsep_x[plcen] + stampsz / 2 )
            ##make ticklabels relative to central star
            xticklabels = np.arange(-1*nticks, nticks+1)*0.25+plsep_x[plcen]*platescale
            yticklabels = np.arange(-1*nticks, nticks+1)*0.25+plsep_y[plcen]*platescale
            xticklabels_str = ["{:4.2f}".format(lab)+'\"' for lab in xticklabels]
            yticklabels_str = ["{:4.2f}".format(lab)+'\"' for lab in yticklabels]
            markctrl=False
        else:
            ylow = int(cen - stampsz / 2 + 1)
            yhigh = int(cen + stampsz / 2 )
            xlow=ylow
            xhigh=yhigh
            markctrl=int(thisd["ctrlrad"])
            xticklabels_str=ticklabels_str
            yticklabels_str=ticklabels_str
        
        print(xlow,xhigh,ylow,yhigh)
        # set up tick labels according to parameter ranges
        plt.setp(ax, xticks=ticks, xticklabels=xticklabels_str, yticks=ticks, yticklabels=yticklabels_str)
     
        if snr==True:
            cbarlabel = 'SNR'
        else:
            cbarlabel = 'Intensity'
        #user defined limits
        if lims != False:
            minm = lims[0]
            sclmax = lims[1]
        #otherwise limits set by default
        else: 
            sclmax = np.nanstd(sdiim[ylow:yhigh, xlow:xhigh])*5
            minm = -1 * sclmax / 2

        im1 = ax.imshow(sdiim[ylow:yhigh, xlow:xhigh], vmin=minm, vmax=sclmax, origin='lower', cmap='magma')
     
        if stampsz<302:
            scalebar = AnchoredSizeBar(ax.transData, 0.1/platescale, '0.1\"', 'lower left', pad=0.1, color='white', 
                frameon=False, size_vertical=1, fontproperties=fontprops)
        else:
            scalebar = AnchoredSizeBar(ax.transData, 0.5/platescale, '0.5\"', 'lower left', pad=0.1, color='goldenrod', 
                frameon=False, size_vertical=1, fontproperties=fontprops)

        ax.add_artist(scalebar)

        print(dset_strs[i])
        if 'combo' in dset_strs[i]:
            title = dset_strs[i].replace('_', ' ')
        else:
            datestr = thisd["date_obj"].strftime("%b %d %Y")
            title = thisd["obj"] + '\n'+ datestr
        ax.set_title(title, **titlestyle)


        #add ctrlrad marker
        if markctrl>0:
            ctrlcirc=patches.Circle((stampcen,stampcen),radius=markctrl, fill=False, ec='white', lw=2, ls='--')
            ax.add_patch(ctrlcirc, )

        #add planet markers
        if plspecs!=False:
            for j in np.arange(len(pllabels)):
                #if candidate, make it a dashed line
                if plcand==True:
                    lsty=':'
                else:
                    lsty='-'
                
                if np.isfinite(plseperr[j]):
                    circrad = float(plseperr[j])
                else:
                    if stampsz<100:
                        circrad = 5
                    else:
                        circrad = 10

                if plcen!=None:
                    circ=patches.Circle((stampcen,stampcen),radius=circrad, fill=False, ec='cyan', lw=2, ls=lsty, alpha=0.5)
                else:
                    circ=patches.Circle((stampcen+plsep_x[j],stampcen+plsep_y[j]),radius=circrad, fill=False, ec='cyan', lw=2, ls=lsty, alpha=0.5)
                circ_label=pllabels[j]
                ax.add_patch(circ, )
                if plspecs[2][j]>180:
                    addx = [circrad+1]
                    align = 'left'
                else:
                    addx = [-1*(circrad+1)]
                    align = 'right'
                if plcen!=None:
                    labelposx = stampcen+addx[0]
                    labelposy = stampcen
                else:
                    labelposx = stampcen+plsep_x[j]+addx[0]
                    labelposy = stampcen+plsep_y[j] 
                ax.text(labelposx,labelposy,circ_label, color='white',fontsize=16, horizontalalignment=align, verticalalignment='center')
        
        if stampsz>302:
            resclr='goldenrod'
        else:
            resclr='white'
        resel = patches.Circle((stampsz-fwhm/2-5,fwhm/2+5 ),radius=fwhm/2, fill=True, color=resclr)
        ax.add_patch(resel, )
        cax = ax.inset_axes([1.04, 0.05, 0.05, 0.9],  transform=ax.transAxes)
        plt.colorbar(im1, ax=ax, cax=cax, orientation='vertical', label=cbarlabel) 


    #plt.tight_layout()

    plt.savefig(str(d["outdir"])+str(d["objlist"])+str(d["theseparams"])+'.png')

    return ()


def plotdict_ctrst(d, maxsep=1,to_mdot=None):
    """

    """
    #if d is string, open it
    if isinstance(d, str):
        d = pickle.load( open(d, "rb" ) )

    totaldsets = d["totaldsets"]

    colors = pl.cm.magma(np.linspace(0.1,0.9,totaldsets))
    platescale = 0.00795

    f,ax = plt.subplots(1, figsize=(8,6), dpi=750)
    dkeys = d.keys()
    dset_strs = []
    for key in dkeys:
        if key!="totaldsets":
            dset_strs.append(key)
    maxctrst=[]
    minctrst = []
    for i in np.arange((totaldsets)):
        thisd = d[dset_strs[i]]
        if "combo" not in thisd["prefix"]:
            seps = thisd["ccurve_seps"]*platescale
            ctrsts = thisd["ccurve_ctrst"]
            IWA = thisd["IWA"]*platescale

            #extrapolate to iwa
            #slope of innermost two points
            m = (np.log10(ctrsts[0])-np.log10(ctrsts[1]))/((seps[1]-seps[0]))
            IWA_y = 10**((seps[0]-IWA)*m+np.log10(ctrsts[0]))

            ax.plot([IWA, seps[0]],[IWA_y, ctrsts[0]],  '--', color=colors[i])
            ax.plot(seps, ctrsts, label = thisd["ccurve_lbl"], color=colors[i])
            maxctrst.append(IWA_y)
            nindlast = len(seps[seps<maxsep])
            minctrst.append(np.nanmin(ctrsts[:nindlast+2]))

    plt.yscale("log")
    plt.xlim(0, maxsep)
    plt.ylim(np.nanmin(minctrst),np.nanmax(maxctrst))
    plt.yscale("log")
    plt.xlabel("distance in arcseconds")
    plt.ylabel("contrast")

    if to_mdot!=None:
        #dist =
        #star_mag = 
        R = 1.6 #jupiter radii 
        units='Msun'
        law='aoyama'
        filtwid = 0.006
        xeropt = 2.339e-5
        mdot = am.contrast_to_MMdot(dist,star_mag,contr,R, M, units=units, law=law,filtwid=filtwid,zeropt=zeropt)


    ax.legend(frameon=False, fontsize=10)
    return()

def indivobj_fig(lineim, contim, sdiim, scale, prefix, title=False, secondscale=False, secondscaleim=False, IWA=0, outputdir='final_ims/', snr=False, stampsz=75, smooth=0, lims = False, plspecs=False, plcand=False, returnfig=False, ax=None, markctrl=None, fwhm=None):
    """
    creates a three panel figure with line, continuum and SDI images

    optional keywords:
    lims = tuple in the form (min, max) to set colorbar limits

    """
    if ax==False:
        f, (ax1,ax2,ax3,ax4) = plt.subplots(1, nplots, sharey=True)
        num=0
        totaldsets=1
    else:
        ax1 = plt.Subplot(ax[0], ax[1][0])
        ax2 = plt.Subplot(ax[0], ax[1][1])
        ax3 = plt.Subplot(ax[0], ax[1][2])  
        ax4 = plt.Subplot(ax[0], ax[1][3])  
        num = ax[2]
        totaldsets = ax[3]

    axs = (ax1,ax2,ax3,ax4)

    if secondscale==True:
        nplots = 4
    else:
        nplots = 3

    platescale = 0.00795
    imsz = lineim.shape[1]

    #mask IWA
    IWAmask = ctrlmask(imsz, imsz, 0, IWA)
    lineim*=IWAmask
    contim*=IWAmask
    sdiim*=IWAmask
    if secondscale!=False:
        secondscaleim*=IWAmask

    ##COORDS RELATIVE TO STAMP CENTER
    stampcen = (stampsz - 1)/2.
    stampsz_asec = stampsz*platescale
    nticks = np.floor(stampsz_asec/2/0.25)
    ticklabels = np.arange(-1*nticks, nticks+1)*0.25
    ticklabels_str = [str(lab)+'\"' for lab in ticklabels]
    ticks = ticklabels/platescale + stampcen

    cen = (imsz - 1) / 2
    low = int(cen - stampsz / 2 + 1)
    high = int(cen + stampsz / 2 )

    if smooth != 0:
        gauss = conv.Gaussian2DKernel(x_stddev=smooth)
        lineim = conv.convolve(lineim, gauss, preserve_nan=True)
        contim = conv.convolve(contim, gauss, preserve_nan=True)
        sdiim = conv.convolve(sdiim, gauss, preserve_nan=True)
        if secondscale!=False:
            secondscaleim = conv.convolve(secondscaleim, gauss, preserve_nan=True)

    # set up tick labels according to parameter ranges
    if secondscale!=False:
        plt.setp(axs, xticks=ticks, xticklabels=ticklabels_str, yticks=ticks, yticklabels=ticklabels_str)
    else:
        plt.setp(axs, xticks=ticks, xticklabels=ticklabels_str, yticks=ticks, yticklabels=ticklabels_str)

    if snr==True:
        cbarlabel = 'SNR'
    else:
        cbarlabel = 'Intensity'
    #user defined limits
    if lims != False:
        minm = lims[0]
        linemax = lims[1]
    #otherwise limits set by default
    else: 
        linemax = np.nanstd(lineim[low:high, low:high])*5
        minm = -1 * linemax / 2

    ##if planets need to be marked, find their coordinates and define patches
    if plspecs!=False:
        plsep_y=[]
        plsep_x=[]
        pllabels=[]
        plseperr=[]

        #apply PA offset for VisAO (Balmer+2022) = 0.497+/-0.192 CCW
        paoff=-0.0497

        for i in np.arange(len(plspecs[0])):
            lb = plspecs[0][i]
            plsep = plspecs[1][i]
            plpa = plspecs[2][i]+paoff
            plsep_y.append(float(plsep)*np.sin((float(plpa)+90)*np.pi/180))
            plsep_x.append(float(plsep)*np.cos((float(plpa)+90)*np.pi/180))
            seperr = plspecs[3][i]
            if fwhm/2>seperr:
                seperr=fwhm/2
            plseperr.append(seperr)
            pllabels.append(str(lb))

    titlestyle=dict(size=16)

    #ax1 = ax[0]
    im1 = ax1.imshow(lineim[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
    im2 = ax2.imshow(contim[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
    im3 = ax3.imshow(sdiim[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
 
    fontprops = fm.FontProperties(size=16)
    if stampsz<250:
        scalebar = AnchoredSizeBar(ax1.transData, 0.1/platescale, '0.1\"', 'lower left', pad=0.1, color='white', 
            frameon=False, size_vertical=1, fontproperties=fontprops)
    else:
        scalebar = AnchoredSizeBar(ax1.transData, 0.5/platescale, '0.5\"', 'lower left', pad=0.1, color='goldenrod', 
            frameon=False, size_vertical=1, fontproperties=fontprops)

    ax1.add_artist(scalebar)

    if num==0:
        ax1.set_title(r'KLIP-ed H$\alpha$ Image', **titlestyle)
        ax2.set_title(r'KLIP-ed Continuum Image',**titlestyle)
        ax3.set_title(r'ASDI Image (H$\alpha$-scale$\times$Cont)',**titlestyle)
    
    ax2.axes.get_yaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)

    if stampsz<300:
        labelstyle = dict(size=16, color='white', weight='bold')
    else:
        labelstyle = dict(size=16, color='goldenrod', weight='bold')        
    ax3.text(0.58,0.93, 'scale='+'{:.2f}'.format(scale), transform=ax3.transAxes,**labelstyle)

    if secondscale!=False:
        im4 = ax4.imshow(secondscaleim[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
        if num==0:
            ax4.set_title(r'ASDI Image (H$\alpha$-Cont)',**titlestyle)
        ax4.text(0.58,0.93, 'scale='+'{:.2f}'.format(secondscale), transform=ax4.transAxes,**labelstyle)
        ax4.axes.get_yaxis().set_visible(False)
        divider = make_axes_locatable(ax4)
 
    if num!=totaldsets-1:
        for a in axs:
            a.axes.get_xaxis().set_visible(False)

    #add ctrlrad marker
    if markctrl>0:
        for a in axs:
            ctrlcirc=patches.Circle((stampcen,stampcen),radius=markctrl, fill=False, ec='white', lw=2, ls='--')
            a.add_patch(ctrlcirc, )

    #add planet markers
    if plspecs!=False:
        for i in np.arange(len(pllabels)):
            #if candidate, make it a dashed line
            if plcand==True:
                lsty=':'
            else:
                lsty='-'
            for a in axs:
                if np.isfinite(plseperr[i]):
                    circrad = float(plseperr[i])
                else:
                    circrad = 5
                circ=patches.Circle((stampcen+plsep_x[i],stampcen+plsep_y[i]),radius=circrad, fill=False, ec='cyan', lw=2, ls=lsty)
                circ_label=pllabels[i]
                a.add_patch(circ, )
                if stampsz<250:
                    addx = [circrad+1 if (plspecs[2][i]>180 or plspecs[2][i]<10) else -3*circrad]
                else:
                    addx = [circrad+10 if plspecs[2][i]>180 else -1*(circrad+15)]
                labelposx = stampcen+plsep_x[i]+addx[0]
                labelposy = stampcen+plsep_y[i]
                a.text(labelposx,labelposy,circ_label, color='white',fontsize=16)
    if stampsz>300:
        resclr ='goldenrod' 
    else:
        resclr='white'
    resel = patches.Circle((stampsz-fwhm/2-5,fwhm/2+5 ),radius=fwhm/2, fill=True, color=resclr)
    ax4.add_patch(resel, )

    if title!=False:
        #plt.suptitle(title, size=22) 
        ax1.text(-1*stampsz/3,stampsz/3.5, title, rotation=90, fontsize=20, weight='bold', multialignment='center')  

    if ax!=False:
        ax[0].add_subplot(ax1)
        ax[0].add_subplot(ax2)
        ax[0].add_subplot(ax3)
        ax[0].add_subplot(ax4)
        cax = ax4.inset_axes([1.04, 0.05, 0.05, 0.9],  transform=ax4.transAxes)
        ax[0].colorbar(im4, ax=ax4, cax=cax, orientation='vertical', label=cbarlabel)    
    
        #ax[0].colorbar(im4, shrink=0.5, orientation='vertical', label=cbarlabel)

    plt.tight_layout() 
    
    if returnfig==True:
        print('done')
    else:
        plt.savefig(outputdir+prefix+'.png')
        plt.show()
        plt.clf()
        return()

def contrast_klcompare(data_str, ha_ctrsts, cont_ctrsts, KLlist, IWA, zone_boundaries, outdir = 'final_ims/'):
    
    """
    """
    #check KL list
    kldim = ha_ctrsts.shape[0]
    if kldim != len(KLlist):
        print("KL list is wrong shape relative to contrast files")

    platescale = 0.00795
    zone_boundaries_arcsec = [x*platescale for x in zone_boundaries]

    #make figure
    plt.figure(figsize=(7,4), dpi=750)
    for klctr in np.arange(len(KLlist)):
        j = klctr / len(KLlist)
        ha_contrast_seps = ha_ctrsts[klctr,0,:]
        ha_corrected_curve = ha_ctrsts[klctr,2,:]
        ha_contrasts = ha_ctrsts[klctr,3:,:]     
        plt.plot(ha_contrast_seps * platescale, ha_corrected_curve, label="H-alpha"+str(KLlist[klctr])+' KL', color=cm.plasma(j))

    floor = np.log10(np.nanmin(ha_contrasts)) - 0.2
    plt.yscale("log")
    plt.title(data_str.replace('_',' '))
    plt.xlim(0, 1)
    plt.ylim(10 ** floor, 1e-1)
    plt.yscale("log")
    plt.xlim(0, 1)
    plt.xlabel("distance in arcseconds")
    plt.ylabel("contrast")
    plt.xlabel("distance in arcseconds")
    plt.ylabel("contrast")
    if IWA > 0:
        plt.plot((IWA * platescale, IWA * platescale), (1e-5, 1e-1), 'k', label='IWA')
    for bd in zone_boundaries_arcsec:
        if bd == zone_boundaries_arcsec[0]:
            plt.plot((bd, bd), (0, 1), 'k--', lw=1, label='zone boundary')
        else:
            plt.plot((bd, bd), (0, 1), 'k--', lw=1)
    plt.legend()
    # write out in data directory
    plt.savefig(outdir+data_str+'Hacontrasts_byKLmode.jpg')
    plt.show()
    plt.clf()

    plt.figure(figsize=(7,4), dpi=750)
    for klctr in np.arange(len(KLlist)):
        j = klctr / len(KLlist)
        cont_contrast_seps = cont_ctrsts[klctr,0,:]
        cont_corrected_curve = cont_ctrsts[klctr,2,:]
        cont_contrasts = cont_ctrsts[klctr,3:,:]        
        plt.plot(cont_contrast_seps * platescale, cont_corrected_curve, label='Continuum'+str(KLlist[klctr])+' KL', color=cm.plasma(j))

    plt.yscale("log")
    plt.title(data_str.replace('_',' '))
    plt.xlim(0, 1)
    plt.ylim(10 ** floor, 1e-1)
    plt.yscale("log")
    plt.xlim(0, 1)
    plt.xlabel("distance in arcseconds")
    plt.ylabel("contrast")
    plt.xlabel("distance in arcseconds")
    plt.ylabel("contrast")
    if IWA > 0:
        plt.plot((IWA * platescale, IWA * platescale), (1e-5, 1e-1), 'k', label='IWA')
    for bd in zone_boundaries_arcsec:
        if bd == zone_boundaries_arcsec[0]:
            plt.plot((bd, bd), (0, 1), 'k--', lw=1, label='zone boundary')
        else:
            plt.plot((bd, bd), (0, 1), 'k--', lw=1)
    plt.legend()
    # write out in data directory
    plt.savefig(outdir+data_str+'Contcontrasts_byKLmode.jpg')
    plt.show()
    plt.clf()


def final_contrast_fig(data_str, ha_ctrsts, cont_ctrsts, IWA, zone_boundaries, outdir = 'final_ims/', point=False):
    
    """
    point = overplot a detection point with a label. Give in form [sep (in arcsec), contrast, label]
    """

    kldim = ha_ctrsts.shape[0]

    ha_contrast_seps = ha_ctrsts[0,0,:]
    ha_corrected_curve = ha_ctrsts[0,2,:]
    ha_contrasts = ha_ctrsts[0,3:,:]
    cont_contrast_seps = cont_ctrsts[0,0,:]
    cont_corrected_curve = cont_ctrsts[0,2,:]
    cont_contrasts = cont_ctrsts[0,3:,:]

    min_ctrst_ha = np.nanmin(ha_contrasts, axis=0)
    max_ctrst_ha = np.nanmax(ha_contrasts, axis=0)

    min_ctrst_cont = np.nanmin(cont_contrasts, axis=0)
    max_ctrst_cont = np.nanmax(cont_contrasts, axis=0)

    platescale = 0.00795
    zone_boundaries_arcsec = [x*platescale for x in zone_boundaries]

    #make figure
    plt.figure(figsize=(7,4), dpi=750)
    plt.plot(ha_contrast_seps * platescale, ha_corrected_curve, label="H-alpha", color='blue')
    plt.fill_between(ha_contrast_seps * platescale, min_ctrst_ha, max_ctrst_ha, color='blue', alpha=0.5)
    plt.plot(cont_contrast_seps * platescale, cont_corrected_curve, label='Continuum', color='red')
    plt.fill_between(cont_contrast_seps * platescale, min_ctrst_cont, max_ctrst_cont, color='red', alpha=0.5)
    floor = np.log10(np.nanmin(ha_contrasts)) - 0.2
    plt.yscale("log")
    plt.title(data_str.replace('_',' '))
    plt.xlim(0, 1)
    plt.ylim(10 ** floor, 1e-1)
    plt.yscale("log")
    plt.xlim(0, 1)
    plt.xlabel("distance in arcseconds")
    plt.ylabel("contrast")
    plt.xlabel("distance in arcseconds")
    plt.ylabel("contrast")
    if IWA > 0:
        plt.plot((IWA * platescale, IWA * platescale), (1e-5, 1e-1), 'k', label='IWA')
    for bd in zone_boundaries_arcsec:
        if bd == zone_boundaries_arcsec[0]:
            plt.plot((bd, bd), (0, 1), 'k--', lw=1, label='zone boundary')
        else:
            plt.plot((bd, bd), (0, 1), 'k--', lw=1)
    if point != False:
        plt.plot(point[0],point[1], 'c*', label=point[3])
    plt.legend()
    # write out in data directory
    plt.savefig(outdir+data_str+'_contrast_curve.jpg')
    plt.show()
    plt.clf()
