from astropy.io import fits
import os
import pandas as pd
import numpy as np
import gaussfitter
import astropy.convolution as conv
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
import matplotlib.patches as patches
import SNRMap_new as snr
from importlib import reload
reload(snr)

def SliceCube(imfile, rotfile, dir='./', slicedir='sliced/'):
    """
    Slices a 3D cube into individual 2D images numbered sequentially
    dir: directory containing the 3D image
    imfile: name of the 3D image
    rotfile: name of the cube containing the rotoff values needed to rotate images
    slicedir: name of the directory to write the sliced images into
    """
    origdir=os.getcwd()
    os.chdir(dir)

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
    if not os.path.exists(slicedir):
        os.makedirs(slicedir)
    os.chdir(slicedir)

    #loop to save 2D images into the sliced directory
    for z in range(nims):
        single = cube[z,:,:]
        head["rotoff"]=rotoffs[z]
        fits.writeto("sliced_"+str(z+1)+".fits", single, head, overwrite=True)

    print("Done creating individual images for KLIP. These live in the", slicedir, 'directory')
    os.chdir(origdir)

    return

def get_cuts_df(dfname):
    # define dataframe if doesn't already exist
    try:
        df=pd.read_csv(dfname)
        print("found existing data frame") #with the following contents: \n", df)
    except:
        print("creating new data frame")
        #if none found, make one
        df_cols = ['Dataset', 'wl', 'pctcut', 'nims', 'minpeak', 'medfwhm','rotrange']
        # make a blank pandas dataframe
        df = pd.DataFrame({})
        for name in df_cols:
            df[name] = []
    return(df)

def peak_cut(subdir, wl, dfname='cuts.csv', imstring='_clip451_flat_reg_nocosmics',
             debug=False, ghost=False, pctcuts=[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90], electrons=False):
    """
    PURPOSE
    Makes new images and rotoff cubes containing only images that have peaks greater
    than or equal to the percentile threshhold(s) specified with the pctcuts keyword.
    Writes files to disk with appropriate suffixes and slices them for KLIP. Records
    peaks, nimages, etc. in dataframe.

    REQUIRED INPUTS
    subdir = directory holding images that you want to reduce
    wl = "Line" or "Cont" depending on which you are reducing
    df = pandas dataframe to store output

    OPTIONAL INPUTS
    imstr = name string for 3D (x,y,t) image holding data
    ghost = if set to True, will do cuts on ghost peak and write ghost and estimate of star
            peak to header
    pctcuts = a list of percentile threshholds
    electrons = keyword to indicate that units are electrons (just used for sanity check of peak values)

    OUTPUTS
    df = populated pandas dataframes

    written by Kate Follette June 2019
    """

    df = get_cuts_df(subdir+dfname)

    if wl == 'Line':
        rotstring = 'rotoff_noLinecosmics'
    if wl == 'Cont':
        rotstring = 'rotoff_noContcosmics'

    origdir = os.getcwd()
    os.chdir(subdir)

    #read in image, header, and rotoffs

    if os.path.exists(rotstring + '.fits'):
        rotoffs = fits.getdata(rotstring + '.fits')
        imcube = fits.getdata(wl + imstring + '.fits')
        head = fits.getheader(wl + imstring + '.fits')
    else:
        rotoffs = fits.getdata(rotstring + '_0pctcut.fits')
        imcube = fits.getdata(wl + imstring + '_0pctcut.fits')
        head = fits.getheader(wl + imstring + '_0pctcut.fits')

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
    fwhmlist = []

    #size of image stamp
    width = 61
    stmpsz = int((width - 1) / 2)

    # loop through images fitting peaks and storing them
    for i in np.arange(0, nims):
        im = imcube[i, ycen - stmpsz-1:ycen + stmpsz, xcen - stmpsz-1:xcen + stmpsz]
        p = gaussfitter.moments(im, circle=True, rotate=False, vheight=False, estimator=np.nanmedian)
        fullfit.append(p)
        peaklist.append(p[0])
        fwhmlist.append(p[3]*2.)
        if electrons == True:
            # if units are electrons, raise the peak threshhold
            if (p[0] < 0) or (p[0] > 50000):
                print("warning: unphysical peak value of", p[0], 'for image', i + 1)
        if electrons == False:
            if (p[0] < 0) or (p[0] > 17000):
                print("warning: unphysical peak value of", p[0], 'for image', i + 1)

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
        imname = wl + imstring + '_' + str(pctcuts[j]) + 'pctcut.fits'
        hdu.writeto(imname, overwrite=True)
        rotnew = rotoffs[goodims]

        #stupid catch for rot angles that pass through zero
        if rotnew[0] - rotnew[-1] >= 180:
            rotnew[-1]+=360
            rotrange = (rotnew[0] - rotnew[-1]) / ((rotoffs[0]) - (rotoffs[-1]+360)) * 100
        else:
            rotrange = (rotnew[0] - rotnew[-1]) / ((rotoffs[0]) - (rotoffs[-1])) * 100
        hdu2 = fits.PrimaryHDU(rotnew)
        rotname = rotstring + '_' + str(pctcuts[j]) + 'pctcut.fits'
        hdu2.writeto(rotname, overwrite=True)

        if debug == True:
            print(pctcuts[j], "% cut")
            print("Peak > ", cut)
            print("check:", len(goodims) / nims * 100, "pct of images")
            print("shape of new data cube is", newim.shape)
            print(rotrange, "pct of rotation")
            print("slicing cube", imname)

        # slice the newly created file in preparation for KLIP
        outdir = wl + '_' + str(pctcuts[j]) + 'pctcut_sliced/'
        SliceCube(imname, rotname, slicedir=outdir)

        # and add star peak to headers of sliced files
        if ghost == True:
            pykh.addstarpeak(outdir, debug=debug, mask=True, ghost=True, wl=wl)
        else:
            pykh.addstarpeak(outdir, debug=debug, mask=True)

        # record relevant quantities for this data cut to dataframe
        datalist = pd.DataFrame([[subdir, wl, pctcuts[j], new_nims, cuts[j], np.nanmedian(goodfwhm),
                                  rotrange]], columns=['Dataset', 'wl', 'pctcut', 'nims', 'minpeak', 'medfwhm','rotrange'])
        df = df.append(datalist)
        j += 1

    os.chdir(origdir)
    df.to_csv(subdir+dfname, index=False)

    return


def make_prefix(subdir, wl, cut):
    if subdir != './':
        prefix = subdir.replace('/', '_') + wl + '_' + str(cut) + 'pctcut'
    else:
        prefix = wl + '_' + str(cut) + 'pctcut'
    return(prefix)

def compute_thrpt(subdir, wl, cut, numann=3, movm=4, KLlist=[10], IWA=0,
                  contrast=1e-2, theta=0., clockang=85, debug=False, record_seps=[0.1, 0.25, 0.5, 0.75, 1.0],
                  ghost=False, savefig=False, iterations=3,dfname='cuts.csv'):
    """
    PURPOSE
    Injects false planets separated by the measured fwhm of the dataset in radius and the clockang parameter
    in azimuth and calculates the  throughput.

    REQUIRED INPUTS
    subdir = directory holding images that you want to reduce
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

    #read in or make data frame
    df = get_cuts_df(subdir+dfname)

    #add contrast to df
    if "ctrst_fkpl" not in df:
        df["ctrst_fkpl"] = np.nan
        print("creating column", "ctrst_fkpl")
    df["ctrst_fkpl"]=contrast

    # set up directories and naming
    origdir = os.getcwd()
    os.chdir(subdir)
    outputdir = 'contrastcurves/'
    # if contrast curve directory doesn't already exist, create it
    if os.path.exists(outputdir) == False:
        os.mkdir(outputdir)
    prefix = make_prefix(subdir, wl, cut)

    ###load data
    filelist = glob.glob(wl + '_' + str(cut) + 'pctcut_sliced' + "/sliced*.fits")

    # sorts file list so it's sequential
    filelist.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    ### pull the values of the star peak from the headers
    starpeak = []
    for i in np.arange(len(filelist)):
        head = fits.getheader(filelist[i])
        starpeak.append(head["STARPEAK"])

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
    first = IWA + fwhm
    n_planets = (dataset.input.shape[1] / 2. - first) / fwhm
    n_planets = int(n_planets)
    n_planets -= 1
    if debug == True:
        print('I will inject ', n_planets, 'planets')

    # calculates separations (in pixels) where planets will be injected
    thrpt_seps = []
    sep = first
    for i in np.arange(n_planets):
        thrpt_seps.append(sep)
        sep += fwhm

    if debug == True:
        print("injecting planets at the following radii:", thrpt_seps)

    size = dataset.input.shape
    nims = size[0]
    xcen = int((size[1] - 1) / 2)
    ycen = int((size[2] - 1) / 2)

    #save a version of the input dataset without any fakes
    dataset_copy = np.copy(dataset.input)
    imsz = dataset.input.shape[1]

    if ghost == False:
        # multiply starstamp image by contrast to create your false planets
        to_inject = contrast * dataset.input

    thrpts = np.zeros((iterations, len(thrpt_seps)))

    prefix_fakes = prefix +'_initPA'+str(theta)+'_CA'+str(clockang)+'_ctrst'+str(contrast)+'_'+str(iterations)+'FAKES'

    dataset_prefix = prefix_fakes +  '_a' + str(numann) + 'm' + str(movm) + 'iwa'+str(IWA)

    tpt_fname = outputdir + prefix_fakes +  '_a' + str(numann) + 'm' + str(movm) + 'iwa'+str(IWA) + '_throughputs.fits'

    calcthrpt=True
    
    if os.path.exists(tpt_fname):
        print("throughput file", tpt_fname, "already exists. reading in values.")
        calcthrpt=False
        thrpt_cube = fits.getdata(tpt_fname)
        thrpt_seps = thrpt_cube[0]
        thrpt_avgs = thrpt_cube[1]
        thrpts = thrpt_cube[2:]

    # full sequence of planet injection, recovery and throughput calculation as many times as the iterations keyword, with
    # starting planet locations clocked by 75 degrees each cycle
    for iter in np.arange(iterations):
        #clock starting planet location by 75deg each cycle
        theta=75.*iter
        
        #make prefixes for output files
        pfx=prefix_fakes + '_a' + str(numann) + 'm' + str(movm) + 'iwa'+str(IWA) + '_set' + str(iter+1)

        #check whether fake files or throughputs with these parameters have already been calculated
        runfakes=True

        if os.path.exists(outputdir+pfx+'-KLmodes-all.fits'):
            print("file with prefix", pfx, "already exists")
            runfakes=False

        #special case where pipeline broke in between the two steps - needs whole cube in memory to do throughput calc.
        if runfakes == False and calcthrpt == True:
            runfakes = True
            print("but I can't find a throughput file so I am regenerating it.")

        if runfakes==True:

            #pull the clean input data every time
            dataset.input = np.copy(dataset_copy)

            # replace any nans in image with zeros (edges, padded) for fake injection
            dataset.input[np.isnan(dataset.input) == True] = 0.

            #inject planets in raw images
            for sep in thrpt_seps:
                if ghost == True:
                    fakes.inject_planet(dataset.input, dataset.centers, np.repeat(contrast, nims)[0], dataset.wcs, sep, theta,
                                    fwhm=fwhm)
                else:
                    fakes.inject_planet(dataset.input, dataset.centers, to_inject, dataset.wcs, sep, theta, fwhm=fwhm,
                                    stampsize=imsz)  # , thetas=thetas)
                theta += clockang

            # put NaNs back
            dataset.input[np.isnan(dataset_copy) == True] = np.nan

            # KLIP dataset with fake planets. Highpass filter here.
            parallelized.klip_dataset(dataset, outputdir=outputdir, fileprefix=pfx, algo='klip', annuli=numann,
                                      subsections=1, movement=movm, numbasis=KLlist, calibrate_flux=False, mode="ADI",
                                      highpass=True, save_aligned=False, time_collapse='median')

        # reset initial theta for recovery loop
        theta = 75.*iter

        # create a list that will store throughputs for this iteration
        thrpt_list = []

        # loop counter
        i = 0

        #calculate a few fixed quantities
        annspacing = (imsz / 2. - dataset.IWA) / numann
        zone_boundaries = np.arange(1, numann) * annspacing + dataset.IWA

        if calcthrpt==True:

            nimages = dataset.output.shape[1]
            fake_fluxes = np.zeros((n_planets, nimages))

            #planet recovery loop
            for sep in thrpt_seps:
                # thetas for retrieve planet call are measured from +x axis, so should always add 90 to pa for thetas keyword
                # dataset shape is [KL modes, n images, 1(wl dim), xdim, ydim]
                fake_flux = fakes.retrieve_planet_flux(dataset.output[0, :, 0, :, :], dataset.centers, dataset.wcs, sep, theta,
                                                       searchrad=8, guessfwhm=fwhm,
                                                       guesspeak=contrast * np.median(dataset.star_flux),
                                                       thetas=np.repeat(theta + 90, dataset.output.shape[1]), refinefit=True)
                fake_fluxes[i,:] = fake_flux
                newthpt = np.nanmedian(fake_flux / (contrast))
                if (newthpt > 1) or (newthpt < 0):
                    print('invalid throughput value of', newthpt, 'for planet at separation', sep, '. Replacing with NaN')
                    newthpt = np.nan
               
                #check for bad value - throughput should increase with distance. decrease suggests inner point may not have been a true recovery
                if len(thrpt_list) > 1:
                    if newthpt < thrpt_list[-1]:
                        if abs(newthpt - thrpt_list[-1]) > 0.1:
                            print('suspicious trend in throughput. Value at', sep, 'is', newthpt, 'but was', thrpt_list[-1], 'in previous iteration.')
                            #if planet is within 10 pixels of zone boundary, relax the requirement that the trend has to be upward
                            if len(zone_boundaries) > 0:
                                if abs(min(zone_boundaries - sep)) < 10:
                                    print(min(zone_boundaries - sep), 'pixels from zone boundary. Keeping.')
                                else:
                                    print('planet is', min(zone_boundaries - sep), 'pixels from nearest zone boundary. Setting to NaN.')
                                    newthpt=np.nan

                thrpt_list.append(newthpt)
                
                #if debug == True:
                    #print("fake planet at ", sep, " pixels has a mean throughput of", np.nanmean(fake_flux / (contrast)),
                      #"median of", np.nanmedian(fake_flux / (contrast)), " and stdev of ",
                      #np.nanstd(fake_flux / (contrast)))
                theta += clockang
                i += 1
            thrpts[iter, :] = thrpt_list

    thrpt_avgs=np.nanmean(thrpts,axis=0)

    #up to 10 iterations. will break if do more. 
    cx = ['rx','gx','kx','cx','mx','rs','gs','ks','cs','ms']
    #if on last iteration, make plot and save
    if (savefig == True) and (iter==iterations-1):
        #plot the individual points
        plt.figure(figsize=(10, 5), dpi=750)
        for iter in np.arange(iterations):
            plt.plot(thrpt_seps, thrpts[iter,:], cx[iter], label="set"+str(iter+1))
        # plot the throughput averages (should all be <1 and should increase outward until they hit a zone boundary)
        plt.plot(thrpt_seps, thrpt_avgs, 'bo', label="average")
        for bd in zone_boundaries:
            plt.plot((bd, bd), (0, 1), '--', label='zone boundary')
        plt.xlabel("separation in pixels")
        plt.ylabel("throughput")
        plt.title(dataset_prefix+" Throughput")
        plt.legend()
        plt.savefig(outputdir + dataset_prefix + '_throughput.jpg')
        plt.show()
        plt.clf()

    # locations to record throughputs in table
    platescale = 0.0078513

    # blank array to store interpolated throughputs
    thrpt_table_vals = []
    # loop through locations, find throughput and record
    for loc in record_seps:
        # interpolate over throughputs to find value at loc
        if loc < dataset.input.shape[1]/2*platescale:
            tpt = np.interp(loc / platescale, thrpt_seps, thrpt_avgs)
        else:
            tpt = np.nan
        thrpt_table_vals.append(tpt)
        # if df keyword is set, records in dataframe.
        if len(df) > 1:
            colname = "tpt_" + str(loc)
            # if column doesn't already exist, create and fill with nans
            if colname not in df:
                df[colname] = np.nan
                print("creating column", colname)
            df[colname].loc[(df.Dataset == subdir) & (df.pctcut == cut)] = tpt

    thrpt_out=np.zeros((2+iterations, len(thrpt_seps)))
    thrpt_out[0]=thrpt_seps
    thrpt_out[1]=thrpt_avgs
    thrpt_out[2:]=thrpts
    np.array([thrpt_seps,thrpt_avgs,thrpts])
    print(thrpt_out.shape)

    head["NPLANET"]=n_planets
    head["CTRST"]=contrast
    head["ITERS"]=iterations
    head["PASTART"]=0
    head["CLOCKANG"]=clockang
    head["ROW1"]="separations (pix)"
    head["ROW2"]="average throughput"
    head["OTHROWS"]="individual throughputs"
        
    fits.writeto(tpt_fname, thrpt_out, header=head, overwrite=True)

    os.chdir(origdir)

    df.to_csv(subdir+dfname, index=False)

    return (thrpt_out, zone_boundaries, df, dataset_prefix)


def make_contrast_curve(subdir, wl, cut, thrpt_out, dataset_prefix, numann=3,
                        movm=4, KLlist=[10], IWA=0, dfname='cuts.csv', record_seps=[0.1, 0.25, 0.5, 0.75, 1.0], 
                        savefig=False, debug=False):
    """
    PURPOSE
    calculates raw contrast by running KLIP on images and then corrects for throughput
    to generate a calibrated contrast curve.

    REQUIRED INPUTS
    subdir = directory holding images that you want to reduce
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
    df = get_cuts_df(subdir+dfname)

    iterations = len(thrpt_out[:,0])-2
    platescale = 0.0078513

    # set up directories and naming
    prefix = dataset_prefix
    print(prefix)
    origdir = os.getcwd()
    os.chdir(subdir)
    outputdir = 'contrastcurves/'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    #pull throughput info 
    thrpt_seps = thrpt_out[0]
    thrpt_avgs = thrpt_out[1]
    thrpts = thrpt_out[2:]

    if os.path.exists(outputdir+prefix+'-KLmodes-all.fits'):
        print('KLIPed image without fakes already exists for these parameters')
        klcube = fits.getdata(outputdir+prefix+'-KLmodes-all.fits')
        klheader = fits.getheader(outputdir+prefix+'-KLmodes-all.fits')

    else:
        # specify directory and read in data again. This time we will NOT inject fakes
        filelist = glob.glob(wl + '_' + str(cut) + 'pctcut_sliced' + "/sliced*.fits")
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
        parallelized.klip_dataset(dataset, outputdir=outputdir, fileprefix=prefix, annuli=numann, subsections=1,
                                  algo='klip', movement=movm, numbasis=KLlist, calibrate_flux=False,
                                  mode="ADI", highpass=True, time_collapse='median')

        ##pull some needed info from headers
        kl_hdulist = fits.open("{out}/{pre}-KLmodes-all.fits".format(out=outputdir, pre=prefix))
        klcube = kl_hdulist[1].data
        klheader = kl_hdulist[0].header
    
    #read in fwhm
    dataset_fwhm = klheader["0PCTFWHM"]

    #raw contrast is independent of fake injection
    rawc_prefix = make_prefix(subdir, wl, cut)

    #set OWA
    klim = klcube[0, :, :]
    if klim.shape[1]/2 <= 140:
        OWA = klim.shape[1]/2
    else:
        OWA = 140  # a little beyond 1" is the furthest out we will go for computing contrast. Can change if needed.


    #check whether has already been computed
    if os.path.exists(outputdir+rawc_prefix+'_rawcontrast.fits'):
        ctrst = fits.getdata(outputdir+rawc_prefix+'_rawcontrast.fits')
        contrast_seps = ctrst[0,:]
        contrast = ctrst[1,:]
    
    #if not, generate raw contrast and write out as .jpg and .fits
    else:
        # first number is KL mode index of cube. Don't change unless you ran multiple KL modes.
        
        dataset_center = [klheader['PSFCENTX'], klheader['PSFCENTY']]
        contrast_seps, contrast = klip.meas_contrast(klim, IWA, OWA, dataset_fwhm,
                                                     center=dataset_center, low_pass_filter=False)

        if savefig == True:
            plt.figure(figsize=(10, 5), dpi=750)
            imsz = klim.shape[1]
            annspacing = (imsz / 2. - IWA) / numann
            zone_boundaries = np.arange(1, numann) * annspacing + IWA
            plt.plot(contrast_seps * platescale, contrast)
            plt.plot(contrast_seps * platescale, contrast, 'bo')
            plt.yscale("log")
            plt.ylim(np.nanmin(contrast), 1e-1)
            plt.xlim(0, platescale * OWA)
            plt.xlabel("distance in arcseconds")
            plt.ylabel("contrast")
            if IWA > 0:
                plt.plot((IWA * platescale, IWA * platescale), (1e-5, 1e-1), label='IWA')
            for bd in zone_boundaries * platescale:
                if bd < OWA * platescale:
                    plt.plot((bd, bd), (0, 1), '--', label='zone boundary')
            plt.legend()
            plt.title(rawc_prefix+" Raw Contrast")
            plt.savefig(outputdir + rawc_prefix + '_rawcontrast.jpg')
            plt.show()
            plt.clf()  # clear figure

        contrast_out=np.zeros((2, len(contrast_seps)))
        contrast_out[0,:]=contrast_seps
        contrast_out[1,:]=contrast
        fits.writeto(outputdir+rawc_prefix+'_rawcontrast.fits', contrast_out, overwrite=True)

    ## corrected contrast figure
    corrected_contrast_curve = np.copy(contrast)

    interp_thrpts = []

    for i, sep in enumerate(contrast_seps):
        thrpt_interp = np.interp(sep, thrpt_seps, thrpt_avgs)
        if debug == True:
            closest_throughput_index = np.argmin(np.abs(thrpt_seps - sep))
            print('for separation', sep, " closest throughput is at separation ", thrpt_seps[closest_throughput_index])
            print('interpolated throughput is', thrpt_interp, 'for separation', sep)
        corrected_contrast_curve[i] /= thrpt_interp
        interp_thrpts.append(thrpt_interp)

    if debug==True:
        #check throughput interpolations
        plt.plot(contrast_seps, interp_thrpts, 'k-', label="interpolated")
        plt.plot(thrpt_seps, thrpt_avgs, 'r--', label = "measured")
        plt.legend()
        plt.show()
        plt.clf()

    if savefig == True:
        plt.figure(figsize=(10, 5), dpi=750)
        plt.plot(contrast_seps * platescale, corrected_contrast_curve, label='corrected 5$\sigma$ contrast - average')
        cx = ['r--','g--','k--','c--','m--','r:','g:','k:','c:','m:']
        #plot the individual points
        for iter in np.arange(iterations):
            thrpts_thisset = np.interp(contrast_seps, thrpt_seps, thrpts[iter,:])
            plt.plot(contrast_seps * platescale, contrast/thrpts_thisset, cx[iter], label="set"+str(iter+1))
        plt.plot(contrast_seps * platescale, contrast, label='raw 5$\sigma$ contrast', color='gray')
        plt.yscale("log")
        plt.ylim(np.nanmin(contrast), 1e-1)
        plt.xlabel("distance in arcseconds")
        plt.ylabel("contrast")
        if IWA > 0:
            plt.plot((IWA * platescale, IWA * platescale), (1e-5, 1e-1), 'k--', label='IWA')
        # for bd in zone_boundaries*platescale:
        #   plt.plot((bd,bd),(0,1),'--',label='zone boundary')
        plt.legend()
        plt.savefig(outputdir + prefix + '_contrastcurve.jpg')
        plt.show()
        plt.clf()  # clear figure

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
        if len(df) > 1:
            colname = "ctrst_" + str(loc)
            # if column doesn't already exist, create and fill with nans
            if colname not in df:
                df[colname] = np.nan
                print("creating column", colname)
            df[colname].loc[(df.Dataset == subdir) & (df.pctcut == cut)] = contrast

    df.to_csv(dfname, index=False)
    contrast_out=np.array([contrast_seps,corrected_contrast_curve])
    fits.writeto(outputdir + prefix +  '_contrast.fits', contrast_out, overwrite=True)
    os.chdir(origdir)
    
    return (contrast_seps, corrected_contrast_curve, df, OWA)

def cut_comparison(subdir, wl, pctcuts=[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90], record_seps=[0.1, 0.25, 0.5, 0.75, 1.0],
                   contrast=1e-2, numann=3, movm=4, KLlist=[10], IWA=0, savefig=False, ghost=False, dfname="cuts.csv", debug=False,
                   iterations=3):
    """
    PURPOSE
    loop through data quality cuts and compile corrected contrast curves into single array

    REQUIRED INPUTS
    subdir = directory holding images that you want to reduce
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
    """
    #read in or make data frame
    df = get_cuts_df(dfname)

    ##loop through data quality cuts
    for cut in pctcuts:

        ##add logic to find existing files if they've already been computed!
        prefix = make_prefix(subdir, wl, cut)
        if os.path.exists(subdir+'contrastcurves/'+ prefix + '_a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA) +  '_throughputs.fits'):
            print ('found existing throughput file', subdir+'contrastcurves/'+ prefix +  '_a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA) + '_throughputs.fits')
            thrpt = fits.getdata(subdir+'contrastcurves/'+ prefix + '_a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA) + '_throughputs.fits')
            print(thrpt.shape)
            thrpt_seps = thrpt[0,:]
            thrpt_list = thrpt[1,:]

        else:
            print('computing throughputs for', cut, 'pct cut')
            thrpt_out, zone_boundaries, df, dataset_prefix = compute_thrpt(subdir, wl, cut,
                                                                    savefig=savefig, ghost=ghost, contrast=contrast,
                                                                    record_seps=record_seps,
                                                                    #KLIP parameters
                                                                    numann=numann, movm=movm,
                                                                    KLlist=KLlist, IWA=IWA, dfname=dfname, 
                                                                    debug=debug, iterations=iterations)

        if os.path.exists(subdir+'contrastcurves/'+ prefix +  '_a' + str(numann) + 'm' + str(movm) + '_contrast.fits'):
            print ('found existing contrast curve', subdir+'contrastcurves/'+ prefix +  '_a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA) + '_contrast.fits')
            curve = fits.getdata(subdir+'contrastcurves/'+ prefix +  '_a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA) + '_contrast.fits')
            contrast_seps = curve[0,:]
            corrected_curve = curve[1,:]

        else:
            print('computing contrasts for', cut, 'pct cut')
            contrast_seps, corrected_curve, df, OWA = make_contrast_curve(subdir, wl, cut,
                                                                 thrpt_out, dataset_prefix, record_seps=record_seps,
                                                                 savefig=savefig,
                                                                 ##KLIP parameters
                                                                 numann=numann, movm=movm, KLlist=KLlist, IWA=IWA,
                                                                 dfname=dfname, debug=debug)

        # compile contrasts for all cuts into array
        if cut == 0:
            # create a blank array to store contrast curves
            contrasts = np.zeros([len(pctcuts), len(corrected_curve)])
            contrasts[0, :] = corrected_curve
            i = 0
        else:
            contrasts[i, :] = corrected_curve

        # loop counter
        i += 1

    df.to_csv(subdir+dfname, index=False)

    return (contrast_seps, contrasts, zone_boundaries, IWA, df, OWA, dataset_prefix)


def contrastcut_fig(subdir, wl, contrast_seps, contrasts, zone_boundaries, dataset_prefix, OWA=225, IWA=0,
                    pctcuts=[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]):
    """
    PURPOSE
    Generate figure for comparison of contrasts with different image quality cuts

    REQUIRED INPUT
    subdir = directory holding images that you want to reduce
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
    outputdir = subdir + 'contrastcurves/'
    platescale = 0.0078513

    if subdir != './':
        prefix = subdir.replace('/', '_') + wl
    else:
        prefix = wl

    if OWA*platescale<1:
        OWA_asec = OWA*platescale
    else:
        OWA_asec=1

    # plt.style.use('seaborn-colorblind')
    plt.figure(figsize=(10, 5), dpi=750)
    for i in np.arange(len(pctcuts)):
        cut = pctcuts[i]
        j = i / len(pctcuts)
        if i % 2 == 0:
            linesty = '-'
        else:
            linesty = '--'
        plt.plot(contrast_seps * platescale, contrasts[i, :],
                 label=str(cut) + ' cut', linestyle=linesty, color=cm.plasma(j))
    floor = np.log10(np.nanmin(contrasts)) - 0.2
    plt.yscale("log")
    plt.title(prefix)
    plt.xlim(0, OWA_asec)
    plt.ylim(10 ** floor, 1e-1)
    plt.xlabel("distance in arcseconds")
    plt.ylabel("contrast")
    if IWA > 0:
        plt.plot((IWA * platescale, IWA * platescale), (1e-5, 1e-1), 'k--', label='IWA')
    for bd in zone_boundaries * platescale:
        if bd < 1:
            plt.plot((bd, bd), (0, 1), 'k-', lw=2, label='zone boundary')
    plt.legend()
    # write out in data directory
    plt.savefig(outputdir + dataset_prefix + '_contrastsbycut.jpg')
    plt.show()
    plt.clf()
    return


def clean_cuts(wl, subdir, keeplist, pctcuts=[0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]):
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


def inject_fakes(subdir, cut, IWA, wl='Line', numann=3, movm=5, KLlist=[10],
                 contrast=1e-2, seps=[10, 10, 10], thetas=[0, 120, 240], debug=True,
                 ghost=False, mask=[3, 15], slicefakes=True):
    """
    PURPOSE
    Injects false planets at specified locations, runs through KLIP, calculates throughput and
    displays a test image. This is largely the same as compute_thrpt, but for specific injected planet locations

    REQUIRED INPUTS
    subdir = directory holding images that you want to reduce
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

    OUTPUT
    separations where planets were injected = separations in pixels where throughputs were computed
    thrpt_list = computed throughputs at those separations (injected/recovered fake planet brightness)

    written by Kate Follette June 2019
    """

    # set up directories and naming
    outputdir = subdir+'fakes/'
    # if contrast curve directory doesn't already exist, create it
    if os.path.exists(outputdir) == False:
        os.mkdir(outputdir)
    prefix=subdir.replace('/', '_') + wl + '_' + str(cut) + 'pctcut_'
    strklip = 'a' + str(numann) + 'm' + str(movm) + 'iwa' + str(IWA) + 'KL' + str(KLlist[0])
    
    sepstr=''
    thetastr = ''
    for sep in seps:
        sepstr+=str(sep)+'_'
    for theta in thetas:
        thetastr+=str(theta)+'_'
    prefix_fakes = prefix +'thetas_'+thetastr+'seps_'+sepstr+'ctrst'+str(contrast)+'_FAKES'

    ###load data
    filelist = glob.glob(subdir + wl + '_' + str(cut) + 'pctcut_sliced' + "/sliced*.fits")
    # sorts file list so it's sequential
    filelist.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    ### pull the values of the star peak from the headers
    starpeak = []
    for i in np.arange(len(filelist)):
        head = fits.getheader(filelist[i])
        starpeak.append(head["STARPEAK"])
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

    if ghost == False:
        # where can, use actual image as "stamp" to avoid edge effects
        stmpsz = imsz
        starstamp = dataset.input[:, int(ycen - stmpsz / 2):int(ycen + stmpsz / 2),
                    int(xcen - stmpsz / 2):int(xcen + stmpsz / 2)]

        # multiply starstamp image by contrast to create your false planets
        to_inject = contrast * starstamp

    # now inject planets
    i = 0
    for sep in seps:
        if ghost == True:
            fakes.inject_planet(dataset.input, dataset.centers, np.repeat(contrast, nims)[0], dataset.wcs, sep,
                                thetas[i], fwhm=fwhm)
        else:
            fakes.inject_planet(dataset.input, dataset.centers, to_inject, dataset.wcs, sep, thetas[i], fwhm=fwhm,
                                stampsize=imsz)  # , thetas=thetas)
        i += 1

    # put NaNs back
    dataset.input[np.isnan(dataset_copy) == True] = np.nan

    #write out fakes cube
    fits.writeto(subdir + prefix_fakes + '.fits', dataset.input, overwrite=True)

    #slice the final fake dataset in preparation for parameter exploration
    if slicefakes==True:
        outdir=subdir+wl+'_'+str(cut)+'pctcut_FAKES_sliced/'
        if os.path.exists(outdir)== False:
            os.mkdir(outdir)
        for i in np.arange(0,len(dataset.prihdrs)):
            fits.writeto(outdir+'sliced_'+str(i)+'.fits', dataset.input[i,:,:], dataset.prihdrs[i], overwrite=True)

    # KLIP dataset with fake planets. Highpass filter here.
    parallelized.klip_dataset(dataset, outputdir=outputdir, fileprefix=prefix_fakes+strklip, algo='klip', annuli=numann,
                              subsections=1, movement=movm, numbasis=KLlist, calibrate_flux=False, mode="ADI",
                              highpass=True, save_aligned=False, time_collapse='median')

    # read in the KLIP cube that was just created
    klcube = fits.getdata("{out}/{pre}-KLmodes-all.fits".format(out=outputdir, pre=prefix_fakes+strklip))
    # compile specs for the injected planets
    planetspecs = (seps, thetas, mask)
    # make a SNRMap cube
    origdir = os.getcwd()
    os.chdir(outputdir)
    Output, snrs, snr_sums, snr_spurious, maskedims = snr.create_map(klcube, fwhm, saveOutput=True, outputName=prefix_fakes+strklip + '_SNRMap.fits',
                                       planets=planetspecs, checkmask=True, method='med')
    os.chdir(origdir)
    # create a list that will store throughputs
    thrpt_list = []

    # loop through planet locations and retrieve flux
    i = 0
    nimages = dataset.output.shape[1]
    fake_fluxes = np.zeros((n_planets, nimages))

    for sep in seps:
        theta = thetas[i]
        # thetas for retrieve planet call are measured from +x axis, so should always add 90 to pa for thetas keyword
        # dataset shape is [KL modes, n images, 1(wl dim), xdim, ydim]
        fake_flux = fakes.retrieve_planet_flux(dataset.output[0, :, 0, :, :], dataset.centers, dataset.wcs, sep, theta,
                                               searchrad=8, guessfwhm=fwhm,
                                               guesspeak=contrast * np.median(dataset.star_flux),
                                               thetas=np.repeat(theta + 90, dataset.output.shape[1]), refinefit=True)
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
    print("max snr for all planets is", snrs)

    return (dataset.input, dataset.prihdrs)


def collapsekl(subdir, fname, kllist, snrmeth='absmed', writestr='test'):
    """
    Reads in a parameter explorer file and collapses it in the KLmode dimension (axis 3)

    PARAMETERS:
    subdir: location of files
    fname: name of paramexplore file
    kllist: list of KL modes you want to collapse over
    snrmeth: one of two ways of computing SNR - possibilities are stdev and absmed (for median absolute value)

    RETURNS:
    avgkl:  array containing averages over the specified KL modes
    stdevkl: array containing standard deviations over the specified KL modes
    """
    # read in image and header
    klcube = fits.getdata(subdir + fname)
    head = fits.getheader(subdir + fname)

    #pull only the SNR map slice with the matching snrmeth value
    if snrmeth == "absmed":
        slice = 1
    if snrmeth == "stdev":
        slice = 0

    klcube=klcube[slice,:,:,:,:]

    dims = klcube.shape

    # set up a blank array to be filled
    klkeep = np.zeros([dims[2], dims[3], len(kllist)])

    # pull KL modes of parameter explorer from the header
    allkl = list(map(int, head['KLMODES'][1:-1].split(",")))

    # find only the KL slices in the list of modes to be collapsed and make an array
    i = 0;
    j = 0;
    keepind = []
    for kl in allkl:
        if kl in kllist:
            keepind.append(i)
            klkeep[:, :, j] = klcube[0, i, :, :]
            j += 1
        i += 1

    # replacing -250 (value when planet too close to annulus) with 0 for ease of display
    klkeep[klkeep == -1000] = np.nan

    # make mean and stdev arrays
    avgkl = np.mean(klkeep, axis=2)
    stdevkl = np.std(klkeep, axis=2)

    # update header to reflect KL modes used
    head["KLMODES"] = str(kllist)

    # write arrays
    fits.writeto(subdir+writestr + '_avgkl.fits', avgkl, head, overwrite=True)
    fits.writeto(subdir+writestr + '_stdevkl.fits', stdevkl, head, overwrite=True)
    return (avgkl, stdevkl)


def find_best(subdir, fname, kllist, snrmeth='absmed', writestr='test', weights=[1,1,0.5,0.5]):
    """
    collapses parameter explorer file and extracts the optimal parameter value

    PARAMETERS:
    subdir: location of files
    fname: name of paramexplore file
    kllist: list of KL modes you want to collapse over
    weights: weights for parameter quality metric. [SNR, SNR neighbors, stdev, stdev neighbors]

    RETURNS:


    """
    #collapse KL mode cube
    avgkl, stdevkl = collapsekl(subdir, fname, kllist, snrmeth=snrmeth, writestr=writestr)
    # stdevkl[stdevkl < 0.01] = np.nan

    #finds location of peak
    maxind = np.where(avgkl == np.nanmax(avgkl))

    #translates to x and y coordinates
    xcoord = maxind[0][0]
    ycoord = maxind[1][0]
    print("peak is at coordinates:", xcoord, ycoord)

    #extracts parameter explorer inputs from the filename (assumes was autonamed by parameter explorer)
    namelist = fname.split('_')
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

    #normalize SNR map
    snr_norm = avgkl / np.nanmax(avgkl)
    # stdevkl[avgkl<3]=np.nan
    # stdev_norm = 1/(stdevkl/np.nanmin(stdevkl))
    stdev_norm = stdevkl / avgkl
    stdev_norm /= np.nanmax(stdev_norm)
    stdev_norm = 1 - stdev_norm

    #computes neighbor quality by smoothing with Gaussian
    kern = conv.Gaussian2DKernel(stddev=2)
    nq_snr = conv.convolve(snr_norm, kern, preserve_nan=True)
    nq_stdev = conv.convolve(stdev_norm, kern, preserve_nan = True)

    #normalizes neighbor quality
    nq_snr /= np.nanmax(nq_snr)
    nq_stdev /= np.nanmax(nq_stdev)

    #calculate parameter quality metric by summing the four quantities
    agg = weights[0] * snr_norm + weights[1] * nq_snr + weights[2] * stdev_norm + weights[3] * nq_stdev

    ##find location or peak of parameter quality metric and print info
    ind = np.where(agg == np.nanmax(agg))

    #extract metric scores for this location
    metric_scores = [snr_norm[ind][0], nq_snr[ind][0], stdev_norm[ind][0], nq_stdev[ind][0], agg[ind][0]]

    #translate to annuli and mocement values
    ann_val = ymin + ind[0] * ystep
    movm_val = xmin + ind[1] * xstep

    print(
    'peak is at', [ind[0][0], ind[1][0]], 'corresponding to annuli', ann_val, ' and movement', movm_val, 'and snr of',
    avgkl[ind])
    print('metric scores for peak (snr, snr neigbors, stdev, stdev neighbors, agg) are:', metric_scores)
    return (snr_norm, nq_snr, stdev_norm, nq_stdev, agg, ann_val, movm_val, metric_scores)


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

def paramexplore_fig(subdir, fname, kllist, snrmeth='absmed', writestr='test', weights=[1,1,0.5,0.5], lowsnr=False):
    snr_norm, nq_snr, stdev_norm, nq_stdev, agg, ann_val, movm_val, metric_scores = \
        find_best(subdir, fname, kllist, snrmeth=snrmeth, writestr=writestr, weights=weights)
    namelist = fname.split('_')
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
    fig = plt.figure(tight_layout=True, figsize=(12,5))
    gs = fig.add_gridspec(2, 4)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    ax5 = fig.add_subplot(gs[:,2:])

    plt.setp((ax5), xticks=np.arange(nstepx + 1), xticklabels=np.arange(xmin, xmax + 1, xstep),
             yticks=np.arange(nstepy + 1), yticklabels=np.arange(ymin, ymax + 1, ystep))

    plt.setp((ax1, ax2, ax3, ax4), xticks=np.arange(nstepx + 1), yticks=np.arange(nstepy + 1), xticklabels=[], yticklabels=[])

    im1 = ax1.imshow(snr_norm, origin='lower', cmap='magma', vmin=0, vmax=1)

    ax1.set_xlabel("movement parameter")
    ax1.set_ylabel("annuli parameter")
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im1, cax=cax, orientation='vertical', label="Average Signal-to-Noise ratio")

    im2 = ax2.imshow(nq_snr, origin='lower',
                     cmap='magma', vmin=0, vmax=1)
    ax2.set_xlabel("movement parameter")
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im2, cax=cax, orientation='vertical', label="SNR Neighbor Quality")

    im3 = ax3.imshow(stdev_norm, origin='lower',
                     cmap='magma', vmin=0, vmax=1)
    ax3.set_xlabel("movement parameter")
    ax3.set_ylabel("annuli parameter")
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im3, cax=cax, orientation='vertical', label="1 - $\sigma$/$\mu$ metric (normalized)")

    im4 = ax4.imshow(nq_stdev, origin='lower',
                     cmap='magma', vmin=0, vmax=1)
    ax4.set_xlabel("movement parameter")
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im4, cax=cax, orientation='vertical', label="Stdev Neighbor Quality")

    # plot metric
    im5 = ax5.imshow(agg, origin='lower', vmin=0, vmax=np.nanmax(agg))
    ax5.set_ylabel("annuli parameter")
    ax5.set_xlabel("movement parameter")
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im5, cax=cax, orientation='vertical', label="Parameter Quality Metric")

    ind = np.where(agg == np.nanmax(agg))
    label_text = 'a' + str(ann_val[0]) + 'm' + str(movm_val[0])
    rect = patches.Rectangle((ind[1][0] - 0.5, ind[0][0] - 0.5), 1, 1, linewidth=2, edgecolor='r', facecolor='none')
    ax5.add_patch(rect)
    ax5.text(ind[1][0] + 0.75, ind[0][0], label_text, color='red')

    plt.savefig(writestr+'_paramqual.png')
    return(ann_val, movm_val)

def get_pe_df(dfname):
    # define dataframe if doesn't already exist
    try:
        df=pd.read_csv(dfname)
        print("found existing data frame") # with the following contents: \n", df)
        df_cols = df.columns.tolist()
    except:
        print("creating new data frame")
        #if none found, make one
        df_cols =['subdir','Dataset', 'Date', 'subset','pctcut', 'optimal movement', 'optimal annuli', 'IWA', 'KL combo used',
        'SNR value', 'SNR metric', 'SNR neighbors metric', 'stdev metric', 'stdev neighbors metric', 'total metric']
        # make a blank pandas dataframe
        df = pd.DataFrame({})
        for name in df_cols:
            df[name] = []
    return(df, df_cols)

def add_to_pe_df(subdir, pedir, fname, kllist, dfname='parameters.csv'):
    df, df_cols =get_pe_df(dfname)
    avgkl, stdevkl = collapsekl(subdir+pedir, fname, kllist)
    snr_norm, nq_snr, stdev_norm, nq_stdev, agg, ann_val, movm_val, metric_scores = find_best(subdir+pedir, fname, kllist)
    ind=np.where(agg == np.nanmax(agg))
    dir_info=re.split('/|_',subdir)
    objname = dir_info[0]
    date = dir_info[1]
    try:
        if (dir_info[2] in ['Line', 'Cont', 'Haplus']):
            subset = ''
        elif (dir_info[2]=='merged'):
            if dir_info[3] not in ['Line', 'Cont']:
                subset = dir_info[3]
            else:
                subset = ''
        else:
            subset = dir_info[2]
    except:
        subset = ''
    try:
        head = fits.getheader(pedir+fname)
        print(head)
        IWA = head["IWA"]
    except:
        IWA = int(input("You are using an old Parameter Explorer file. Please enter IWA:"))
    pedir_info = re.split('/|_', pedir)
    cut=[s for s in pedir_info if 'cut' in s]
    cut=(cut[0].split('pctcut'))[0]
    datalist=pd.DataFrame([[subdir,objname,date,subset,cut,movm_val[0],ann_val[0],IWA,str(kllist),avgkl[ind][0],snr_norm[ind][0], nq_snr[ind][0],
                            stdev_norm[ind][0], nq_stdev[ind][0], agg[ind][0]]], columns=df_cols)
    df = df.append(datalist)
    df.to_csv(dfname, index=False)
    return(df)

def clean_pe(subdir, keepparams):
    """
    Deletes superfluous parameter explorer that are not going to be used for analysis
    prefix: file prefix
    subdir: path to directory where files live
    keeplist: list of files to keep in [annuli, movem] form
    """
    fulllist = glob.glob(subdir+'*.fits')
    keeplist=[]
    for keep in keepparams:
        list = glob.glob(subdir+'*a'+str(keep[0])+'m'+str(keep[1])+'.*.fits')
        keeplist+=list
    pefiles = glob.glob(subdir+'paramexplore*.fits')
    keeplist+=pefiles
    for file in fulllist:
        if file not in keeplist:
            os.remove(file)

    return


def get_klip_inputs(subdir, pe_dfname='parameters.csv', cuts_dfname='cuts.csv'):
    """
    pulls info for naming and KLIP parameters from two data frames storing this info
    pe_dfname: dataframe containing info about optimal parameters
    subdir: directory containing data
    :return:
    """
    df, df_cols =get_pe_df(pe_dfname)
    df2=get_cuts_df(subdir+cuts_dfname)

    dir_info = re.split('/|_', subdir)
    objname = dir_info[0]
    date = dir_info[1]
    match = [direc for direc in df.subdir if subdir in direc]
    if len(match) > 1:
        print('there is more than one entry in the dataframe for this dataset. please delete the duplicates.')
        return
    else:
        cut = int(df[df["subdir"].values == match]["pctcut"].values[0])
        movm = int(df[df["subdir"].values == match]["optimal movement"].values[0])
        numann = int(df[df["subdir"].values == match]["optimal annuli"].values[0])
        fwhm = df2[df2["Dataset"] == subdir]["medfwhm"].values[0]
        IWA = int(df[df["subdir"] == subdir]["IWA"].values[0])
        kllist = df[df["subdir"] == subdir]["KL combo used"].values[0]
        kllist = eval(kllist)
    #if subdir == 'HD142527/8Apr14/merged_long/':
        #IWA = 6
    return (objname, date, cut, movm, numann, fwhm, IWA, kllist)


def klip_data(subdir, wl, params=False, fakes=False, imstring='_clip451_flat_reg_nocosmics_'):
    if params == False:
        objname, date, cut, movm, numann, fwhm, IWA, kllist = get_klip_inputs(subdir)
    else:
        [objname, date, cut, movm, numann, fwhm, IWA, kllist] = params
    if fakes==False:
        namestr = 'pctcut_sliced'
    else:
        namestr = 'pctcut_FAKES_sliced'
    #if hasn't already been sliced, slice it
    slicedir = subdir + wl + '_' + str(cut) + namestr+'/'
    if os.path.exists(slicedir) == False:
        print(wl, " image has not yet been sliced. Slicing now.")
        imname = subdir+wl+imstring+str(cut)+'pctcut.fits'
        rotoff_name = subdir+'rotoff_no'+wl+'cosmics_'+str(cut)+'pctcut.fits'
        SliceCube(imname, rotoff_name, slicedir=slicedir)
        pykh.addstarpeak(slicedir, debug=True, mask=True)

    prefix = subdir.replace('/', '_') + wl + '_' + str(cut) + 'pctcut_' + 'a' + str(numann) + 'm' + str(
        movm) + 'iwa' + str(IWA)
    if os.path.exists(subdir+prefix+'-KLmodes-all.fits'):
        klcube=fits.getdata(subdir+prefix+'-KLmodes-all.fits')
    else:   
        filelist = glob.glob(slicedir + "/sliced*.fits")
    
        dataset = MagAO.MagAOData(filelist, highpass=True) 
        parallelized.klip_dataset(dataset, outputdir=subdir, fileprefix=prefix, algo='klip', annuli=numann, subsections=1, movement=movm,
                              numbasis=kllist, calibrate_flux=False, mode="ADI", highpass=False, save_aligned=False, time_collapse='median')
    # fits.writeto('output_test.fits', dataset.output, overwrite=True)
        klcube = fits.getdata("{out}/{pre}-KLmodes-all.fits".format(out=subdir, pre=prefix))
    
    snmap, snrs, snr_sums, snr_spurious = snr.create_map(klcube, fwhm, saveOutput=True, outputName=subdir+prefix + '_SNRMap.fits')
    #snmap, snrs, snr_sums, snr_spurious = snr.create_map(klcube, fwhm, saveOutput=True, outputName=prefix + '_SNRMap.fits')
    return (klcube, snmap, fwhm)


def get_scale_factor(subdir):
    df3 = pd.read_csv('analysis/scale_factors/apphot_scales.csv')
    scale = df3[df3["subdir"] == subdir]["scale"].values[0]
    return (scale)


def run_redx(subdir, scale = False, imstring='_clip451_flat_reg_nocosmics_'):
    wls = ['Line', 'Cont']
    objname, date, cut, movm, numann, fwhm, IWA, kllist = get_klip_inputs(subdir)
    linecube, linesnr, linefwhm = klip_data(subdir, wls[0], imstring=imstring)
    contcube, contsnr, contfwhm = klip_data(subdir, wls[1], imstring=imstring)
    if scale == False:
        scale = get_scale_factor(subdir)
    else:
        scale = float(scale)
    sdicube = linecube - scale * contcube
    prefix = subdir.replace('/', '_') + 'SDI_' + str(cut) + 'pctcut_' + 'a' + str(numann) + 'm' + str(
        movm) + 'iwa' + str(IWA)
    sdisnr, snrs, snr_sums, snr_spurious= snr.create_map(sdicube, (linefwhm + contfwhm) / 2., saveOutput=True, outputName=prefix + '_SNRMap.fits')
    return (linecube, linesnr, contcube, contsnr, sdicube, sdisnr, prefix)


def indivobj_fig(subdir, lineim, contim, sdiim, prefix, stampsz=75, smooth=0, lims = False):
    """
    creates a three panel figure with line, continuum and SDI images

    optional keywords:
    lims = tuple in the form (min, max) to set colorbar limits

    """
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
    f.set_figwidth(15)
    f.set_figheight(10)
    pixscale = 0.0078513

    imsz = lineim.shape[1]

    ticks = (np.arange(11) - 10 / 2.) * 0.25 / pixscale
    ticklabels = str(ticks * pixscale)
    # ticks = np.arange(imsz)-(imsz-1)/2
    # ticklabels = str(ticks*pixscale)+'\"'

    cen = (imsz - 1) / 2
    low = int(cen - stampsz / 2)
    high = int(cen + stampsz / 2 + 1)

    if smooth != 0:
        gauss = conv.Gaussian2DKernel(stddev=smooth)
        lineimsm = conv.convolve(lineim, gauss, preserve_nan=True)
        contimsm = conv.convolve(contim, gauss, preserve_nan=True)
        sdiimsm = conv.convolve(sdiim, gauss, preserve_nan=True)

    # set up tick labels according to parameter ranges
    plt.setp((ax1, ax2, ax3), xticks=ticks, xticklabels=ticklabels,
             yticks=ticks, yticklabels=ticklabels)
    
    #user defined limits
    if lims != False:
        minm = lims[0]
        linemax = lims[1]
    #otherwise limits set by default
    else:
        linemax = np.nanstd(lineim[low:high, low:high])*5
        minm = -1 * linemax / 2
    if smooth > 0:
        im1 = ax1.imshow(lineimsm[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
    else:
        im1 = ax1.imshow(lineim[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
    # ax1.set_xlabel("")
    # ax1.set_ylabel("")
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im1, cax=cax, orientation='vertical', label="Intensity")

    if smooth > 0:
        im2 = ax2.imshow(contimsm[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
    else:
        im2 = ax2.imshow(contim[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
    # ax2.set_xlabel("movement parameter")
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im2, cax=cax, orientation='vertical', label="Intensity")
    # plot metric

    if smooth > 0:
        im3 = ax3.imshow(sdiimsm[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
    else:
        im3 = ax3.imshow(sdiim[low:high, low:high], vmin=minm, vmax=linemax, origin='lower', cmap='magma')
    # ax3.set_xlabel("movement parameter")
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.colorbar(im3, cax=cax, orientation='vertical', label="Intensity")
    plt.savefig(subdir+prefix+'.png')
    plt.show()
    plt.clf()
    return
