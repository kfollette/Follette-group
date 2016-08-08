import glob
from astropy.io import fits
import numpy as np
import math
import scipy
import time
import os
from scipy.ndimage.filters import gaussian_filter
import image_registration as imreg
import scipy.ndimage.interpolation as sci
import scipy.optimize as optimize
from visao_helpers import *

"""
NAME: visao_inventory

PURPOSE: Analyzes the images in a 'raw' folder in current directory and prints out number of science, dark, flat files,
        and the total rotation for the images

DESCRIPTION: Uses visao_getimtypes to obtain a dictionary of lists of all the needed keys words from the fits files
            in a raw folder in your current directory. Uses the data from 'AOLOOPST' to get number of science, dark,
            and flat frames. Looks at exposure times and if more than one exposure time is detected warns user.
            Detects what type of filter was used for the images (H-alpha, SII, OI, zp). If more than one filter is
            detected, warns user. Prints out total rotation for the images using requad_angle to deal with sets of
            images that decrease in angle or cross th 360/0 degree mark.

INPUTS: Optional lists of science, dark, and flat images, list of science image angles

OPTIONAL MAP KEYS: {'wfe': [], 'mag1': [], 'totrot': []}

OUTPUTS:
    return_imlists['sci_imlist'] : list of science images
    return_imlists['dark_imlist'] : list of dark images
    return_imlists['flat_imlist'] : list of flat images
    return_imLists['filt'] = filt
    return_imLists['mag1max'] = mag1max

HISTORY:
    CREATOR: Kate Follette?
    PY TRANS: 2016-07-07 by Wyatt Mullen, wmullen1@stanford.edu
"""
def visao_inventory(sci_imlist=None, dark_imlist=None, flat_imlist=None, rotoff_sciims=None, *keyWords, **keysMap):

    fnames = []
    imtypes = []
    imInfo = visao_getimtypes(fnames, imtypes, 'EXPTIME', 'VFW3POSN', 'VFW2POSN', 'AOLOOPST', 'ROTOFF', 'AVGWFE', 'AM',
                              subdir='raw', region='region')

    num_sciframes = num_darks = num_flats = 0
    dark_indices = []
    flat_indices = []

    for x in range(0, len(fnames)):
        if imtypes[x] == 0 and imInfo['AOLOOPST'][x] == 1:
            num_sciframes += 1
        if imtypes[x] == 2:
            dark_indices.append(x)
            num_darks += 1
        if imtypes[x] == 4:
            flat_indices.append(x)
            num_flats += 1

    if num_sciframes > 0:
        print('Found ' + str(num_sciframes) + ' closed loop science frames.')
    else:
        print('no science frames')

    if num_darks > 0:
        print('Found ' + str(num_darks) + ' darks')
    else:
        print('no dark frames')

    if num_flats > 0:
        print('Found ' + str(num_flats) + ' flats')
    else:
        print('no flat frames')

    exp_sort = sorted(imInfo['EXPTIME'])
    uniq_exp = list(set(exp_sort))
    if len(uniq_exp) > 1: print('Warning - more than one exposure time in this image set - separate before proceeding')

    if 'wfe' in keysMap:
        wfemax = keysMap['wfe']
    else:
        wfemax = 1000

    if 'mag1' in keysMap:
        mag1max = keysMap['mag1']
    else:
        mag1max = 1000

    # not sure why we are doing this loop for multiple exposures if we print out a warning beforehand
    exp_list = []
    exposure_type = {'alpha': [], 'S II': [], 'O I': [], 'z': []}
    for l in range(0, len(uniq_exp)):
        for image_num in range(0, len(fnames)):
            # select only closed loop images with same filter and exposure time, and with wfe cutoff if /WFE keyword is set
            for exp in exposure_type:
                if ((exp in imInfo['VFW3POSN'][image_num]) and (imtypes[image_num] == 0) and
                        (imInfo['AOLOOPST'][image_num] == 1) and (imInfo['EXPTIME'][image_num] == uniq_exp[0]) and
                        (imInfo['AVGWFE'][image_num] < wfemax)):
                    exposure_type[exp].append(image_num)

        if len(exposure_type['alpha']) > 0:
            print('Found ' + str(len(exposure_type['alpha'])) + ' closed loop h alpha with exposure time ' +
                  str(uniq_exp[l]) + ' and wfe< ' + str(wfemax))
            exp_list.append(0)

        if len(exposure_type['S II']) > 0:
            print('Found ' + str(len(exposure_type['S II'])) + ' closed loop [SII] with exposure time ' +
                  str(uniq_exp[l]))
            exp_list.append(1)

        if len(exposure_type['O I']) > 0:
            print('Found ' + str(len(exposure_type['S II'])) + ' closed loop [OI] with exposure time ' +
                  str(uniq_exp[l]))
            exp_list.append(2)

        if len(exposure_type['z']) > 0:
            print('Found ' + str(len(exposure_type['z'])) + ' closed loop zp with exposure time ' +
                  str(uniq_exp[l]))
            exp_list.append(3)

    types = ['alpha', 'S II', 'O I', 'z']
    nfilt = exp_list
    if len(nfilt) > 1:
        print('more than one SDI filter in this image set - separate before proceeding')
        sci_ims, filt = [], []
    elif len(nfilt) == 0:
        print('There are no science frames in this data set.')
        sci_ims, filt = [], []
    else:
        sci_ims = exposure_type[types[nfilt[0]]]
        types = ['Ha', 'SII', 'OI', 'zp']
        filt = types[nfilt[0]]
    #if len(nfilt) > 1: sci_ims = math.nan  # think this is correct, not completely sure yet

    if rotoff_sciims == None:
        rotoff_sciims = [imInfo['ROTOFF'][i] for i in sci_ims]
    if len(rotoff_sciims) != 0:
        rotoff_cont = requad_angles(rotoff_sciims)
        totrot = max(rotoff_cont) - min(rotoff_cont)
        print('Total rotation of this data set is ' + str(totrot) + ' degrees')

    return_imLists = {'sci_imlist': [], 'dark_imlist': [], 'flat_imlist': []}
    return_imLists['sci_imlist'] = [fnames[i] for i in sci_ims]
    print(return_imLists['sci_imlist'])
    return_imLists['dark_imlist'] = [fnames[i] for i in dark_indices]
    return_imLists['flat_imlist'] = [fnames[i] for i in flat_indices]
    return_imLists['filt'] = filt
    return_imLists['mag1max'] = mag1max

    # I think we can just put a stop when running the code by using debugger
    # if 'stp' in keyWords:
    return return_imLists

"""
NAME: visao_dark

PURPOSE: Create a master dark .fits file based on dark files in directory

DESCRIPTION: takes an inventory of all the images in a raw file and uses the lists returned. Looks at all the dark
            images in the dark_image list and takes the median of every pixel in that image cube. That median is then
            written to a new fits file in the current directory with a header that contains date created, exposure
            time, and the number of files from which the median was found.

OPTIONAL INPUTS:
    dark_imlist[]: list of all the dark images in a set, if not provided, obtains from visao_inventory
    master_dark[]: list of master_dark images to be combined into one "master_dark"

Outputs:
    master_dark.fits:   file written to current directory
    dark_imlist:        if entered in function call, returns list

HISTORY:
    CREATOR: 2016 by Kate Follette
    PY TRANS: 2016-07-11 by Wyatt Mullen, wmullen1@stanford.edu
"""
def visao_dark(dark_imlist = [], master_dark = []):
    imLists = visao_inventory()

    #check to make sure dark_imlist hasn't been provided or that there were no dark frames from visao_inventory
    if not dark_imlist and 'dark_imlist' in imLists:
        dark_imlist += imLists['dark_imlist']

    print('Creating master dark from ' + str(len(dark_imlist)) + ' dark frames.')
    dark_cube_list = []
    exp_time = []
    for dark in dark_imlist:
        #may need to include 'raw' here
        dark_fits = fits.open(dark)
        dark_data = dark_fits[0].data
        dark_header = dark_fits[0].header
        dark_cube_list.append(dark_data)
        exp_time.append(dark_header['EXPTIME'])

    darks = np.dstack(dark_cube_list)
    uniq_exp = list(set(exp_time))
    if len(uniq_exp) == 1:
        if len(dark_imlist) == 1:
            master_dark = darks
            fits.writeto('master_dark.fits', master_dark, clobber = True)
        else:
            master_dark = np.median(darks, axis=2) #think this will work
            fits.writeto('master_dark.fits', master_dark, clobber = True)
            f = fits.open('master_dark.fits', mode = 'update')
            dark_changes = f[0].header
            dark_changes['DATE'] = (time.strftime('%Y-%m-%d'), 'Creation UTC (CCCC-MM-DD) date of FITS header')
            dark_changes['EXPTIME'] = exp_time[0]
            dark_changes['NDARKS'] = len(dark_imlist)
            f.close()
    else:
        print('More than one exposure time in dark list - no dark created.')

"""
NAME: visao_separate_sdi

PURPOSE: Separates a series of images into the HA (line) emission and Continuum emission cubes

DISCRIPTION: Splits the given images vertically in half (in our case new dimensions are (512,1024)) so that they
        are split between the HA image and the Continuum image. The master dark frame created with visao_dark is
        subtracted from every image and then the image is divided by a flat frame if provided. The images are added
        to cubes or if indiv is specified are put into an aligned directory with increasing file names.
        A dictionary is returned that contains lists of values for each image of vwf3posn, ccd_gain, object, expt, as
        well as the list of line images and the list of continuum images.

OUTPUTS: Map of lists of values and images
        {'Line', 'Cont', 'vwf3posn', 'ccd_gain', 'object', 'expt'}

OPTIONAL KEYWORDS:  ['fits', 'indiv']

OPTIONAL MAP KEYS: {'Line': [], 'Cont': [], 'avgwfe'; 0, 'rotoff': 0, 'flat': 'file_name'}

HISTORY:
    CREATED: 2016 by Kate Follette
    PY TRANS: 2016-07-12 by Wyatt Mullen, wmullen1@stanford.edu
"""
def visao_separate_sdi(*keywords, **keysMap):
    imLists = visao_inventory()
    sci_imlist = imLists['sci_imlist']

    if not os.path.isfile('master_dark.fits'):
        print('You need to make a master dark first - please run visao_dark')
    else:
        master_dark = fits.open('master_dark.fits')[0]
    #makes a new aligned directory
    if not os.path.isdir('./aligned'):
        os.makedirs('aligned')
    #opens flat image if specified in function call
    if 'flat' in keysMap:
        flatim = fits.open(keysMap['flat'])[0]

    #dictionary filled with values that is returned from the function
    #this is instead of writing the arrays to fits files
    dataDict = {'expt': [], 'avgwfe': [], 'rotoff': [], 'object': [], 'vfw3posn': [], 'ccd_gain': []}
    #Empty image cubes which are filled and put into fits files
    Line, Cont = [], []
    i = 0
    for image in sci_imlist:
        i += 1
        temp_image = fits.open(image)[0]
        if 'flat' in keysMap:
            raw_im = (temp_image.data - master_dark.data) / flatim.data
        else:
            raw_im = (temp_image.data - master_dark.data)

        #splits images in half with one half going to the line emission cube, and one half to continuum emission
        dim1, dim2 = raw_im.shape
        if imLists['filt'] == 'Ha':
            line_im = raw_im[0:dim1/2, 0:dim2]
            cont_im = raw_im[dim1/2:dim1, 0:dim2]
        else:
            cont_im = raw_im[0:dim1/2, 0:dim2]
            line_im = raw_im[dim1/2:dim1, 0:dim2]
        Line.append(line_im)
        Cont.append(cont_im)

        temp_head = temp_image.header
        dataDict['expt'].append(temp_head['EXPTIME'])
        dataDict['avgwfe'].append(temp_head['AVGWFE'])
        dataDict['rotoff'].append(temp_head['ROTOFF'])
        dataDict['object'].append(temp_head['OBJECT'])
        dataDict['vfw3posn'].append(temp_head['VFW3POSN'])
        dataDict['ccd_gain'].append(temp_head['V47GAIN'])

        #creates two new files for each image (Ha and Cont) and places them in the aligned directory
        if 'indiv' in keywords:
            if 'flat' in keysMap:
                line_name = './aligned/Line_flat_' + str('%04d' % (i,)) + '.fits'
                cont_name = './aligned/Cont_flat_' + str('%04d' % (i,)) + '.fits'
                fits.writeto(line_name, Line[i-1], clobber = True)
                fits.writeto(cont_name, Cont[i-1], clobber = True)
            else:
                line_name = './aligned/Line_' + str('%04d' % (i,)) + '.fits'
                cont_name = './aligned/Cont_' + str('%04d' % (i,)) + '.fits'
                fits.writeto(line_name, Line[i-1], clobber=True)
                fits.writeto(cont_name, Cont[i-1], clobber=True)

    #checks to make sure exposure, objects, filters, and CCD sensitivity are constant
    expt = list(set(dataDict['expt']))
    object = list(set(dataDict['object']))
    vwf3 = list(set(dataDict['vfw3posn']))
    ccd_gain = list(set(dataDict['ccd_gain']))
    if len(expt) > 1: print('WARNING - More than one exposure time in this code!!!')
    if len(object) > 1: print('WARNING - More than one object in this cube!!!')
    if len(vwf3) > 1: print('WARNING - More than one SDI filter in this cube!!!')
    if len(ccd_gain) > 1: print('WARNING - CCD sensitivity is nt same for all images!!!')
    if master_dark.header['EXPTIME'] != expt[0]: print('WARNING - dark does not have same exposure time as images!!!')

    dataDict['Line'] = Line
    dataDict['Cont'] = Cont
    #writes image cube to files
    Line = np.array(Line)
    Cont = np.array(Cont)
    if 'fits' in keywords:
        if ('indiv' not in keywords) and ('flat' in keysMap):
            fits.writeto('Line_flat_preproc.fits', Line, clobber = True)
            fits.writeto('Cont_flat_preproc.fits', Cont, clobber = True)
        else:
            fits.writeto('Line_preproc.fits', Line, clobber = True)
            fits.writeto('Cont_preproc.fits', Cont, clobber = True)

    return dataDict

"""
NAME: visao_reg

PURPOSE: To perform a sub-pixel registration (centering) of both the continuum and line image cubes against a
        specified image index, and if desired to clip the images to nxn dimensions

DESCRIPTION: Uses visao_inventory to get a list of images and some general stats on the image set.
            Accepts an index of an image that it assigns as the reference image. Makes a 2d gaussian and registers the
            reference image against it to get an initially centered image. For every image in each data cube,
            'register_images()' returns an x,y shift to move the data cube image to the center based on a
            specified precision (entering 200 makes precision 1/200). Scipy.shift() then shifts the image by this
            amount and if specified, continuum images are scaled up by a specified scale factor (cont to HA is 1.02).
            Each image is then added to a new cube of dimensions nxn with the star centered and these cubes are written
            to new files.

INPUT:
    ref - number of reference image for which the program will compare all other images to (high quality PSF)

OPTIONAL KEYWORDS: [indiv]
        indiv:  if indiv used in previous steps, function will need to open file with indiv name

OPTIONAL MAP KEYS: {scl, fwhm, clip, sdi}
        scl:    scaling factor used to scale up cont (ex: 1.02)
        fwhm:   full-width half max to be used in creating the original centered gaussian
        clip:   the dimension of the output images in the cube (ex clip=400 would output 400x400 images)
        sdi:    the sdi (spectral differential imaging) scale factor (how much brighter is HA than cont)

OUTPUTS:
    Line/Cont Cubes:    Cubes of centered and clipped images

HISTORY:
    CREATED: 2014-06-17 by Kate Follette, kfollette@as.arizona.edu
        Modified April 2016 to include SDI images, write clip### into filename to avoid overwriting if want multiple FOV
        Modified 4/27/16 to register against a single high quality science image and split out separating channels, dark
            and flat to separate procedure visao_separate_sdi
    PY TRANS: 2016-07-12 by Wyatt Mullen, wmullen1@stanford.edu
"""
def visao_reg(ref, *keywords, **keysMap):

    imLists = visao_inventory()
    sci_imlist = imLists['sci_imlist']
    # variables I will use later on
    nims = len(sci_imlist)
    dummy_im = fits.open(sci_imlist[0])[0]
    dimY, dimX = dummy_im.shape
    xcen = (dimX - 1) / 2
    ycen = (dimY / 2 - 1) / 2

    if not 'indiv' in keywords:
        Line = fits.open('Line_flat_preproc.fits')[0]
        Cont = fits.open('Cont_flat_preproc.fits')[0]
        center_ref_line = Line.data[ref - 1, :, :]
        center_ref_cont = Cont.data[ref - 1, :, :]
        print('Registering against image: ' + sci_imlist[ref - 1])
    else:
        center_ref_line = fits.open('./indiv/Line_flat_' + str('%04d' % (ref,)) + '.fits')
        center_ref_cont = fits.open('./indiv/Cont_flat_' + str('%04d' % (ref,)) + '.fits')
        print('Registering against image: ' + './indiv/Line/Cont_flat_' + str('%04d' % (ref,)) + '.fits')

    if not 'FWHM' in keysMap:
        fwhm = 10
    else:
        fwhm = keysMap['fwhm']
    gauss_cen = make_Gaussian(dimX, dimY / 2, fwhm, (xcen, ycen))
    # smoothing to avoid cosmic rays
    center_ref_line_smooth = gaussian_filter(center_ref_line, sigma=5)
    center_ref_cont_smooth = gaussian_filter(center_ref_cont, sigma=5)

    # this is where we do most of the work
    # usfac is the precision so that if usfac=20, precision is 1/20 of pixel
    temp_arr1, temp_arr2 = np.zeros((Line.shape[1], Line.shape[2])), np.zeros((Line.shape[1], Line.shape[2]))
    dxl, dyl = imreg.register_images(center_ref_line_smooth, gauss_cen, usfac=1000)
    sci.shift(center_ref_line, (dyl, dxl), order=1, output=temp_arr1)
    center_ref_line = temp_arr1
    dxc, dyc = imreg.register_images(center_ref_cont_smooth, gauss_cen, usfac=1000)
    sci.shift(center_ref_cont, (dyc, dxc), order=1, output=temp_arr2)
    center_ref_cont = temp_arr2

    if 'clip' in keysMap:
        clip = keysMap['clip']
    else:
        clip = dimY
    # print(broccoli)
    # Initializing eventual return lists and some other helpful lists
    Lnim = Line.shape[0]
    Line_reg, Cont_reg, SDI_im = np.zeros((Lnim, clip, clip)), np.zeros((Lnim, clip, clip)),\
                                 np.zeros((Lnim, clip, clip))
    Line_smooth, Cont_smooth = np.zeros(Line.shape), np.zeros(Line.shape)
    # main for loop
    for i in range(0, nims):
        Line_smooth[i, :, :] = gaussian_filter(Line.data[i, :, :], sigma=5)
        Cont_smooth[i, :, :] = gaussian_filter(Cont.data[i, :, :], sigma=5)
        temp_arr3, temp_arr4 = np.zeros((Line.shape[1], Line.shape[2])), np.zeros((Line.shape[1], Line.shape[2]))
        shiftXL, shiftYL = imreg.register_images(Line_smooth[i, :, :], center_ref_line, usfac=200)
        sci.shift(Line.data[i, :, :], (shiftYL, shiftXL), order=1, output=temp_arr3)
        Line.data[i, :, :] = temp_arr3
        shiftXC, shiftYC = imreg.register_images(Cont_smooth[i, :, :], center_ref_cont, usfac=200)
        sci.shift(Cont.data[i, :, :], (shiftYC, shiftXC), order=1, output=temp_arr4)
        Cont.data[i, :, :] = temp_arr4

        #Zoom option doesn't work with any nan values, will have to find other option
        if 'scl' in keysMap:
            scl = keysMap['scl']
            temp_arr5 = np.zeros((scl*Line.shape[1], scl*Line.shape[2]))
            sci.zoom(Cont.data[i,:,:], scl, output = temp_arr5, order=1)
            Cont.data[i,:,:] = temp_arr5[ycen * scl - dimY / 4: ycen * scl + dimY / 4,
                                         xcen * scl - dimX / 2: xcen * scl + dimX / 2]

        #print('Line x,y = (' + str(shiftXL) + ',' + str(shiftYL) + ')     Cont x,y = ' + str(shiftXC) + ' ' + str(
            #shiftYC))
        # making sure clip doesn't overshoot frame in some images for line/cont
        if dimX / 2 - math.fabs(shiftXL) < clip / 2:
            print('Line overshoots by ' + str(dimX / 2 - math.fabs(shiftXL) - clip / 2) + ' pixels in x.')
            Line_cropx = math.floor(dimX / 2 - math.fabs(shiftXL))
        else:
            Line_cropx = clip / 2
        if dimY / 4 - math.fabs(shiftYL) < clip / 2:
            print('Line overshoots by ' + str(dimY / 4 - math.fabs(shiftYL) - clip / 2) + ' pixels in y.')
            Line_cropy = math.floor(dimY / 4 - math.fabs(shiftYL))
        else:
            Line_cropy = clip / 2
        if dimX / 2 - math.fabs(shiftXC) < clip / 2:
            print('Cont overshoots by ' + str(dimX / 2 - math.fabs(shiftXC) - clip / 2) + ' pixels in x.')
            Cont_cropx = math.floor(dimX / 2 - math.fabs(shiftXC))
        else:
            Cont_cropx = clip / 2
        if dimY / 4 - math.fabs(shiftYC) < clip / 2:
            print('Cont overshoots by ' + str(dimY / 4 - math.fabs(shiftYC) - clip / 2) + ' pixels in y.')
            Cont_cropy = math.floor(dimY / 4 - math.fabs(shiftYC))
        else:
            Cont_cropy = clip / 2

        # reassigning Line_reg and Cont_reg (yes I know the assignment is confusing, don't think about it)
        Line_reg[i, clip / 2 - Line_cropy:clip / 2 + Line_cropy, clip / 2 - Line_cropx:clip / 2 + Line_cropx] = \
            Line.data[i, ycen - Line_cropy + 1:ycen + Line_cropy + 1, xcen - Line_cropx + 1:xcen + Line_cropx + 1]
        Cont_reg[i, clip / 2 - Cont_cropy:clip / 2 + Cont_cropy, clip / 2 - Cont_cropx:clip / 2 + Cont_cropx] = \
            Cont.data[i, ycen - Cont_cropy + 2:ycen + Cont_cropy + 2, xcen - Cont_cropx + 2:xcen + Cont_cropx + 2]

        # shifting the array if the clip value is odd numbered
        if clip % 2 != 0:
            temp_arr6, temp_arr7 = np.zeros((Line.shape[1], Line.shape[2])), np.zeros((Line.shape[1], Line.shape[2]))
            sci.shift(Line_reg[:, :, i], 0.5, output=temp_arr6)
            Line_reg[:, :, i] = temp_arr6
            sci.shift(Cont_reg[:, :, i], 0.5, output=temp_arr7)
            Cont_reg[:, :, i] = temp_arr7

        print('Processed image ' + str(i + 1) + ' of ' + str(nims))

        if 'sdi' in keysMap:
            sdi = keysMap['sdi']
            SDI_im[:, :, i] = Line_reg[:, :, i] - sdi * Cont_reg[:, :, i]

    if 'clip' in keysMap:
        fits.writeto('Line_clip' + str('%03d' % (clip,)) + '_reg-2.fits', Line_reg, clobber=True)
        fits.writeto('Cont_clip' + str('%03d' % (clip,)) + '_reg-2.fits', Cont_reg, clobber=True)
    else:
        fits.writeto('Line_reg.fits', Line_reg, clobber=True)
        fits.writeto('Cont_reg.fits', Cont_reg, clobber=True)

    if 'sdi' in keysMap:
        if 'clip' in keysMap:
            fits.writeto('SDI_sc' + str('%.2f' % clip) + '_clip' + str('%03d' % (clip,)) + 'reg.fits',
                         SDI_im, clobber=True)
        else:
            fits.writeto('SDI_sc' + str('%.2f' % clip) + 'reg.fits', SDI_im, clobber=True)


"""
NAME: py_idl_compare

PURPOSE: Compares two cubes of images of the same shape (preferably IDL vs Python) to compare coding methods/outputs

DESCRIPTION: Takes in two image cubes of a specified type (for example Line, Cont, small data, etc) and saves a cube
            of the differences between the images by subtracting the second cube from the first. This cube can be used
            to find differences in algorithms or coding methods. Additionally a single median difference image is
            written to a file which is the median image of the cube.

Inputs:
    cube1, cube2:   two image cubes run through different algorithms or languages to compare
    type:           a string specifying the type of image cubes enered (Line, Cont, etc)

Outputs: returns nothing
    'Dif_cube_'type'.fits:      Cube of the differences between each image in the cube
    'Median_dif_'type'.fits:    Median image of the difference cube

HISTORY:
    Created: 2016-07-13 by Wyatt Mullen, wmullen1@stanford.edu
"""
def py_idl_compare(cube1, cube2, type):
    cube1 = fits.open(cube1)[0]
    cube2= fits.open(cube2)[0]

    if cube1.shape != cube2.shape:
        raise RuntimeError('The shape of these cubes is not the same!!!!')

    out_cube = np.zeros(cube1.shape)
    for i in range(0, cube1.shape[0]):
        out_cube[i,:,:] = cube1.data[i,:,:] - cube2.data[i,:,:]

    fits.writeto('Dif_cube_' + type + '.fits', out_cube, clobber=True)

    median_dif = np.median(out_cube, axis=0)
    fits.writeto('Median_dif_' + type + '.fits', median_dif, clobber=True)


"""
NAME: cube_median

PURPOSE/DESCRIPTION: Takes the median of a given cube and writes it to a fits image

INPUTS:
    image_cube: fits image_cube for which we are taking the median
    type:       string for type (Line, IDL_Cont) to be added to output file name

OUTPUTS:
    median image: median image of the image cube

HISTORY:
    CREATED: 2016-07-25 by Wyatt Mullen, wmullen1@stanford.edu
"""
def cube_median(image_cube, type):
    image_cube = fits.open(image_cube)[0]
    med = np.nanmedian(image_cube.data, axis=0)
    fits.writeto('im_cube_med_' + type + '.fits', med, clobber=True)
