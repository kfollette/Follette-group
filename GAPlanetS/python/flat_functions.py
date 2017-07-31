"""
A number of flat processing functions in order to create a single processed flat image from a number of raw fits files
as well as dust-masked images for spots with specific attenuation rates. Some of the functions are out of date or
aren't that useful such as box_blur and flat_analysis (meaning they don't function as well as they could). The meat
of the class?? function?? is the raw_flat_process, the dust_mask, and the batch_process.

Most of these functions are written for specified folder structures, but they could be generalized to work better.
Furthermore if someone could find a way to automate the entire process (including unzipping the fits files, splitting
the darks from the flats from the science) that would make the entire process much quicker. Currently the part that
requires the most human intervention is writing the files and their long, descriptive names. If these names could be
included in the headers (not difficult, but not something I've done yet) it might be easier.

-Wyatt Mullen (wmullen1@stanford.edu)

Contains: 11 functions for processing flat images
Includes:
    raw_flat_process()
    box_blur()
    flat_analysis()
    dust_mask()
    flat_stats()
    count_spots()
    flat_compare()
    dust_vals()
    batch_dark()
    batch_mask()
    batch_compare()
"""

import numpy as np
from astropy.io import fits
import bottleneck as bn
import math
import matplotlib.pyplot as plt
import glob
import visao_process as vis
from scipy.ndimage.filters import uniform_filter
from scipy.ndimage.filters import median_filter
from astropy.stats import sigma_clip
import os

"""
##################################################################################################################
NAME: raw_flat_process

PURPOSE: Takes a directory containing a subdirectory of raw flat files and combines them into a processed flat image
        with a mean pixel value approximately of one.

DESCRIPTION: Given a name string, a directory where the raw files are located, and possibly a dark, will collect all
            the files in that given directory. For each image the specified dark is subtracted from it,
           (due to a difference between sensitivity of bottom and top half of image) a constant is added to
           the bottom amplifier, a uniform filter of x (most likely 41 pixels) smooths a copy, the image is divided
           by this smooth image (highlights the dust spots), and then the image is added to an image cube. Once every
           image has been processed a median is taken of the cube and this cube is written to a fits file.

           Possibly the median filter is a better option (was used by Jared) but takes almost 100 times as long
           (dozens of minutes instead of a few seconds). Also I do not use a 3-sigma clipped mean but rather a median.

INPUTS:
    namestr:    a string added to the name of the written file which identifies the processed flat file
    directory:  a string, specified where the raw flat images are located
    dark:       a bias frame to subtract from each image, either specified or the code may be able to find one,
                    preferable to use a dark from image taken at the same time as the flats (same exp length)
    blur:       integer, area of box blur, a larger number increases the blurriness (dust spots more prominent)

OUTPUTS:
    'visao_flat_namestr_med41.fits':    processed flat file

HISTORY:
    Created: 2016-08-10 by Wyatt Mullen, wmullen1@stanford.edu
"""
def raw_flat_process(namestr, directory='raw', dark=None, blur=41):
    searchstr = directory + '/*.fits'
    files = glob.glob(searchstr)
    if dark == None:
        #imList = vis.visao_inventory()
        #if not imList['dark_imlist']:
            #raise FileNotFoundError('Please specify a master dark to use.')
        #else:
        vis.visao_dark(subdir=directory[:-3] + 'darks', name=namestr)
        dark = fits.open('master_dark.fits')[0]
    else:
        dark = fits.open(dark)[0]

    flat_cube = np.zeros((len(files), dark.shape[0], dark.shape[1]))
    im_num = 0
    for image in files:
        temp_im = fits.open(image)[0]
        #subtracting dark frames (bias)
        raw_im = temp_im.data - dark.data
        #adding constant to bottom amplifier in the image
        bot_med = np.nanmedian(raw_im[511:512, 0:1024])
        top_med = np.nanmedian(raw_im[512:513, 0:1024])
        raw_im[0:512, 0:1024] += top_med - bot_med
        #smoothing by median filter
        smooth_im = uniform_filter(raw_im, size=blur)
        #smooth_im = median_filter(raw_im, size=41)
        #dividing by smoothed image
        raw_im = raw_im / smooth_im
        #adding image to cube of images
        flat_cube[im_num,:,:] = raw_im
        im_num += 1
        print('Completed image ' + str(im_num) + ' of ' + str(len(files)))

    #taking median of image cube for each pixel
    master_flat = np.nanmedian(flat_cube, axis=0)
    fits.writeto('visao_flat_' + namestr + '_med' + str(blur) + '.fits', master_flat, clobber=True)

"""
##################################################################################################################
NAME: box_blur

PURPOSE/DESCRIPTION: An option for blurring images to be used in raw_flat process. Draws a box of size n around each
                    pixel and replaces the pixel at the center with the median value of the pixels in the box. Slower
                    and not as effective as uniform_filter.

INPUT:
    image:      2-dimensional array of values that represent an image
    size:       size of box to be created. For examplse size=20 would create a box 20 pixels on each side. Should pick
                    an odd number to make smoothing easier

OUTPUT:
    2-d smoothed array of the original image

HISTORY:
    Created: 2016-08-07 by Wyatt Mullen, wmullen1@stanford.edu
"""
def box_blur(image, size):
    image = np.array(image)
    temp_image = image
    for y in range(size//2, image.shape[0] - size//2 - 2):
        for x in range(size//2, image.shape[1] - size//2 - 2):
            temp_image[y][x] = np.nanmedian(image[-size//2 + y:size//2 + y, -size//2 + x:size//2 + x])
    return np.array(temp_image).tolist()


"""
##################################################################################################################
NAME: flat_analysis (depreciated)

PURPOSE: Interactive function for analyzing flat before and after being masked. Will create a masked processed flat.

DESCRIPTION: Asks for user to input an unmasked image and appends it to a list of image. The statistics of the flat
            before the mask are printed and then the code asks for the minimum value for the mask, the date string, and
            whether you want to make a histogram of the data. This repeats until you respond you have no other images
            you want to process, then function compares images in the list and prints out the number of differences
            in the location of spots for the images. This functionality does not work very well (or at all), so don't
            pay much attention to it.

INPUTS: asks for: image name, min pixel value, date string, and whether you want a histogram of pixel values

OUTPUTS: prints statistics, plots pixel histogram if desired, returns the location of the spots

HISTORY:
    Created: 2016-08-05 by Wyatt Mullen, wmullen1@stanford.edu
"""
def flat_analysis():
    spotsDict = {}
    images = []
    while True:
        im = input('Name of an unmasked image: ')
        image = fits.open(im)[0]
        images.append(im)
        print('These are the statistics before the mask is applied.')
        spots = flat_stats(image)
        minval = input('Minimum pixel value unmasked in image: ')
        date = input('Date String: ')
        plot = input('Would you like to make a histogram of the data? (y/n): ')
        if len(spots) == 0:
            mask_image = dust_mask(image, float(minval))
            fits.writeto('visao_flat_' + str(date) +'_nansmask.fits', mask_image, clobber=True)
            print('These are the statistics after the mask is applied.')
            if plot == 'y':
                spotsDict[im] = flat_stats(mask_image, True)
            else:
                spotsDict[im] = flat_stats(mask_image, False)
        else:
            print('This image has already been masked.')
        done = input('Do you have any more images? (y/n): ')
        if done == 'n':
            break
    for i in range(0,len(images)):
        if len(spotsDict) > 1 and len(spotsDict) > i + 1:
            im1_spots = spotsDict[images[i]]
            print(im1_spots)
            im2_spots = spotsDict[images[i+1]]
            print(im2_spots)
            print()
            print(images[i+1] + ' has ' + str(len(im2_spots) - len(im1_spots)) + ' more spots than ' + images[i])
            print('Spots in ' + images[i+1] + ' but not in ' + images[i])
            im2_spots = [tuple(np.array(x)//10) for x in im2_spots]
            im1_spots = [tuple(np.array(x)//10) for x in im1_spots]
            numdif = 0
            for i in im2_spots:
                if i not in im1_spots:
                    print(tuple(np.array(i) * 10))
                    numdif += 1
            print('There are ' + str(numdif) + ' spots dif')
            numdif = 0
            for i in im1_spots:
                if i not in im2_spots:
                    print(tuple(np.array(i) * 10))
                    numdif += 1
            print('There are ' + str(numdif) + ' spots dif')

    return spots


"""
##################################################################################################################
NAME: dust_mask

PURPOSE: Create a mask around dust spots for flat images where the dust spots are masked as NaNs

DESCRIPTION: The function takes an image, a minimum opacity (out of 1.00), and a radius in pixels. The bad values are
            found first and replaced with nans. Then the program loops through the rest of the image and replaces any
            pixels within a circular radius 'rad' with the value NaN. The mask is then returned.

INPUTS:
    image:  flat image to be masked
    minval: minimum value that will ot be masked in the flat image. Default is 0.98
    rad:    radius of pixels around center mask to make sure that dust spots will have minimal impact. Default is 2
    sem:    string, if specified writes a fits file of the mask
    ones:   boolean, if True makes mask ones and NaNs

OUTPUTS:
    mask:   2-d masked image array

HISTORY:
    CREATED: 2016-07-20 by Wyatt Mullen, wmullen1@stanford.edu
    UPDATED: 2016-08-16 by WM -- added writing to file functionality and option to make a ones/NaN mask
"""
def dust_mask(flat_namestr, minval=0.98, rad=2, ones=False):
    #masking all values below the minval
    image = fits.open(flat_namestr+'.fits')[0] #added
    nans = np.where(image.data < minval)
    data = np.where(image.data >= minval)
    dimy = image.shape[0]
    dimx = image.shape[1]
    mask = np.zeros(image.shape)
    temp_mask = np.zeros(image.shape)
    mask[nans] = np.NaN
    temp_mask[nans] = np.NaN
    if ones:
        mask[data] = 1
    else:
        mask[data] = image.data[data]

    #putting pixel radius around every dust spot
    #temp_mask = mask
    for y in range(0,dimy):
        for x in range(0,dimx):
            if np.isnan(temp_mask[y][x]):
                for i in range(-rad, rad+1):
                    for j in range(-rad, rad+1):
                        if 0 <= x + j and x+j < dimx and 0 <= y+i and y+i < dimy:
                            mask[y+i][x+j] = np.NaN

#if sem != '': #create masked file if a semester is specified
        fits.writeto(flat_namestr + '_dustmask_' + str(minval) + '_' + str(rad) + 'pix.fits', mask, clobber=True) #added
    return mask

"""
##################################################################################################################
NAME: flat_stats

PURPOSE: Prints statistics for a single already opened flat image

DESCRIPTION: Prints dimensions, median, mean, std, #nan, and #dust spots to the console. Calculates the number of
            dust spots within roughly 20 pixels of the edge by using recursive backtracking where the loop in this
            functions calls count_spots() for every pixel in the image. Also prints the locations of the lower left
            corner of each dust spot. If plot is True then will create a histgram of all the pixel values in the flat
            mask.

INPUTS:
    image:      first element of hdulist (already opened image) so that it is passed the fits.open()[0]
    plot:       boolean, if True then a popup histgram will emerge of all the pixel values. Should have a definite
                    normal shape.

OUTPUTS:
    dust_spots:     list of integer tuples which are the locations of the dust_spots
    plot:           histogram (see above)

HISTORY:
    Created: 2016-08-05 by Wyatt Mullen, wmullen1@stanford.edu
"""
def flat_stats(image, plot=False):
    print('Dimensions (y,x): ' + str(image.shape))
    print('Median: ' + str(bn.nanmedian(image.data)))
    print('Mean: ' + str(bn.nanmean(image.data)))
    print('Standard Deviation: ' + str(bn.nanstd(image.data)))
    print()

    Nans = np.where(np.isnan(image.data))
    print('Number of NaN: ' + str(Nans[0].size))
    print('Calculating number of dust spots...')
    temp_im = np.array(image.data)[20:1000, 20:1000]
    num_spots = 0
    dust_spots = []
    for y in range(0, temp_im.shape[0]):
        for x in range(0, temp_im.shape[1]):
            if count_spots(temp_im, temp_im[y][x], y, x):
                num_spots += 1
                dust_spots.append((y+20, x+20))
    if num_spots == 0:
        print('This image has no dust spots.')
    else:
        print('Number of dust spots: ' + str(num_spots))
        print('Dust spot lower right corner coordinates (y,x):')
        print(dust_spots)
    #plots a histogram of the pixel counts in the image
    if plot:
        im_list = list(image.data.flatten())
        not_nans = np.where(~np.isnan(im_list))[0].tolist()
        clean_im_list = [im_list[i] for i in not_nans]
        plt.hist(clean_im_list, bins=1000, range=(0.98,1.02), color='blue', histtype='stepfilled')
        plt.show()
    return dust_spots

"""
##################################################################################################################
NAME: count_spots

PURPOSE/DESCRIPTION: recursive function that goes every spot in a masked image and counts the spots. Cannot figure
                    out how to deal with recursion overflow error so if anybody has any ideas I would appreaciate it.
                    If the pixel is not a nan-pixel then the function returns false, otherwise it continues to
                    recurse until it runs out of nans after a nan is found it is replaced with a default value of 100
                    in a temporary image.

INPUTS:
    temp_im:    2-d array full of floats, allows us to replace nans with 100 without changing the image
    pixel:      function is called for each pixel in the image, have to keep track of the value of the pixel
    y,x:        coordinates of the pixel

OUTPUTS:
    returns False if pixel is not a nan and True if pixel is a nan, if one nan is found will always eventually
    return True

HISTORY:
    Created: 2016-08-01 by Wyatt Mullen, wmullen1@stanford.edu
"""
def count_spots(temp_im, pixel, y, x):
    #spot = False
    if not np.isnan(pixel) or y < 0 or x < 0 or y + 2 > temp_im.shape[0] or x + 2 > temp_im.shape[1]:
        return False
    else:
        temp_im[y][x] = 100.0
        for i in range(-1,2):
            for j in range(-1,2):
                count_spots(temp_im, temp_im[y+i][x+j], y + i, x + j)
            #if RecursionError: return True
        return True


"""
##################################################################################################################
NAME: flat_compare

PURPOSE: Compares two flat images and creates a new image where different values highlight differences in the
        in the dust spots between years or semesters.

DESCRIPTION: Takes two images as input (usually designed for two masked flats from the same semester) and compares the
            locations of the nans values where the 'dust' spots have been masked. The two images are opened and where
            nans only show up in the first image become ones, where nans show up in both images become 2s, and where
            nans show up in only the second image become 3s. Every other pixel is a zero and this image is then written
            out.

INPUTS:
    im1/im2:    strings, two masked flat images
    year:       string, written into the header, identifies what year or period the flats are from

OUTPUTS:
    'flat_compare_*year*_.fits':    comparison images made of 0s, 1s, 2s, and 3s

HISTORY:
    Created: 2016-07-23 by Wyatt Mullen, wmullen1@stanford.edu
"""
def flat_compare(im1, im2, im3, year):
    im1 = fits.open(im1)[0]
    im2 = fits.open(im2)[0]
    im3 = fits.open(im3)[0]
    im_out = np.zeros(im1.shape)
    ones = np.where((np.isnan(im1.data))) # & (~np.isnan(im2.data))
    twos = np.where((~np.isnan(im1.data)) & (np.isnan(im2.data)))
    threes = np.where((~np.isnan(im2.data)) & (~np.isnan(im1.data)) & (np.isnan(im3.data)))
    im_out[ones] = 1
    im_out[twos] = 2
    im_out[threes] = 3

    fits.writeto('flat_compare_' + str(year) + '.fits', im_out, clobber=True)

"""
##################################################################################################################
NAME: dust_vals

PURPOSE: Simply prints the total number of pixels that are less than or greater certain values

DESCRIPTION: Helps print flat statistics on an image by taking the values specified in the function. For values
            greater than 1, prints the number of pixels greater than that value. For values less than 1, prints
            the number of pixels less than that value.

INPUTS:
    im:     string, name of flat image to be analyzed

OUTPUTS:    prints # of pix fitting certain parameters

HISTORY:
    Created:    2016-07-20 by Wyatt Mullen, wmullen1@stanford.edu
"""
def dust_vals(im):
    im = fits.open(im)[0]
    pixVals = [1.06, 1.04, 1.02, 0.98, 0.96, 0.94, 0.92, 0.9, 0.85, 0.8, 0.75, 0.7]
    pixList = []
    for val in pixVals:
        if val > 1.00:
            a = np.where(im.data > val)
            #print('NumPix greater than ' + str(val) + ': ' + str(len(a[0])))
            pixList.append(len(a[0]))
        else:
            a = np.where(im.data < val)
            #print('NumPix less than ' + str(val) + ': ' + str(len(a[0])))
            pixList.append(len(a[0]))

    return pixList

"""
##################################################################################################################
NAME: batch_dark

PURPOSE/DESCRIPTION: goes through all the visao_* flat folders in a certain semester and for each folder, goes into
                    all of its subfolders that contain a 'sec' exposure time in the heading and runs visao_dark. If no
                    darks exist in the folder does not create a dark, otherwise creates a dark with the name of the
                    folder and the directory inside the directory. Somewhat specialized use, but allows many folders to
                    be checked quickly to see whether they contain darks or flats or science images.

INPUTS:     None

OUTPUTS:    From visao_dark (tells you whether a dark has been created or not in each directory)

HISTORY:
    Created:    2016-08-15 by Wyatt Mullen, wmullen1@stanford.edu
"""
def batch_dark():
    dirs = glob.glob('visao_*')
    for d in dirs:
        os.chdir(d)
        sub = glob.glob('*sec*')
        for s in sub:
            string = s + '/raw'
            vis.visao_dark(subdir=string, name=d[6:14] + '_' + s)
        os.chdir('..')

"""
##################################################################################################################
NAME: batch_mask

PURPOSE: Quickly masks a number of processed flat files in a directory and writes a file with statistics on number
        of pixels by value and number of spots.

DESCRIPTION: Takes a list of names and a min value for which to mask the flats at. For every file in the current
            directory with a name of 'visao_flat_2*', the program collects the files and begins to process them.

            IMPORTANT: Run batch mask first without the list of file names so that you can determine the order
            of files that the program will use from the print statement. After the program crashes the first time,
            you can use that knowledge to put your list of names in the correct order. This was the best way I could
            think of having some control over the names instead of flat_1, flat_2, ...

            ALSO: Make sure that in your file you have a directory named "minval"_level, otherwise the program will
            give you an error. For example if you are masking at the 0.95 level make a directory named '0.95_level'.

            Similar to the implementation of flat_stats, this function writes the number of pixels above values greater
            than one and number of pix below values less than one, but it writes it to a
            file named 'flat_stats_min.txt'. These statistics include the median, mean, and std of the unmasked flat.
            recursive backtracking is used with count_spots to count the number of spots and this is also added to
            the text file (this will be the main difference between txt files of the same semester but different min
            values).

            WARNING: The count spots method does not work if there are too many spots. The program will reach it's
            maximum recursion depth and will crash. If this happens, don't worry. Simply remove the offending image
            from the current directory (place it in a directory called bad images or something like that) and make sure
            to also remove the name from the list of names you entered as parameters. This may happen multiple times.
            The masked image will still be made, but the final flat_stats_min.txt will not include info about it.

INPUTS:
    names:  list of strings, see the IMPORTANT above, however make sure string list and order of files matches up
            (it is not the same order of files you see in your terminal/finder)
    min:    float, minimum value to mask, fewer spots will be masked at smaller values, standard values: 0.96, 0.98

OUPUTS:     does not return anything
    'flat_stats_min.txt':                       statistics for all the images at different specified min values
    'min_level/visao_flat_mask_name_min.fits':  masked files, however many files are in the directory

HISTORY:
    Created: 2016-08-17 by Wyatt Mullen, wmullen1@stanford.edu

BUGS: Would like to make it easier to name files and not have my spot counting algorithm crash when it reaches max
        recursion depth. Instead it might print 'too many spots to count' or something similar.
"""
def batch_mask(names, min=0.96):
    flat_file = open('flat_stats_' + str(min) + '.txt', 'w')
    pixVals = [1.06, 1.04, 1.02, 0.98, 0.96, 0.94, 0.92, 0.9, 0.85, 0.8, 0.75, 0.7]
    files = glob.glob('visao_flat_2*')
    print(files)

    for i, f in enumerate(files):
        pixList = dust_vals(f)
        flat_file.write(f + '\n')
        flat_file.write(str(pixVals) + '\n')
        flat_file.write(str(pixList) + '\n')
        mask = dust_mask(f, min, 2)
        fits.writeto(str(min) + '_level/visao_flat_mask_' + names[i] + '_min' + str(min) + '.fits', mask, clobber=True)
        flat_file.write('Median: ' + str(bn.nanmedian(mask)) + '\t')
        flat_file.write('Mean: ' + str(bn.nanmean(mask)) + '\t')
        flat_file.write('Standard deviation: ' + str(bn.nanstd(mask)) + '\n')

        Nans = np.where(np.isnan(mask))
        flat_file.write('Num NaN: ' + str(Nans[0].size) + '\n')
        print('Processing flat ' + str(i + 1) + ' of ' + str(len(files)))
        temp_im = np.array(mask)[20:1000, 20:1000]
        num_spots = 0
        dust_spots = []
        for y in range(0, temp_im.shape[0]):
            for x in range(0, temp_im.shape[1]):
                if count_spots(temp_im, temp_im[y][x], y, x):
                    num_spots += 1
                    dust_spots.append((y + 20, x + 20))
        flat_file.write('Number of dust spots: ' + str(num_spots) + '\n\n')

    flat_file.close()


"""
##################################################################################################################
NAME: batch_compare

PURPOSE/DESCRIPTION: Takes a list of masked flat images and numbers NaN pixels based on their frequency. For example
                    a pixel that appears in 5 of the images will be given 5 counts, a pixel that only appears in one
                    image will be given one count. This helps the study of the stability of the flat spots.

INPUTS:
    imList:     list of strings, names of images to analyze
    name:       string, name to save final fits image as

OUTPUTS:
    'spot_propagation_name.fits':   a fits image made up of predominantly zeros with other small values to describe
                                    stability of the spots

HISTORY:
    Created: 2016-08-19 by Wyatt Mullen, wmullen1@stanford.edu
"""
def batch_compare(imList, name):
    temp = fits.open(imList[0])
    outArr = np.zeros((temp[0].shape))
    temp.close()
    for i, im in enumerate(imList):
        im = fits.open(im)
        nans = np.where(np.isnan(im[0].data))
        outArr[nans] += 1
        im.close()

    fits.writeto('spot_propagation_' + name + '.fits', outArr, clobber = True)
