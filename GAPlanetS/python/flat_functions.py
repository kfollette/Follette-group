import numpy as np
from astropy.io import fits
import bottleneck as bn
import math
import matplotlib.pyplot as plt
import glob
import visao_attempt1 as vis
from scipy.ndimage.filters import uniform_filter
from scipy.ndimage.filters import median_filter
from astropy.stats import sigma_clip

"""
NAME: raw_flat_process



"""
def raw_flat_process(directory='raw', dark=None, blur=41):
    searchstr = directory + '/*.fits'
    files = glob.glob(searchstr)
    if dark == None:
        imList = vis.visao_inventory()
        if not imList['dark_imlist']:
            raise FileNotFoundError('Please specify a master dark to use.')
        else:
            vis.visao_dark()
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
        bot_med = np.median(raw_im[511:512, 0:1024])
        top_med = np.median(raw_im[512:513, 0:1024])
        raw_im[0:512, 0:1024] += top_med - bot_med
        #smoothing by median filter
        smooth_im = uniform_filter(smooth_im, size=blur)
        #smooth_im = median_filter(raw_im, size=41)
        #dividing by smoothed image
        raw_im = raw_im / smooth_im
        #adding image to cube of images
        flat_cube[im_num,:,:] = raw_im
        im_num += 1
        print('Completed image ' + str(im_num) + ' of ' + str(len(files)))

    #taking median of image cube for each pixel
    master_flat = np.median(flat_cube, axis=0)
    fits.writeto('visao_flat_' + files[0][8:16] + '_median' + str(blur) + '-2.fits', master_flat, clobber=True)


"""

        #print('After bot-top subtraction median: ' + str(np.median(raw_im)) + ' mean: ' + str(np.mean(raw_im)) +
              #' std: ' + str(np.std(raw_im)))
        #smooth_im = raw_im
                if im_num < 4:
            fits.writeto('smooth_im-' + str(im_num + 1) + '.fits', smooth_im, clobber=True)
                #smooth_im = uniform_filter(smooth_im, size=blur)
        #smooth_im = sci.gaussian_filter(smooth_im, sigma=10)
        #print('After smooth median: ' + str(np.median(smooth_im)) + ' mean: ' + str(np.mean(smooth_im)) +
              #' std: ' + str(np.std(smooth_im)))
                      if im_num < 3:
            fits.writeto('divided_im' + str(im_num + 1) + '.fits', raw_im, clobber=True)

    #master_flat = np.mean(sigma_clip(flat_cube, sigma=3), axis=0)
    #print(line)
    #print(np.mean(line, axis=0))
    #clip1 = np.where(master_flat > np.mean(master_flat) + 3 * np.std(master_flat) or
                     #master_flat < np.mean(master_flat) + 3 * np.std(master_flat))
   # master_flat[clip1] = np.nan

"""
def box_blur(image, size):
    image = np.array(image)
    temp_image = image
    for y in range(size//2, image.shape[0] - size//2 - 2):
        for x in range(size//2, image.shape[1] - size//2 - 2):
            temp_image[y][x] = np.nanmedian(image[-size//2 + y:size//2 + y, -size//2 + x:size//2 + x])
            #print(np.nanmedian(image[-size//2 + y:size//2 + y, -size//2 + x:size//2 + x]))
        #print(temp_image[y][x])
    print(cloriform)
    return np.array(temp_image).tolist()

"""


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
            fits.writeto('visao_flat_' + str(date) +'_median41_nansmask-2.fits', mask_image, clobber=True)
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
NAME: dust_mask

PURPOSE: Create a mask around dust spots for flat images where the dust spots are masked as NaNs

DESCRIPTION: The function takes an image, a minimum opacity (out of 1.00), and a radius in pixels. The bad values are
            found first and replaced with nans. Then the program loops through the rest of the image and replaces any
            pixels within a circular radius 'rad' with the value NaN. The mask is then returned.

INPUTS:
    image:  flat image to be masked
    minval: minimum value that will ot be masked in the flat image. Default is 0.98
    rad:    radius of pixels around center mask to make sure that dust spots will have minimal impact. Default is 2

OUTPUTS:
    mask:   2-d masked image array

HISTORY:
    CREATED: 2016-07-20 by Wyatt Mullen, wmullen1@stanford.edu
"""
def dust_mask(image, minval=0.98, rad=2, sem='2014A'):
    #masking all values below the minval
    image = fits.open(image)[0] #added
    nans = np.where(image.data < minval)
    data = np.where(image.data >= minval)
    dimy = image.shape[0]
    dimx = image.shape[1]
    mask = np.zeros(image.shape)
    temp_mask = np.zeros(image.shape)
    mask[nans] = np.NaN
    temp_mask[nans] = np.NaN
    #mask[data] = image.data[data]
    mask[data] = 1 #added

    #putting pixel radius around every dust spot
    #temp_mask = mask
    for y in range(0,dimy):
        for x in range(0,dimx):
            if np.isnan(temp_mask[y][x]):
                for i in range(-2, 3):
                    for j in range(-2, 3):
                        if 0 <= x + j and x+j < dimx and 0 <= y+i and y+i < dimy:
                            mask[y+i][x+j] = np.NaN
    fits.writeto(sem + '_mask_ones_' + str(minval) + '.fits', mask, clobber=True) #added
    #return mask

"""


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



def count_spots(temp_im, pixel, y, x):
    #spot = False
    if not np.isnan(pixel) or y < 0 or x < 0 or y + 2 > temp_im.shape[0] or x + 2 > temp_im.shape[1]:
        return False
    else:
        temp_im[y][x] = 100.0
        for i in range(-1,2):
            for j in range(-1,2):
                count_spots(temp_im, temp_im[y+i][x+j], y + i, x + j)
        return True

def flat_compare(im1, im2, year):
    im1 = fits.open(im1)[0]
    im2 = fits.open(im2)[0]
    im_out = np.zeros(im1.shape)
    ones = np.where((np.isnan(im1.data)) & (~np.isnan(im2.data)))

    twos = np.where((np.isnan(im1.data)) & (np.isnan(im2.data)))
    print(twos)
    threes = np.where((np.isnan(im2.data)) & (~np.isnan(im1.data)))
    im_out[ones] = 1
    im_out[twos] = 2
    im_out[threes] = 3

    fits.writeto('flat_compare_' + str(year) + '.fits', im_out, clobber=True)

def dust_vals(im):
    im = fits.open(im)[0]
    pixVals = [1.06, 1.04, 1.02, 0.98, 0.96, 0.94, 0.92, 0.9, 0.85, 0.8, 0.75, 0.7]
    for val in pixVals:
        if val > 1.00:
            a = np.where(im.data > val)
            print('NumPix greater than ' + str(val) + ': ' + str(len(a[0])))
        else:
            a = np.where(im.data < val)
            print('NumPix less than ' + str(val) + ': ' + str(len(a[0])))