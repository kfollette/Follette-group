#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

import numpy as np
import scipy as sp
import scipy.ndimage
from astropy.io import fits
import sys
import numpy.ma as ma
import math
import pidly
import os
import glob

"""
fitsinfo_vlt.py: Script to extract/view FITS header information from a list of files. 

Usage: python fitsinfo_vlt.py

History:
(v1) 2018-06-01 - KWD: Updated from original script for Python 3 compatibility and astropy updates.
"""


# Create a list of all of the fits files in the current directory:
caliblist = glob.glob('*.fits')




def fileinfo(caliblist):
	'''
    fileinfo
    ----------
    Routine to return information for each FITS file in a given list of filenames. Function assumes that 
    the header keywords are specific to VLT/NaCo data.
        
    Inputs
    ----------
    caliblist      : (list) list of fits filenames, each in string format
    
    Returns
    ----------
    (print output) : header keywords and corresponding values from the input fits file, printed to terminal
    
    headerinfo.txt : (txt file) contains output of keywords and values 

    Dependents
    ----------
    None
    '''

	n = len(caliblist)

	# Create a new text file to contain the output of the fileinfo function.
	outputfile = open('headerinfo.txt', 'w')

	outputfile.write("Filename; ESOFITSname; Object; Type; Exptime; Filter 1; Filter 2; Filter 3; Filter 4; Camera; ImageDimensions; PI-Name; ProgID; Date" + "\n")

	print("Filename; ESOFITSname; Object; Type; Exptime; Filter 1; Filter 2; Filter 3; Filter 4; Camera; ImageDimensions; PI-Name; ProgID; Date")

	for i in range(0,n):
		fitsfilename = caliblist[i].rstrip()
		fitsfile = fits.open(caliblist[i])
		header = fitsfile[0].header
		im = fitsfile[0].data
		shape = im.shape
		fitstype = header['HIERARCH ESO DPR TYPE']
		fitsname = header['ARCFILE']
		fitstime = header['EXPTIME']
		fitsfilter1 = header['HIERARCH ESO INS OPTI3 NAME']
		fitsfilter2 = header['HIERARCH ESO INS OPTI4 NAME']
		fitsfilter3 = header['HIERARCH ESO INS OPTI5 NAME']
		fitsfilter4 = header['HIERARCH ESO INS OPTI6 NAME']
		fitscamera = header['HIERARCH ESO INS OPTI7 ID']
		fitsPI = header['HIERARCH ESO OBS PI-COI NAME']
		fitsprogid = header ['HIERARCH ESO OBS PROG ID']
		obsdate = header['DATE']
		#fitsxval = header['CRPIX1A']
		#fitsyval = header['CRPIX2A']
		filename = caliblist[i].rstrip()
		#starx = header['CRPIX1A']
		#fitsrot = header['HIERARCH ESO ADA POSANG']
		#print fitsname,";", fitstype,";", fitstime,";", fitsfilter1,";", fitsfilter2,";", fitsfilter3,";", fitsfilter4,";", fitscamera,";",filename,";",shape
		#print fitsname,";", fitstype,";", fitstime,";", fitsfilter1,";", fitsfilter2,";", fitsfilter3,";", fitsfilter4,";", fitscamera,";", fitsxval,";", fitsyval,";",filename,";",shape
		
		# if fits object name exists:
		fitsob = header['OBJECT']
		print(fitsfilename,";",fitsname,";", fitsob,";", fitstype,";", fitstime,";", fitsfilter1,";", fitsfilter2,";", fitsfilter3,";", fitsfilter4,";", fitscamera,";",shape,";",fitsPI,';',fitsprogid,';',obsdate)

		# if fits object name does not exist:
		#print fitsname,";", fitstype,";", fitstime,";", fitsfilter1,";", fitsfilter2,";", fitsfilter3,";", fitsfilter4,";", fitscamera,";",filename,";",shape

		#print fitsname, fitstype, fitstime, fitsrot, fitsob

		# write output to file:
		outputfile.write(fitsfilename + ";" + fitsname + ";" + fitsob + ";" + fitstype + ";" + str(fitstime) + ";" + fitsfilter1 + ";" + fitsfilter2 + ";" + fitsfilter3 + ";" + fitsfilter4 + ";" + fitscamera + ";" + str(shape) + ";" + fitsPI + ';' + fitsprogid + ';' + obsdate + '\n')

fileinfo(caliblist)



