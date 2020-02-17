import sys
from math import *
import numpy as np
from matplotlib import *
from pylab import *
import pandas as pd
from astropy.io import fits
from scipy.interpolate import interp1d 


"""
# Actual usage:
# python pradialprofile.py image parallax(mas) pixelscale(arcsec/px) filter AbsMag 

N.B. Eventually will have to add xcenter and ycenter keywords, but for now grab those vals from CFHT pixel headers (change for VLT)

Last updated: 
2020-02-17 - (KWD) General cleaning and revising.
2019-07-21 - (KWD) Quick fixes for print statements for Python 3.
2017-08-14 - (KWD) Cleaned up some of the code and old comments, path to evolutionary models no longer hard-coded, but 
				text file needs to be in the same directory to work. (MassMagRelation_5Gyr_JHK.txt) 
2013-02-02 - (KWD) Started code (approximately). 

"""



im = fits.getdata(sys.argv[1]+'.fits')
fitsheader = fits.getheader(sys.argv[1]+'.fits')


pixelscale = float(sys.argv[3])

# distance in parsecs
dist = 1000.0/float(sys.argv[2])

# use set 1 or set 2 depending on origin of data (CFHT and MMT has CRPIX1A vals, VLT does not)
xcenter = fitsheader['CRPIX1A']
ycenter = fitsheader['CRPIX2A']

#xcenter = float(sys.argv[6])
#ycenter = float(sys.argv[7])

print("\n Target is",sys.argv[1])
print("Xcenter is",xcenter)
print("Ycenter is",ycenter)

# Define box in which to search for stellar peak
xmin = int(np.floor(xcenter-5))
xmax = int(np.ceil(xcenter+5))
ymin = int(np.floor(ycenter-5))
ymax = int(np.ceil(ycenter+5))

# Extract peak pixel value in image using Python (as opposed to pyraf...)
peak = np.nanmax(im[ymin:ymax,xmin:xmax])
peakcoords = np.where(im==peak)

print(f"The peak pixel count is {peak}")

# Define new xcenter, ycenter from that peak value
xcenter = peakcoords[1][0]
ycenter = peakcoords[0][0]


# Work out maximum radius for this frame
s = np.shape(im)
x = np.zeros((s[0],s[1]))
y = np.zeros((s[0],s[1]))

for i in np.arange(0,s[0]):
	y[i,:] = i

for i in np.arange(0,s[1]):
	x[:,i] = i


# Define x coordinate relative to star coordinate in original image
x = x-xcenter

# Define y coordinate relative to star coordinate in original image
y = y-ycenter

# Define r as look up table of radii for each pixel from star center
r = np.sqrt((x**2)+(y**2))


# Retrieve maximum radius (out to edge of detector) from command line
maxrad = np.max(r)-10.0

# --------- Calculation of the separation coverage in AU --------- #
# Define a large image to determine the radial extent of the sensitivity:
x_lg = np.zeros((s[0]+1000,s[1]+1000))
y_lg = np.zeros((s[0]+1000,s[1]+1000))
for i in np.arange(0,s[0]+1000):
	y_lg[i,:] = i

for i in np.arange(0,s[1]+1000):
	x_lg[:,i] = i

x_lg = x_lg - xcenter - 500
y_lg = y_lg - ycenter - 500

r_lg = np.sqrt((x_lg**2)+(y_lg**2))

foo = np.zeros(np.shape(r_lg))
foo[:] = np.nan 
foonew = foo
# embed the image in a sea of NaNs
foonew[500:(500+np.shape(im)[0]),500:(500+np.shape(im)[1])] = im

minarea = np.min(np.shape(im))

#for i in np.linspace(np.min(np.shape(im))-200,np.max(r_lg),50):
for i in np.linspace(300,np.max(r_lg),100):
	contained = np.where(r_lg < (i+5))
	num_elements = np.shape(contained)[1]
	num_of_not_nans = float(np.sum(np.isfinite(foonew[contained])))
	#print i
	if (num_of_not_nans/num_elements) < 0.95:
		maxpixeldistance = i 
		break
	#break

maxseparation_as = pixelscale*maxpixeldistance
maxseparation_au = maxseparation_as*dist

print(f"The maximum radius with {num_of_not_nans*100/num_elements} coverage is {maxpixeldistance} pixels, corresponding to {maxseparation_as} arcseconds or {maxseparation_au} AU") 

del foonew,r_lg,x_lg,y_lg,foo









#print "Calculating radial profile... Please wait..."

# Measure radial profile though stsdas.plot's pradprof command
#radprofdata = iraf.plot.pradprof("%s" % sys.argv[1],xcenter,ycenter,radius=maxrad,Stdout=1)

# Read in data from iraf and find the sorting order of the radius column 
#data = np.loadtxt(radprofdata,dtype=float)
#print data
#order  = data[:,0].argsort()

# Sort data array by radius
#sorteddata = np.take(data,order,0)
#sorteddata = data

#lenarray = len(sorteddata)
#peak = float(sys.argv[2])

# rename things
#newsortdata = sorteddata

# convert second column of array -- counts -- into delta magnitudes, using peak pixel value of star provided on command line
#newsortdata[:,1] = 2.5*np.log10(peak/sorteddata[:,1])

# round each value in radius column up to whole pixel value
#sorteddata[:,0] = np.ceil(sorteddata[:,0])


# Define empty lists to place calculated values
radiuslist = list()
medianlist = list()
fivesigma = list()
threesigma = list()

print("Calculating sensitivity limits for each annulus from core...")

nanindex = np.isnan(im)
im[nanindex] = 0

nonzeroind=(im != 0)

for i in np.arange(0,np.floor(maxrad)):
	#index = np.where((r >= i) & (r < (i+5))) 
	#len1 = index
	index = np.where((r >= i) & (r < (i+5)) & nonzeroind) #explicitly ignore NaN pixels
	#print i
	sigma = np.std(im[index],ddof=1)
	median = np.median(im[index])
	new5sigma = 2.5*np.log10(peak/(5.0*sigma))
	new3sigma = 2.5*np.log10(peak/(3.0*sigma))
	if (np.isnan(new5sigma) | np.isinf(new5sigma)):
		new5sigma = 0
	if (np.isnan(new3sigma) | np.isinf(new3sigma)):
		new3sigma = 0
	radiuslist.append(i+2.5)
	medianlist.append(median)
	fivesigma.append(new5sigma)
	threesigma.append(new3sigma)

'''for i in np.arange(0,maxrad):
	index = np.where((newsortdata[:,0] >= i) & (newsortdata[:,0] < (i+5)))
	sigma = np.std(newsortdata[:,1][index],ddof=1)
	median = np.median(newsortdata[:,1][index])
	#print i+2.5,median,sigma
	# i+2.5 is the middle of the annulus; i is the interior edge; i+5 is the exterior edge
	new5sigma = 2.5*np.log10(peak/(5.0*sigma))
	new3sigma = 2.5*np.log10(peak/(3.0*sigma))
	radiuslist.append(i+2.5)
	medianlist.append(median)
	fivesigma.append(new5sigma)
	threesigma.append(new3sigma)'''


print('Finished calculating sensitivity!')
	


#threesig = [x*3.0 for x in sigmalist]



# convert list of radii from pixels to arcseconds
seplist = [x*pixelscale for x in radiuslist]

PrimAbsMag = float(sys.argv[5])
massmag = pd.read_csv("MassMagRelation_5Gyr_JHK.txt",delim_whitespace=True)
AbsMag5 = [x+PrimAbsMag  for x in fivesigma]
AbsMag3 = [x+PrimAbsMag  for x in threesigma]



# Calculate values at bottom of main sequence for various filters at 3 sigma confidence
if sys.argv[4] == 'K':
	MSdeltamag = 10 - float(sys.argv[5])
	print("Delta mag. required to reach bottom of Main Sequence:",MSdeltamag)
	f_K = interp1d(massmag['Mag'],massmag['KMass'],kind='cubic')
	prim_mass = float(f_K(PrimAbsMag))
	massK_3 = f_K(AbsMag3)
	massK_5 = f_K(AbsMag5)
	massratio_3 = [x/prim_mass for x in massK_3]
	massratio_5 = [x/prim_mass for x in massK_5]
	Msol_3 = massK_3
	Msol_5 = massK_5

if sys.argv[4] == 'H':
	MSdeltamag = 10.3 - float(sys.argv[5])
	print("Delta mag. required to reach bottom of Main Sequence:",MSdeltamag)
	f_H = interp1d(massmag['Mag'],massmag['HMass'],kind='cubic')
	prim_mass = float(f_H(PrimAbsMag))
	massH_3 = f_H(AbsMag3)
	massH_5 = f_H(AbsMag5)
	massratio_3 = [x/prim_mass for x in massH_3]
	massratio_5 = [x/prim_mass for x in massH_5]
	Msol_3 = massH_3
	Msol_5 = massH_5


if sys.argv[4] == 'J':
	MSdeltamag = 11.0 - float(sys.argv[5])
	print("Delta mag. required to reach bottom of Main Sequence:",MSdeltamag)
	f_J = interp1d(massmag['Mag'],massmag['JMass'],kind='cubic')
	prim_mass = float(f_J(PrimAbsMag))
	massJ_3 = f_J(AbsMag3)
	massJ_5 = f_J(AbsMag5)
	massratio_3 = [x/prim_mass for x in massJ_3]
	massratio_5 = [x/prim_mass for x in massJ_5]
	Msol_3 = massJ_3
	Msol_5 = massJ_5

# from 5 Gyr Baraffe model
if sys.argv[4] == 'L':
	MSdeltamag = 9.97 - float(sys.argv[5])
	print("Delta mag. required to reach bottom of Main Sequence:",MSdeltamag)
	f_L = interp1d(massmag['Mag'],massmag['LMass'],kind='cubic')
	prim_mass = float(f_L(PrimAbsMag))
	massL_3 = f_L(AbsMag3)
	massL_5 = f_L(AbsMag5)
	massratio_3 = [x/prim_mass for x in massL_3]
	massratio_5 = [x/prim_mass for x in massL_5]
	Msol_3 = massL_3
	Msol_5 = massL_5

max_xval = max(seplist)+10.0



# ----------- Plot 1: Delta Mag vs. Sep (arcsec) ----------- #
#plt.errorbar(seplist,medianlist,yerr=threesig,fmt='bo')
plt.semilogx(seplist,fivesigma,'b--', label=r'$5 \sigma$ Sensitivity')
plt.semilogx(seplist,threesigma,'r-.',label=r'$3 \sigma$ Sensitivity')
#plt.xscale('log')


plt.xlabel(r'Radius from star center (arcsec)')
plt.ylabel(r'$\Delta$ mag')
plt.ylim(np.ceil(max(threesigma)),0)

# Change upper limit for various instruments
plt.xlim(0.1,50.0)
plt.legend(numpoints=1)
plt.title(sys.argv[1])
plt.hlines(MSdeltamag,0.0001,100)
plt.text(1,MSdeltamag-0.25,r'Bottom of Main Sequence: $\Delta$M = %.3f' % MSdeltamag)
name = sys.argv[1]

plt.gcf()
plt.savefig("%s"%name + ".contrast_arcsec.png") 

#plt.show()



plt.clf()

# ----------- Plot 2: Delta Mag vs. Sep (AU) ----------- #
# Second plot in AU instead


AUlist = [x*pixelscale*dist for x in radiuslist]

plt.semilogx(AUlist,fivesigma,'b--', label=r'$5 \sigma$ Sensitivity')
plt.semilogx(AUlist,threesigma,'r-.',label=r'$3 \sigma$ Sensitivity')
#plt.xscale('log')

plt.xlabel(r'Radius from star center (AU)')
plt.ylabel(r'$\Delta$ mag')
plt.ylim(np.ceil(max(threesigma)),0)

# Change upper limit for various instruments
plt.xlim(0.5,1000.0)
plt.legend(numpoints=1)
plt.title(sys.argv[1])
plt.hlines(MSdeltamag,0.0001,1000)
plt.text(1,MSdeltamag-0.25,r'Bottom of Main Sequence: $\Delta$M = %.3f' % MSdeltamag)
name = sys.argv[1]

plt.gcf()
plt.savefig("%s"%name + ".contrast_AU.png") 

#plt.show()
plt.clf()


# ----------- Plot 3: Mass in Msol vs. Sep (AU)----------- #
# Third plot in MASS vs AU sep
plt.semilogx(AUlist,Msol_5,'b--', label=r'$5 \sigma$ Sensitivity')
plt.semilogx(AUlist,Msol_3,'r-.',label=r'$3 \sigma$ Sensitivity')
#plt.xscale('log')

plt.xlabel(r'Projected radius from star center (AU)')
plt.ylabel(r'Mass ($M_{\odot}$)')
#plt.ylim(0,1)

# Change upper limit for various instruments
plt.xlim(0.5,1000.0)
plt.legend(numpoints=1)
plt.title(sys.argv[1])
#plt.hlines(MSdeltamag,0.0001,1000)
#plt.text(1,MSdeltamag-0.25,r'Bottom of Main Sequence: $\Delta$M = %.3f' % MSdeltamag)
name = sys.argv[1]

plt.gcf()
plt.savefig("%s"%name + ".contrast_MassMsol_AU.png") 

plt.clf()


# ----------- Plot 4: Mass Ratio vs. Sep (AU) ----------- #
# Fourth plot in mass ratio vs AU sep
plt.semilogx(AUlist,massratio_5,'b--', label=r'$5 \sigma$ Sensitivity')
plt.semilogx(AUlist,massratio_3,'r-.',label=r'$3 \sigma$ Sensitivity')
#plt.xscale('log')

plt.xlabel(r'Radius from star center (AU)')
plt.ylabel(r'Mass Ratio $M_{2} / M_{1}$')
plt.ylim(0,1)

# Change upper limit for various instruments
plt.xlim(0.5,1000.0)
plt.legend(numpoints=1)
plt.title(sys.argv[1])
#plt.hlines(MSdeltamag,0.0001,1000)
#plt.text(1,MSdeltamag-0.25,r'Bottom of Main Sequence: $\Delta$M = %.3f' % MSdeltamag)
name = sys.argv[1]

plt.gcf()
plt.savefig("%s"%name + ".contrast_MassRatioAU.png") 


# ----------- Values for the bottom of the Main Sequence ----------- #
number = len(threesigma)

# array of values for bottom of MS in the particular filter of observation
MSmag = np.zeros(number) + MSdeltamag


# find approx position in delta mag array where bottom of main sequence is reached
pos_5 = min(range(len(fivesigma)), key=lambda i: abs(fivesigma[i]-MSdeltamag))

pos_3 = min(range(len(threesigma)), key=lambda i: abs(threesigma[i]-MSdeltamag))

# find the value in AU at this position (5 sigma confidence)
val_5 = AUlist[pos_5]

# find the value in AU at this position (3 sigma confidence)
val_3 = AUlist[pos_3]

MSseparation_5 = np.zeros(number) + val_5
MSseparation_3 = np.zeros(number) + val_3

maxseparation_as = np.zeros(number) + maxseparation_as
maxseparation_au = np.zeros(number) + maxseparation_au




# "fivesigma" and "threesigma" are both in delta magnitudes
fileoutput1 = zip(seplist,AUlist,medianlist,fivesigma,threesigma,MSmag,MSseparation_5,MSseparation_3,Msol_5,Msol_3,massratio_5,massratio_3,maxseparation_as,maxseparation_au)

df = pd.DataFrame(fileoutput1,columns=['Radius(arcsec)','Radius(AU)','MedianPxCt','FiveSigma','ThreeSigma','DeltaMagtoMS','AUSepAtMS_5sig','AUSepAtMS_3sig','MSol_5sig','MSol_3sig','MassRatio_5sig','MassRatio_3sig','MaxSepArcSec','MaxSepAU'])

df.to_csv("%s"%name + '.contrast_wmass.csv')


