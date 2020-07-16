import numpy as np
import matplotlib.pyplot as plt 
from scipy import optimize
from scipy import stats
from scipy.interpolate import UnivariateSpline
import sys


# How to run: python Herczeg_Hillenbrand_SpTy_Teff_Relation.py "on"(or "off") listofSpTys_toconvert
plotcheck = str(sys.argv[1])

targets_to_convert = str(sys.argv[2])

# Now read from the list of targets file; should be in a tab-delimited format of:
# TargetName SpTy SpTyNum
#
# (these are calculated in /Users/Kim/Research/ALMA/TaurusContinuumDraft/TempSpTyRelation.xlsx)


TargetName, TargSpTy, TargSpTyNum = np.genfromtxt(targets_to_convert, skip_header=1, dtype='str', delimiter='\t').T

# Make the SpTyNum into an array
TargSpTyNum = [float(x) for x in TargSpTyNum]
TargSpTyNum = np.array(TargSpTyNum)

TargSpTyNum_earlylim = TargSpTyNum-0.5

TargSpTyNum_latelim = TargSpTyNum+0.5

# Mix of temperature scales from Herczec & Hillenbrand 2014
SpTyNum, Teff, SpTy = np.genfromtxt('HerczegHillenbrand_SpTyTeff_Numericized.txt', skip_header=1, dtype='str').T 

SpTyNum, Teff = [float(x) for x in SpTyNum], [float(y) for y in Teff]

# To decode the numerical spectral types to actual types:
# SpectralClass   Value   Prefix
# B   0-9 0
# A   10-19   1
# F   20-29   2
# G   30-39   3
# K   40-48   4
# M   50-59   5


# We'll do a spline interpolation over these:

spl = UnivariateSpline(SpTyNum, Teff) 


# And apply the interpolation to the dataset in question:
derived_teffs = spl(TargSpTyNum)


# Apply the interpolation to the lower and upper limits for the possible spectral types:
derived_teffs_earlylim = spl(TargSpTyNum_earlylim)
derived_teffs_latelim = spl(TargSpTyNum_latelim)

# Write to a file:
to_txt = np.c_[(TargetName, TargSpTy, TargSpTyNum, derived_teffs, derived_teffs_earlylim, derived_teffs_latelim)]
np.savetxt(targets_to_convert.split('.')[0] + '_DerivedTeff_Results_HH14.txt', to_txt, delimiter=';', fmt='%s')


# and plot this, if desired:
if plotcheck == "on":
    plt.clf()

    type_spacing = np.arange(9, max(SpTyNum), 0.1)

    newteff = spl(type_spacing)
    plt.plot(SpTyNum, Teff, 'k.', type_spacing, newteff, 'b-')
    plt.xlabel('Spectral Type')
    plt.ylabel(r'T$_{eff}$ (K)')
    plt.xticks(SpTyNum[::3], SpTy[::3])
    plt.xlim(23,60)
    plt.ylim(2000,8000)
    plt.savefig('All_Spty_Temps_HH14.png', dpi=200)
    plt.show()


