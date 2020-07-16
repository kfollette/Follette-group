import numpy as np
import matplotlib.pyplot as plt 
from scipy import optimize
from scipy import stats
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import UnivariateSpline
import sys
import glob


# How to run: python ComboMassTeff_Interp_B98-B15-MESA_1myr-2myr.py on filename_with_teffULs+LLs 2 
plotcheck = str(sys.argv[1])

targets_to_convert = str(sys.argv[2])

ages = str(sys.argv[3])

models = glob.glob('./Baraffe*txt')
mesamodels = glob.glob('./MESA_*.txt') #format is the same: Mass    Teff    LogLum  Logg    Radius
for ii in mesamodels:
    models.append(ii)


print(models)

#targets_to_convert = 'AndrewsClassII_SpectralTypes_DerivedTeff_Results.txt'


##########
# Unpack the target file, which is an output from Luhman2003-KenyonHartmann_SpTy_Teff_relation.py
# That output file should be in the following form:
# Target SpTy NumericalSpTy Teff
##########

TargetName, TargSpTy, TargSpTyNum, TargTeff, TargTeff_earlylim, TargTeff_latelim = np.genfromtxt(targets_to_convert, skip_header=1, delimiter=';', dtype='str').T

# Make the Teff into an array of floats
TargTeff = [float(x) for x in TargTeff]
TargTeff = np.array(TargTeff)

#Apply interpolation for all values within the 2015 Baraffe models (1.4Msun -- some of the Teffs are higher than this and will crash the task, so we'll now restrict it to everything less than 4658 K)
TargTeff_lim = TargTeff[np.where(TargTeff < 4658.)]

# And match the target name, etc, arrays to correspond to these values:
TargetName_lim = TargetName[np.where(TargTeff < 4658.)]
TargSpTy_lim = TargSpTy[np.where(TargTeff < 4658.)]
TargSpTyNum_lim = TargSpTyNum[np.where(TargTeff < 4658.)]

plt.clf()


# loop over the five models (three if excluding 1 myr popn)

if ages == "2":
    linestyles = ['-', '--', '-.']
    models = ['Baraffe1998_2Myr.txt', 'Baraffe2015_2Myr.txt', 'MESA_2Myr.txt']
elif ages == "1":
    linestyles = ['-', '--', '-.']
    models = ['Baraffe1998_1Myr.txt', 'Baraffe2015_1Myr.txt', 'MESA_1Myr.txt'] 
else:
    linestyles = ['-', '--', '-', '--', '-.']

#reverse the model list so the MESA models are plotted first
models.reverse()

for idx, modelname in enumerate(models):

    ##########
    # Import mass and teff data from isochrones
    ##########
    modelshortname = modelname.split('.')[0]
    labelname = modelshortname.split('_')[0] + " " +modelshortname.split('_')[1]

    mass, teff, loglum, logg, radius = np.loadtxt(modelname, skiprows=1).T


    ##########
    # Interpolating the data 
    ##########

    if modelname == "MESA_2Myr.txt":
        spl = LSQUnivariateSpline(teff, radius, [3401, 4617, 5052, 6028, 6350, 11550, 11912, 12110, 14108, 14500])
        spl_mesa = LSQUnivariateSpline(teff, radius, [3410, 4617, 5052, 6028, 6350, 11550, 11912, 12110, 14108, 14500])
        spl_lum = LSQUnivariateSpline(teff, loglum, [3298, 3355, 3462, 4068, 4646, 5219, 6277, 6319, 11738, 11943, 14123, 14281])
        spl_mesa_lum = LSQUnivariateSpline(teff, loglum, [3298, 3355, 3462, 4068, 4646, 5219, 6277, 6319, 11738, 11943, 14123, 14281])
        linecolor = 'green'    

    elif modelname == "MESA_1Myr.txt":
        spl = LSQUnivariateSpline(teff, radius, [2980, 3101, 3619, 3995, 4257, 4512, 4691, 5183])
        #spl = LSQUnivariateSpline(teff, mass, [3427, 3927, 4088, 4667, 5016, 5316, 5373, 5400, 5589][2:-1])
        spl_mesa = LSQUnivariateSpline(teff, radius, [3427, 3927, 4088, 4667, 5016, 5316, 5373, 5400, 5589])
        spl_lum = LSQUnivariateSpline(teff, loglum, [3298, 3355, 3462, 4068, 4646, 5219, 6277, 6319])
        spl_mesa_lum = LSQUnivariateSpline(teff, loglum, [3298, 3355, 3462, 4068, 4646, 5219, 6277, 6319])
        linecolor = 'green'                 

    elif modelname == "Baraffe2015_1Myr.txt":    
        #spl = LSQUnivariateSpline(teff, mass, [2864, 2896, 2897, 2898, 2908]) # these are the "knots" around 0.08 Msun, where the model is discontinuous
        spl = LSQUnivariateSpline(teff, radius, [2864, 2896, 2897, 2898, 2908]) # these are the "knots" around 0.08 Msun, where the model is discontinuous
        spl_b15_1myr = LSQUnivariateSpline(teff, radius, [2864, 2866, 2896, 2897, 2900, 2908])
        spl_lum = LSQUnivariateSpline(teff, loglum, [2599, 2869, 2898, 2902, 2911, 3437])
        spl_b15_lum = LSQUnivariateSpline(teff, loglum, [2599, 2869, 2898, 2902, 2911, 3437])
        linecolor = 'r'

    elif modelname == "Baraffe1998_1Myr.txt":    
        #spl = LSQUnivariateSpline(teff, radius, [2947, 3024, 3073, 3119, 3427, 3555, 3898])
        #spl_b98_1myr = LSQUnivariateSpline(teff, radius, [2947, 3024, 3073, 3119, 3427, 3555, 3898])
        spl_lum = LSQUnivariateSpline(teff, loglum, [2863, 3078, 3166, 3767])
        #spl_b98_lum = LSQUnivariateSpline(teff, loglum, [2545, 2707, 2796, 2863, 2921, 3005, 3078, 3166, 3509, 3767])
        spl_b98_lum = LSQUnivariateSpline(teff, loglum, [2863, 3078, 3166, 3767])
        linecolor = 'b'

    elif modelname =="Baraffe1998_2Myr.txt":
        #spl = LSQUnivariateSpline(teff, radius, [2928, 2937, 2947, 3038, 3427, 3555, 3898])
        #spl_b98_2myr = LSQUnivariateSpline(teff, radius, [2928, 2937, 2947, 3038, 3427, 3555, 3898])
        spl_lum = LSQUnivariateSpline(teff, loglum, [2975, 3050, 4000])
        spl_b98_lum = LSQUnivariateSpline(teff, loglum, [2975, 3050, 4000])
        linecolor = 'b'

    elif modelname == "Baraffe2015_2Myr.txt":
        spl = LSQUnivariateSpline(teff, radius, [2864, 2896, 2897, 2898, 2908, 2935]) # these are the "knots" around 0.08 Msun, where the model is discontinuous
        spl_b15_2myr = LSQUnivariateSpline(teff, radius, [2864, 2896, 2897, 2898, 2908, 2935]) 
        spl_lum = LSQUnivariateSpline(teff, loglum, [2458, 2714, 2913, 2974, 3231])
        spl_b15_lum = LSQUnivariateSpline(teff, loglum, [2458, 2714, 2913, 2974, 3231])
        linecolor = 'r'           

    else:
        print("error!")

    DerivedRadii_lim = spl(TargTeff_lim)


    newteff = np.arange(min(teff),max(teff),0.1)
    newradius = spl(newteff)
    newloglum = spl_lum(newteff)

    plt.figure(1, figsize = (8,8))
    plt.subplot(211)
    plt.plot(teff, radius, 'x', mfc=linecolor, mec=linecolor)
    plt.plot(newteff, newradius, linestyle=linestyles[idx], color=linecolor, label = labelname)
    plt.xlim(1500, 10000)
    plt.ylim(0, 20)
    plt.xlabel(r'T$_{eff}$ (K)')
    plt.ylabel(r'Radius (R$_{J}$)')
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend(loc=1, numpoints=1)
    

    plt.subplot(212)
    plt.plot(teff, loglum, 'x', mfc=linecolor, mec=linecolor)
    plt.plot(newteff, newloglum, linestyle=linestyles[idx], color=linecolor, label = labelname)
    plt.xlabel(r'T$_{eff}$ (K)')
    plt.ylabel(r'Log Luminosity (L$_{\odot}$)')
    plt.ylim(-3, 2.5)
    plt.xlim(1500, 10000)
    plt.xscale('log')
    plt.legend(loc=4, numpoints=1)
    


plt.subplots_adjust(hspace=0.35)

plt.savefig('B15-B98-MESA_Mass_Teff_Relation_%.0fMyr.png' % float(ages), dpi=200)

if plotcheck == "on":
    plt.show()




