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
        spl = LSQUnivariateSpline(teff, mass, [3401, 4617, 5052, 6028, 6350, 11550, 11912, 12110, 14108, 14500])
        spl_mesa = LSQUnivariateSpline(teff, mass, [3410, 4617, 5052, 6028, 6350, 11550, 11912, 12110, 14108, 14500])
        spl_lum = LSQUnivariateSpline(teff, loglum, [3298, 3355, 3462, 4068, 4646, 5219, 6277, 6319, 11738, 11943, 14123, 14281])
        spl_mesa_lum = LSQUnivariateSpline(teff, loglum, [3298, 3355, 3462, 4068, 4646, 5219, 6277, 6319, 11738, 11943, 14123, 14281])
        linecolor = 'green'    

    elif modelname == "MESA_1Myr.txt":
        spl = LSQUnivariateSpline(teff, mass, [2980, 3101, 3619, 3995, 4257, 4512, 4691, 5183])
        #spl = LSQUnivariateSpline(teff, mass, [3427, 3927, 4088, 4667, 5016, 5316, 5373, 5400, 5589][2:-1])
        spl_mesa = LSQUnivariateSpline(teff, mass, [3427, 3927, 4088, 4667, 5016, 5316, 5373, 5400, 5589])
        spl_lum = LSQUnivariateSpline(teff, loglum, [3298, 3355, 3462, 4068, 4646, 5219, 6277, 6319])
        spl_mesa_lum = LSQUnivariateSpline(teff, loglum, [3298, 3355, 3462, 4068, 4646, 5219, 6277, 6319])
        linecolor = 'green'                 

    elif modelname == "Baraffe2015_1Myr.txt":    
        #spl = LSQUnivariateSpline(teff, mass, [2864, 2896, 2897, 2898, 2908]) # these are the "knots" around 0.08 Msun, where the model is discontinuous
        spl = LSQUnivariateSpline(teff, mass, [2864, 2896, 2897, 2898, 2908]) # these are the "knots" around 0.08 Msun, where the model is discontinuous
        spl_b15_1myr = LSQUnivariateSpline(teff, mass, [2864, 2866, 2896, 2897, 2900, 2908])
        spl_lum = LSQUnivariateSpline(teff, loglum, [2599, 2869, 2898, 2902, 2911, 3437])
        spl_b15_lum = LSQUnivariateSpline(teff, loglum, [2599, 2869, 2898, 2902, 2911, 3437])
        linecolor = 'r'

    elif modelname == "Baraffe1998_1Myr.txt":    
        spl = LSQUnivariateSpline(teff, mass, [2947, 3024, 3073, 3119, 3427, 3555, 3898])
        spl_b98_1myr = LSQUnivariateSpline(teff, mass, [2947, 3024, 3073, 3119, 3427, 3555, 3898])
        spl_lum = LSQUnivariateSpline(teff, loglum, [2863, 3078, 3166, 3767])
        #spl_b98_lum = LSQUnivariateSpline(teff, loglum, [2545, 2707, 2796, 2863, 2921, 3005, 3078, 3166, 3509, 3767])
        spl_b98_lum = LSQUnivariateSpline(teff, loglum, [2863, 3078, 3166, 3767])
        linecolor = 'b'

    elif modelname =="Baraffe1998_2Myr.txt":
        spl = LSQUnivariateSpline(teff, mass, [2928, 2937, 2947, 3038, 3427, 3555, 3898])
        spl_b98_2myr = LSQUnivariateSpline(teff, mass, [2928, 2937, 2947, 3038, 3427, 3555, 3898])
        spl_lum = LSQUnivariateSpline(teff, loglum, [2975, 3050, 4000])
        spl_b98_lum = LSQUnivariateSpline(teff, loglum, [2975, 3050, 4000])
        linecolor = 'b'

    elif modelname == "Baraffe2015_2Myr.txt":
        spl = LSQUnivariateSpline(teff, mass, [2864, 2896, 2897, 2898, 2908, 2935]) # these are the "knots" around 0.08 Msun, where the model is discontinuous
        spl_b15_2myr = LSQUnivariateSpline(teff, mass, [2864, 2896, 2897, 2898, 2908, 2935]) 
        spl_lum = LSQUnivariateSpline(teff, loglum, [2458, 2714, 2913, 2974, 3231])
        spl_b15_lum = LSQUnivariateSpline(teff, loglum, [2458, 2714, 2913, 2974, 3231])
        linecolor = 'r'           

    else:
        print("error!")

    DerivedMasses_lim = spl(TargTeff_lim)


    newteff = np.arange(min(teff),max(teff),0.1)
    newmass = spl(newteff)
    newloglum = spl_lum(newteff)

    plt.figure(1, figsize = (8,8))
    plt.subplot(211)
    plt.plot(teff, mass, 'x', mfc=linecolor, mec=linecolor)
    plt.plot(newteff, newmass, linestyle=linestyles[idx], color=linecolor, label = labelname)
    plt.xlim(1500, 10000)
    plt.ylim(0, 5.0)
    plt.xlabel(r'T$_\textnormal{eff}$ (K)')
    plt.ylabel(r'Mass (M$_{\odot}$)')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=4, numpoints=1)
    

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





'''
# Assuming we want to look at the 2 Myr models:
'''
if ages == "2":
    # Now write the new Baraffe+1998/2015 2-myr models and MESA 2-myr models to a new output file for central object mass:


    ##########
    # Unpack the target file, which is an output from Luhman2003-KenyonHartmann_SpTy_Teff_relation.py
    # That output file should be in the following form:
    # Target SpTy NumericalSpTy Teff
    ##########

    TargetName, TargSpTy, TargSpTyNum, TargTeff, TargTeff_earlylim, TargTeff_latelim = np.genfromtxt(targets_to_convert, delimiter=';', dtype='str').T

    TargetName_lim_combo, TargSpTy_lim_combo, TargSpTyNum_lim_combo, TargTeff_lim_combo, DerivedMasses_lim_combo, DerivedMasses_minlim_combo, DerivedMasses_maxlim_combo = [],[],[],[],[],[],[]

    # Make the Teff and its min and max values into an array of floats
    TargTeff = [float(x) for x in TargTeff]
    TargTeff = np.array(TargTeff)

    TargTeff_max = [float(x) for x in TargTeff_earlylim] # assuming earlier spectral types have larger Teffs...
    TargTeff_max = np.array(TargTeff_max)

    TargTeff_min = [float(x) for x in TargTeff_latelim] # assuming later spectral types have smaller Teffs...
    TargTeff_min = np.array(TargTeff_min)

    '''# First Baraffe 2015 only!'''
    #Apply interpolation for all values within the Baraffe models 
    TargTeff_lim = TargTeff[np.where(TargTeff < 4211.)]

    TargTeff_maxlim = TargTeff_max[np.where(TargTeff < 4211.)]

    TargTeff_minlim = TargTeff_min[np.where(TargTeff < 4211.)]

    # And match the target name, etc, arrays to correspond to these values:
    TargetName_lim = TargetName[np.where(TargTeff < 4211.)]
    TargSpTy_lim = TargSpTy[np.where(TargTeff < 4211.)]
    TargSpTyNum_lim = TargSpTyNum[np.where(TargTeff < 4211.)]

    DerivedMasses_lim = spl_b15_2myr(TargTeff_lim) # for the values
    DerivedMasses_maxlim = spl_b15_2myr(TargTeff_maxlim) # for the upper limit (earlier spty - hotter)
    DerivedMasses_minlim = spl_b15_2myr(TargTeff_minlim) # for the lower limit (later spty - cooler)

    DerivedLuminosities_lim = spl_b15_lum(TargTeff_lim)
    DerivedLuminosities_maxlim = spl_b15_lum(TargTeff_maxlim)
    DerivedLuminosities_minlim = spl_b15_lum(TargTeff_minlim)


    '''# Then Baraffe 1998 only!'''
    #Apply interpolation for all values within the Baraffe models (older)

    DerivedMasses_lim_b98 = spl_b98_2myr(TargTeff_lim) # for the values 
    DerivedMasses_maxlim_b98 = spl_b98_2myr(TargTeff_maxlim) # for the upper limit (earlier spty - hotter)
    DerivedMasses_minlim_b98 = spl_b98_2myr(TargTeff_minlim) # for the lower limit (later spty - cooler)

    DerivedLuminosities_lim_b98 = spl_b98_lum(TargTeff_lim)
    DerivedLuminosities_maxlim_b98 = spl_b98_lum(TargTeff_maxlim)
    DerivedLuminosities_minlim_b98 = spl_b98_lum(TargTeff_minlim)


    '''# Then MESA only!'''

    #Apply interpolation for higher mass targets:
    TargTeff_lim_mesa = TargTeff[np.where(TargTeff > 4211.)]

    TargTeff_maxlim_mesa = TargTeff_max[np.where(TargTeff > 4211.)]

    TargTeff_minlim_mesa = TargTeff_min[np.where(TargTeff > 4211.)]

    # And match the target name, etc, arrays to correspond to these values:
    TargetName_lim_mesa = TargetName[np.where(TargTeff > 4211.)]
    TargSpTy_lim_mesa = TargSpTy[np.where(TargTeff > 4211.)]
    TargSpTyNum_lim_mesa = TargSpTyNum[np.where(TargTeff > 4211.)]

    DerivedMasses_lim_mesa = spl_mesa(TargTeff_lim_mesa) # for the values
    DerivedMasses_maxlim_mesa = spl_mesa(TargTeff_maxlim_mesa) # for the upper limit (earlier spty - hotter)
    DerivedMasses_minlim_mesa = spl_mesa(TargTeff_minlim_mesa) # for the lower limit (later spty - cooler)

    DerivedLuminosities_lim_mesa = spl_mesa_lum(TargTeff_lim_mesa) # for the values
    DerivedLuminosities_maxlim_mesa = spl_mesa_lum(TargTeff_maxlim_mesa) # for the upper limit (earlier spty - hotter)
    DerivedLuminosities_minlim_mesa = spl_mesa_lum(TargTeff_minlim_mesa) # for the lower limit (later spty - cooler)


    '''# combine everything to write to output file'''

    TargetName_lim_combo = np.concatenate([TargetName_lim, TargetName_lim_mesa])
    TargSpTy_lim_combo = np.concatenate([TargSpTy_lim, TargSpTy_lim_mesa])
    TargSpTyNum_lim_combo = np.concatenate([TargSpTyNum_lim,TargSpTyNum_lim_mesa])
    TargTeff_lim_combo = np.concatenate([TargTeff_lim, TargTeff_lim_mesa])
    DerivedMasses_lim_combo = np.concatenate([DerivedMasses_lim, DerivedMasses_lim_mesa])
    DerivedMasses_minlim_combo = np.concatenate([DerivedMasses_minlim, DerivedMasses_minlim_mesa])
    DerivedMasses_maxlim_combo = np.concatenate([DerivedMasses_maxlim, DerivedMasses_maxlim_mesa])

    DerivedMasses_lim_b98_combo = np.concatenate([DerivedMasses_lim_b98, DerivedMasses_lim_mesa])
    DerivedMasses_minlim_b98_combo = np.concatenate([DerivedMasses_minlim_b98, DerivedMasses_minlim_mesa])
    DerivedMasses_maxlim_b98_combo = np.concatenate([DerivedMasses_maxlim_b98, DerivedMasses_maxlim_mesa])

    DerivedLuminosities_lim_combo = np.concatenate([DerivedLuminosities_lim, DerivedLuminosities_lim_mesa])
    DerivedLuminosities_minlim_combo = np.concatenate([DerivedLuminosities_minlim, DerivedLuminosities_minlim_mesa])
    DerivedLuminosities_maxlim_combo = np.concatenate([DerivedLuminosities_maxlim, DerivedLuminosities_maxlim_mesa])

    DerivedLuminosities_lim_b98_combo = np.concatenate([DerivedLuminosities_lim_b98, DerivedLuminosities_lim_mesa])
    DerivedLuminosities_minlim_b98_combo = np.concatenate([DerivedLuminosities_minlim_b98, DerivedLuminosities_minlim_mesa])
    DerivedLuminosities_maxlim_b98_combo = np.concatenate([DerivedLuminosities_maxlim_b98, DerivedLuminosities_maxlim_mesa])


    # and write the combination to file:

    to_txt = zip(TargetName_lim_combo, TargSpTy_lim_combo, TargSpTyNum_lim_combo, TargTeff_lim_combo, DerivedMasses_lim_combo, DerivedMasses_minlim_combo, DerivedMasses_maxlim_combo, DerivedMasses_lim_b98_combo, DerivedMasses_minlim_b98_combo, DerivedMasses_maxlim_b98_combo, DerivedLuminosities_lim_combo, DerivedLuminosities_minlim_combo, DerivedLuminosities_maxlim_combo, DerivedLuminosities_lim_b98_combo, DerivedLuminosities_minlim_b98_combo, DerivedLuminosities_maxlim_b98_combo)


    np.savetxt(targets_to_convert.split('.')[0] + '_DerivedMass_Results_B98+B15+MESA_2Myr.txt', to_txt, delimiter=';', fmt='%s', header = 'Name; SpTy; SpTyNum; Teff; Mstar_b15+mesa; Mstar_ll; Mstar_ul; Mstar_b98+mesa; Mstar_b98_ll; Mstar_b98_ul; LogLstar_b15+mesa; LogLstar_ll; LogLstar_ul; LogLstar_b98+mesa; LogLstar_b98_ll; LogLstar_b98_ul')








'''
# Assuming we want to look at the 1 Myr models:
'''

if ages == "1":
    # Now write the new Baraffe+1998/2015 1-myr models and MESA 1-myr models to a new output file for central object mass:


    ##########
    # Unpack the target file, which is an output from Luhman2003-KenyonHartmann_SpTy_Teff_relation.py
    # That output file should be in the following form:
    # Target SpTy NumericalSpTy Teff
    ##########

    TargetName, TargSpTy, TargSpTyNum, TargTeff, TargTeff_earlylim, TargTeff_latelim = np.genfromtxt(targets_to_convert, delimiter=';', dtype='str').T

    TargetName_lim_combo, TargSpTy_lim_combo, TargSpTyNum_lim_combo, TargTeff_lim_combo, DerivedMasses_lim_combo, DerivedMasses_minlim_combo, DerivedMasses_maxlim_combo = [],[],[],[],[],[],[]

    # Make the Teff and its min and max values into an array of floats
    TargTeff = [float(x) for x in TargTeff]
    TargTeff = np.array(TargTeff)

    TargTeff_max = [float(x) for x in TargTeff_earlylim] # assuming earlier spectral types have larger Teffs...
    TargTeff_max = np.array(TargTeff_max)

    TargTeff_min = [float(x) for x in TargTeff_latelim] # assuming later spectral types have smaller Teffs...
    TargTeff_min = np.array(TargTeff_min)

    '''# First Baraffe 2015 only!'''
    #Apply interpolation for all values within the Baraffe models 
    TargTeff_lim = TargTeff[np.where(TargTeff < 4211.)]

    TargTeff_maxlim = TargTeff_max[np.where(TargTeff < 4211.)]

    TargTeff_minlim = TargTeff_min[np.where(TargTeff < 4211.)]

    # And match the target name, etc, arrays to correspond to these values:
    TargetName_lim = TargetName[np.where(TargTeff < 4211.)]
    TargSpTy_lim = TargSpTy[np.where(TargTeff < 4211.)]
    TargSpTyNum_lim = TargSpTyNum[np.where(TargTeff < 4211.)]

    DerivedMasses_lim = spl_b15_1myr(TargTeff_lim) # for the values
    DerivedMasses_maxlim = spl_b15_1myr(TargTeff_maxlim) # for the upper limit (earlier spty - hotter)
    DerivedMasses_minlim = spl_b15_1myr(TargTeff_minlim) # for the lower limit (later spty - cooler)

    DerivedLuminosities_lim = spl_b15_lum(TargTeff_lim)
    DerivedLuminosities_maxlim = spl_b15_lum(TargTeff_maxlim)
    DerivedLuminosities_minlim = spl_b15_lum(TargTeff_minlim)


    '''# Then Baraffe 1998 only!'''
    #Apply interpolation for all values within the Baraffe models (older)

    DerivedMasses_lim_b98 = spl_b98_1myr(TargTeff_lim) # for the values 
    DerivedMasses_maxlim_b98 = spl_b98_1myr(TargTeff_maxlim) # for the upper limit (earlier spty - hotter)
    DerivedMasses_minlim_b98 = spl_b98_1myr(TargTeff_minlim) # for the lower limit (later spty - cooler)

    DerivedLuminosities_lim_b98 = spl_b98_lum(TargTeff_lim)
    DerivedLuminosities_maxlim_b98 = spl_b98_lum(TargTeff_maxlim)
    DerivedLuminosities_minlim_b98 = spl_b98_lum(TargTeff_minlim)


    '''# Then MESA only!'''

    #Apply interpolation for higher mass targets:
    TargTeff_lim_mesa = TargTeff[np.where(TargTeff > 4211.)]

    TargTeff_maxlim_mesa = TargTeff_max[np.where(TargTeff > 4211.)]

    TargTeff_minlim_mesa = TargTeff_min[np.where(TargTeff > 4211.)]

    # And match the target name, etc, arrays to correspond to these values:
    TargetName_lim_mesa = TargetName[np.where(TargTeff > 4211.)]
    TargSpTy_lim_mesa = TargSpTy[np.where(TargTeff > 4211.)]
    TargSpTyNum_lim_mesa = TargSpTyNum[np.where(TargTeff > 4211.)]

    DerivedMasses_lim_mesa = spl_mesa(TargTeff_lim_mesa) # for the values
    DerivedMasses_maxlim_mesa = spl_mesa(TargTeff_maxlim_mesa) # for the upper limit (earlier spty - hotter)
    DerivedMasses_minlim_mesa = spl_mesa(TargTeff_minlim_mesa) # for the lower limit (later spty - cooler)

    DerivedLuminosities_lim_mesa = spl_mesa_lum(TargTeff_lim_mesa) # for the values
    DerivedLuminosities_maxlim_mesa = spl_mesa_lum(TargTeff_maxlim_mesa) # for the upper limit (earlier spty - hotter)
    DerivedLuminosities_minlim_mesa = spl_mesa_lum(TargTeff_minlim_mesa) # for the lower limit (later spty - cooler)


    '''# combine everything to write to output file'''

    TargetName_lim_combo = np.concatenate([TargetName_lim, TargetName_lim_mesa])
    TargSpTy_lim_combo = np.concatenate([TargSpTy_lim, TargSpTy_lim_mesa])
    TargSpTyNum_lim_combo = np.concatenate([TargSpTyNum_lim,TargSpTyNum_lim_mesa])
    TargTeff_lim_combo = np.concatenate([TargTeff_lim, TargTeff_lim_mesa])
    DerivedMasses_lim_combo = np.concatenate([DerivedMasses_lim, DerivedMasses_lim_mesa])
    DerivedMasses_minlim_combo = np.concatenate([DerivedMasses_minlim, DerivedMasses_minlim_mesa])
    DerivedMasses_maxlim_combo = np.concatenate([DerivedMasses_maxlim, DerivedMasses_maxlim_mesa])

    DerivedMasses_lim_b98_combo = np.concatenate([DerivedMasses_lim_b98, DerivedMasses_lim_mesa])
    DerivedMasses_minlim_b98_combo = np.concatenate([DerivedMasses_minlim_b98, DerivedMasses_minlim_mesa])
    DerivedMasses_maxlim_b98_combo = np.concatenate([DerivedMasses_maxlim_b98, DerivedMasses_maxlim_mesa])

    DerivedLuminosities_lim_combo = np.concatenate([DerivedLuminosities_lim, DerivedLuminosities_lim_mesa])
    DerivedLuminosities_minlim_combo = np.concatenate([DerivedLuminosities_minlim, DerivedLuminosities_minlim_mesa])
    DerivedLuminosities_maxlim_combo = np.concatenate([DerivedLuminosities_maxlim, DerivedLuminosities_maxlim_mesa])

    DerivedLuminosities_lim_b98_combo = np.concatenate([DerivedLuminosities_lim_b98, DerivedLuminosities_lim_mesa])
    DerivedLuminosities_minlim_b98_combo = np.concatenate([DerivedLuminosities_minlim_b98, DerivedLuminosities_minlim_mesa])
    DerivedLuminosities_maxlim_b98_combo = np.concatenate([DerivedLuminosities_maxlim_b98, DerivedLuminosities_maxlim_mesa])


    # and write the combination to file:

    to_txt = zip(TargetName_lim_combo, TargSpTy_lim_combo, TargSpTyNum_lim_combo, TargTeff_lim_combo, DerivedMasses_lim_combo, DerivedMasses_minlim_combo, DerivedMasses_maxlim_combo, DerivedMasses_lim_b98_combo, DerivedMasses_minlim_b98_combo, DerivedMasses_maxlim_b98_combo, DerivedLuminosities_lim_combo, DerivedLuminosities_minlim_combo, DerivedLuminosities_maxlim_combo, DerivedLuminosities_lim_b98_combo, DerivedLuminosities_minlim_b98_combo, DerivedLuminosities_maxlim_b98_combo)


    np.savetxt(targets_to_convert.split('.')[0] + '_DerivedMass_Results_B98+B15+MESA_%.0fMyr.txt' % float(ages), to_txt, delimiter=';', fmt='%s', header = 'Name; SpTy; SpTyNum; Teff; Mstar_b15+mesa; Mstar_ll; Mstar_ul; Mstar_b98+mesa; Mstar_b98_ll; Mstar_b98_ul; LogLstar_b15+mesa; LogLstar_ll; LogLstar_ul; LogLstar_b98+mesa; LogLstar_b98_ll; LogLstar_b98_ul')
