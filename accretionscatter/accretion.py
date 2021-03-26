''''
    
    AUTHOR:
        Joe Palmo

'''





################### IMPORT STATEMENTS #########################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import seaborn as sb
from astropy import constants as const
import random
import astropy.constants as const
import math
from tqdm import tqdm
import extinction as ex
import pdb
import glob
import scipy.optimize as optimization
import scipy.interpolate as sinterp
from sklearn.neighbors import KernelDensity
import pickle

################### Mdot Calculations #########################

def unlog(x):
    return 10**x    
    
def flux_to_luminosity(flux, distance):
    '''
    This function will turn a flux estimate (erg/s/cm^2) into a luminosity using distance
    '''
    lum = flux * (4*np.pi*(distance*const.pc.to('cm').value)**2)
    return lum
    

def luminosity_to_flux(lum, dist):
    '''
    This function will turn a flux estimate (erg/s/cm^2) into a luminosity using distance
    '''
    flux = lum / (4*np.pi*(dist*const.pc.to('cm').value)**2)
    return flux
    
    
def Lacc_to_Mdot(Lacc, mass, radius, Rin=5):
    '''
    This function will turn an accretion luminosity into a mass accretion rate estimate using the widely
    accepted relationship. The inputs of the function are:
    
    Lacc (float) : Accretion luminosity [Lsun]
    mass (float) : mass of object [Msun]
    radius (float) : radius of object [Rsun]
    
    Optional:
    Rin (float) : magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    Mdot (float) : mass accretion rate [Msun/yr]
    '''
    Mdot = (1-(radius/Rin))**-1*(radius*const.R_sun.value*Lacc*const.L_sun.value)/(const.G.value*mass*const.M_sun.value) * (365*24*3600) / (const.M_sun.value)


    return Mdot
    
    
def Mdot_to_Lacc(Mdot, mass, radius, Rin=5):
    '''
    This function will turn a mass accretion rate estimate into an accretion luminosity using the widely
    accepted relationship.
    
    Inputs:
    Mdot (float) : mass accretion rate [Msun/yr]
    mass (float) : mass of object [Msun]
    radius (float) : radius of object [Rsun]
    
    Optional:
    Rin (float) : magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    Lacc (float) : Accretion luminosity [Lsun]
    '''
    Lacc = Mdot*(1-(radius/Rin))*(const.G.value*mass*const.M_sun.value)*(const.M_sun.value)/(radius*const.R_sun.value)/(365*24*3600)/const.L_sun.value


    return Lacc
    
    
def UVExcess_to_Mdot(UVexcess, bc, dist, mass, radius, Av, Rin=5):
    '''
    This function will transform a UV Excess flux value into a mass accretion rate estimate by following 
    the process described in Herczeg 2008. 
    
    Inputs:
    UVexcess - UV continuum excess flux [erg/(s*cm^2)]
    bc - bolometric correction
    dist - distance to object [pc]
    mass - mass of object [Msun]
    radius - radius of object [Rsun]
    
    Optional:
    Rin - magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    Mdot - mass accretion rate [Msun/yr]
    '''
    #Extinction correction of flux
    deredUV = ex.remove(Av, UVexcess)
    
    #use the bolometric correction factor to scale the UV excess flux to total accretion flux
    total_excess = deredUV*bc
    
    #accretion flux to accretion luminosity
    Lacc = flux_to_luminosity(total_excess, dist)
    
    #convert accretion luminosity to solar luminosity
    Lacc = Lacc / const.L_sun.to('erg/s').value
    
    #accretion luminosity to Mdot
    Mdot = Lacc_to_Mdot(Lacc, mass, radius, Rin)
    
    return Mdot
    
    
def Mdot_to_UVExcess(Mdot, bc, dist, mass, radius, Av, Rin=5):
    '''
    This function will transform a mass accretion rate estimate to a UV Excess flux value by following 
    the process described in Herczeg 2008. 
    
    Inputs:
    Mdot - mass accretion rate [Msun/yr]
    bc - bolometric correction
    dist - distance to object [pc]
    mass - mass of object [Msun]
    radius - radius of object [Rsun]
    
    Optional:
    Rin - magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    UVexcess - UV continuum excess flux [erg/(s*cm^2)]
    '''
    #Mdot to accretion luminosity
    Lacc = Mdot_to_Lacc(Mdot, mass, radius, Rin)
    
    #convert accretion luminosity to erg/s
    Lacc = Lacc * const.L_sun.to('erg/s').value
    
    #accretion luminosity to accretion flux
    total_excess = luminosity_to_flux(Lacc, dist)
    
    #use the bolometric correction factor to scale the total accretion flux to UV excess flux
    UVexcess = total_excess / bc
    
    #Extinction correction of flux
    redUV = ex.apply(Av, UVexcess)
    
    return redUV
    

def UbandExcess_to_Mdot(Uexcess, dist, mass, radius, Av, Rin=5, unc=False):
    '''
    This function will transform a U band flux value into a mass accretion rate estimate by following 
    the process described in Robinson 2019. 
    
    Inputs:
    Uexcess - U band continuum excess flux [erg/(s*cm^2)]
    dist - distance to object [pc]
    mass - mass of object [Msun]
    radius - radius of object [Rsun]
    
    Optional:
    Rin - magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    Mdot - mass accretion rate [Msun/yr]
    '''
    #Extinction correction of flux
    deredU = ex.remove(Av, Uexcess) 
    
    #U-band flux to Lu
    Lu = flux_to_luminosity(deredU, dist)
    
    #convert Lu to solar luminosity
    Lu = Lu / const.L_sun.to('erg/s').value
    
    #Lu to Lacc using Robinson paper -- natural logarithms
    #uncertainties 0.03 for each constant
    if unc == False:
        logLacc = 0.93 *np.log(Lu) + 0.5 
    else:
        logLacc = (0.93+np.random.normal(0.03)) *np.log(Lu) + (0.5+np.random.normal(0.03)) 
    Lacc = np.exp(logLacc)
    
    #accretion luminosity to Mdot
    Mdot = Lacc_to_Mdot(Lacc, mass, radius, Rin)
    
    return Mdot
    


def lineflux_to_Mdot(flux, dist, mass, radius, Av, Rin=5, line=None, A=None, B=None):
    '''
    This function will turn a line flux into a mass accretion rate estimate using the Lacc-Lline fits derived
    by Alcala et al 2017. 
    
    Inputs:
    flux (float) : line flux [erg/(s*cm^2)]?
    dist (float) : distance to object [pc]
    mass (float) : mass of object [Msun]
    radius (float) : radius of object [Rsun]
    
    Optional:
    Rin (float) : magnetospheric radius (default 5 [Rsun])
    line (str) : type of line. for now acceptable inputs are H-alpha: 'Ha', Pa-beta: 'Pab', and Br-gamma: 'Brg'
    A (float) : If you want to input the parameters for your own line flux vs Lacc relationship
    B (float) : If you want to input the parameters for your own line flux vs Lacc relationship
    
    Outputs:
    Mdot (float) : mass accretion rate [Msun/yr]
    '''
    
    #a & b values pulled directly from the paper
    
    if line == None:
        a = A
        b = B
    elif line == 'Ha':
        a=1.13 #+/- 0.05
        b=1.74  #+/- 0.19
    elif line == 'Pab':
        a=1.06 #+/- 0.07
        b=2.76  #+/- 0.34
    elif line == 'Brg':
        a=1.19 #+/- 0.10
        b=4.02  #+/- 0.51
    else:
        print('Line not found.')
        return
    
    #extinction correction
    deredflux = ex.remove(Av, flux)
    
    #find Lline in erg/s
    Lline = deredflux * (4*np.pi*(dist*const.pc.to('cm').value)**2)

    #convert to solar luminosity
    Lline = Lline / const.L_sun.to('erg/s').value
    
    #Find Lacc using Alcala relationships
    logLacc = a*np.log10(Lline)+b
    #solar luminosity
    Lacc = unlog(logLacc)
    
    Mdot = Lacc_to_Mdot(Lacc, mass, radius, Rin=Rin)
    
    return Mdot
    
    
    
def Mdot_to_lineflux(Mdot, dist, mass, radius, Av, Rin=5, line=None, A=None, B=None):
    '''
    This function will turn a mass accretion rate estimate into a line flux using the Lacc-Lline fits derived
    by Alcala et al 2017. 
    
    Inputs:
    Mdot (float) : mass accretion rate [Msun/yr]
    dist (float) : distance to object [pc]
    mass (float) : mass of object [Msun]
    radius (float) : radius of object [Rsun]
    
    Optional:
    Rin (float) : magnetospheric radius (default 5 [Rsun])
    line (str) : type of line. for now acceptable inputs are H-alpha: 'Ha', Pa-beta: 'Pab', and Br-gamma: 'Brg'
    A (float) : If you want to input the parameters for your own line flux vs Lacc relationship
    B (float) : If you want to input the parameters for your own line flux vs Lacc relationship
    
    Outputs:
    flux (float) : line flux [erg/(s*cm^2)]
    '''
    
    #a & b values pulled directly from the paper
    
    if line == None:
        a = A
        b = B
    elif line == 'Ha':
        a=1.13 #+/- 0.05
        b=1.74  #+/- 0.19
    elif line == 'Pab':
        a=1.06 #+/- 0.07
        b=2.76  #+/- 0.34
    elif line == 'Brg':
        a=1.19 #+/- 0.10
        b=4.02  #+/- 0.51
    else:
        print('Line not found.')
        return
    
    # Mdot to Lacc
    Lacc = Mdot_to_Lacc(Mdot, mass, radius, Rin=Rin)
    
    #log(Lacc)
    logLacc = np.log10(Lacc)
    
    #Find log(Lline) (in solar luminosity) using Alcala relationships
    logLline = (logLacc-b)/a
    
    #unlog, convert Lline to erg/s
    Lline = unlog(logLline)
    Lline = Lline * const.L_sun.to('erg/s').value
    
    #convert Lline to flux
    flux = Lline / ((4*np.pi*(dist*const.pc.to('cm').value)**2))
    
    #extinction correction
    redflux = ex.apply(Av, flux)
    
    return redflux
    

def empiricalMdot(mass, scalefactor=1.79941029, intercept=7.99351629):
    '''
    This function will empirically estimate a mass accretion rate value using the 
    
    Inputs:
    mass (float) : mass [Msun]
    
    Optional:
    scalefactor (float) : slope of empirical power law fit between mass and mass accretion rate in log space 
                            (current value is derived from AccDB)
    intercept (str) :  intercept of empirical power law fit between mass and mass accretion rate in log space 
                            (current value is derived from AccDB)
    
    Outputs:
    Mdot (float) : empirical mass accretion rate [Msun/yr]
    '''
    #Herczeg relation: scalefactor=1.87, intercept=7.73
    Mdot = (mass**scalefactor) / (10**intercept)
    
    return (Mdot)


#Load in age to intercept function
with open('ageinterceptfunc.pickle', 'rb') as f:
    age_to_intercept_func = pickle.load(f)
    

def getIntercept(age):
    '''
    This function will modify the intercept of the empirical relationship based on the age of the object
    
    Inputs:
    age (float) : age [Myr] 
    
    Outputs:
    intercept (float) : intercept of empirical power law fit between mass and mass accretion rate in log space based on the age
                        of the object
    '''
    return age_to_intercept_func(age)
    
def age_intercept_exponential(age, popt=np.array([ 1.48749094,  0.50602283, -8.61059045]), 
    pcov=np.array([[ 9.66619837e-04,  3.79765397e-04, -2.74220068e-06],
                   [ 3.79765397e-04,  3.08006780e-04,  1.19773965e-04],
                   [-2.74220068e-06,  1.19773965e-04,  1.16096209e-04]])):
    '''
    This function will modify the intercept of the empirical relationship based on the age of the object, and take into account
    the uncertainties of the exponential fit between age and intercept
    
    Inputs:
    age (float) : age [Myr] 
    
    Optional:
    popt (array-like) : empirically derived values for a, b, and c in an exponential fit between age and intercept
    pcov (array-like) :  covariance matrix of the exponential fit
    
    Outputs:
    intercept (float) : intercept of empirical power law fit between mass and mass accretion rate in log space based on the age
                        of the object
    '''
    
    a_random, b_random, c_random = np.random.multivariate_normal(popt, pcov)
    
    return a_random * np.exp(-b_random * age) + c_random



################### Object Parameter Estimation #########################


#Need to have these models handy to derive stellar parameters
models = glob.glob('StellarParams/Baraffe*txt')
mesamodels = glob.glob('StellarParams/MESA_*.txt')


def Num_to_SpTy(spectypenum):
    '''
    This function takes a spectral type in the form - 'letternumber', i.e. 'M5'.
    It then translates that spec type into a spectral type number identifier used for interpolation.
    '''
    spty_dict = {0 : 'B',
                 1 : 'A',
                 2 : 'F',
                 3 : 'G',
                 4 : 'K',
                 5 : 'M'}
    
    if int(spectypenum/10) > 5:
        letter = 'M'
        number = int(spectypenum%10)+10
    else:
        letter = spty_dict[int(spectypenum/10)]
        number = int(spectypenum%10)
    
    spty = letter+str(number)
    
    return spty

def to_SpTyNum(spectype):
    '''
    This function takes a spectral type in the form - 'letternumber', i.e. 'M5'.
    It then translates that spec type into a spectral type number identifier used for interpolation.
    '''
    spty_dict = {'B' : 0,
                 'A' : 1,
                 'F' : 2,
                 'G' : 3,
                 'K' : 4,
                 'M' : 5}
    
    letter = spectype[0]
    number = spectype[1:]
    
    sptynum = spty_dict[letter]*10 + int(number)
    
    return sptynum
    
    
def SpTy_to_Teff(spectype):
    '''
    This function will take a numerical spectral type identifier, and using interpolation from the tables in
    Herczeg and Hillenbrand 2014 it will calculate an effective temperature.
    '''
    SpTyNum, Teff, SpTy = np.genfromtxt('StellarParams/HerczegHillenbrand_SpTyTeff_Numericized.txt', skip_header=1, dtype='str').T
    SpTyNum, Teff = [float(x) for x in SpTyNum], [float(y) for y in Teff]
    
    spl = sinterp.UnivariateSpline(SpTyNum, Teff)
    
    teff = spl(spectype)
    
    return teff
    
    
def Teff_to_params(Teff, age):
    '''
    This function will take an effective temperature and an age, and by interpolating the Baraffe 2015 models 
    it will return a mass and radius estimate for the object.
    
    Inputs:
    Teff (float) - effective temperature in Kelvin
    age (float) - age of object in Myr
    
    Outputs:
    m (float) - interpolated mass of object in solar masses
    r (floar) - interpolated radius of object in solar radii
    
    '''
    #Find the most accurate model given age
    ages = [1, 2, 3, 8, 10]
    age_diff = []
    for a in ages:
        age_diff.append(np.abs(age-a))
    closest_age = ages[age_diff.index(min(age_diff))]
    
    model = None
    for m in models:
        if (str(closest_age)+'Myr') in m:
            model = m
            
    mass, teff, loglum, logg, radius = np.loadtxt(model, skiprows=1).T
    
    #Polynomial Fit -- doesn't work outside of model grid boundaries
    #p_mass = np.polyfit(teff, mass, 10)
    #f_mass = np.poly1d(p_mass)
    f_mass = sinterp.interp1d(teff, mass, fill_value='extrapolate', kind='linear')
    
    #Polynomial Fit -- doesn't work outside of model grid boundaries
    #p_radius = np.polyfit(teff, radius, 10)
    #fradius = np.poly1d(p_radius)
    #newteff = np.linspace(min(teff), max(teff), 10000)
    f_radius = sinterp.interp1d(teff, radius, fill_value='extrapolate', kind='linear')
    
    m = f_mass(Teff)
    r = f_radius(Teff)
    
    return m, r
    
    
def mass_to_Teff(massobj, age):
    '''
    This function will take a mass and an age, and by interpolating the Baraffe 2015 models 
    it will return a temperature estimate for the object.
    '''
    #Find the most accurate model given age
    ages = [1, 2, 3, 8, 10]
    age_diff = []
    for a in ages:
        age_diff.append(np.abs(age-a))
    closest_age = ages[age_diff.index(min(age_diff))]
    
    model = None
    for m in models:
        if (str(closest_age)+'Myr') in m:
            model = m
            
    mass, teff, loglum, logg, radius = np.loadtxt(model, skiprows=1).T
    
    #Polynomial Fit -- doesn't work outside of model grid boundaries
    #p_teff = np.polyfit(mass, teff, 10)
    #f_teff = np.poly1d(p_teff)
    f_teff = sinterp.interp1d(mass, teff, fill_value='extrapolate', kind='linear')
    
    t = f_teff(massobj)
    
    return t
    
    
def Teff_to_SpTy(teff):
    '''
    This function will take an effective temperature, and using interpolation from the tables in
    Herczeg and Hillenbrand 2014 it will calculate a numerical spectral type identifier.
    '''
    SpTyNum, Teff, SpTy = np.genfromtxt('StellarParams/HerczegHillenbrand_SpTyTeff_Numericized.txt', skip_header=1, dtype='str').T
    SpTyNum, Teff = [float(x) for x in SpTyNum], [float(y) for y in Teff]
    SpTyNum.reverse()
    Teff.reverse()
    
    spl = sinterp.UnivariateSpline(Teff, SpTyNum)
    
    spty = spl(teff)
    
    return spty
    

'''
########### 2D Interpolation for Object Parameters ############
###### code shown below, functions saved as pickle files ######
##### This didn't end up being used, the 2D interpolation failed in certain ranges. #####

from astropy.io import ascii, fits
#Load in the isochrones
iso = ascii.read('BHAC15.dat')
#Interpolate between the models. This may take a minute.
    
print('Creating interpolation function for T')
MtoTfunc = sinterp.interp2d(iso['M'], iso['logT'], iso['Teff'], kind = 'linear')
    
print('Creating interpolation function for R')
MtoRfunc = sinterp.interp2d(iso['M'], iso['logT'], iso['R'], kind = 'linear')   
    
print('Creating interpolation function for M')
TtoMfunc = sinterp.interp2d(iso['Teff'], iso['logT'], iso['M'], kind = 'linear')
    
print('Creating interpolation function for R')
TtoRfunc = sinterp.interp2d(iso['Teff'], iso['logT'], iso['R'], kind = 'linear')

with open('MtoTfunc.pickle', 'wb') as f:
    pickle.dump(MtoTfunc, f)

with open('MtoRfunc.pickle', 'wb') as f:
    pickle.dump(MtoRfunc, f)

with open('TtoMfunc.pickle', 'wb') as f:
    pickle.dump(TtoMfunc, f)
    
with open('TtoRfunc.pickle', 'wb') as f:
    pickle.dump(TtoRfunc, f)
    
######## load in MtoTfunc.pickle, MtoRfunc.pickle, TtoMfunc.pickle, TtoRfunc.pickle to use #########
    
with open('MtoTfunc.pickle', 'rb') as f:
    T = pickle.load(f)
    
with open('MtoRfunc.pickle', 'rb') as f:
    R = pickle.load(f)
    
with open('TtoMfunc.pickle', 'rb') as f:
    TtoMfunc = pickle.load(f)
    
with open('TtoRfunc.pickle', 'rb') as f:
    TtoRfunc = pickle.load(f)
'''    
    

################### Uncertainty Distributions #########################


def gaussian(x, mu, sigma):
    '''
    Builds a gaussian distribution for a parameter using the formula for a gaussian.
    
    Inputs:
    x (array-like) - the x-values used to generate f
    mu (float) - mean value of the parameter
    sigma (float) - uncertainty of the parameter
    
    Outputs:
    g (array-like) - a gaussian curve
    '''
    g = (1/(sigma*np.sqrt(2*np.pi)))*np.exp((-1/2)*(x-mu)**2/sigma)
    return g
    
def cdf(g, dx):
    cdf = np.cumsum(g) * dx
    '''
    Sums all of the values of a gaussian probability distribution function (cdf) to compute its cumulative distribution function 
    (cdf)
    
    Inputs:
    g (array-like) - a gaussian curve
    
    Outputs:
    dx (float) - the step size of the numerical integration
    '''
    return cdf
    

def inverse_transform_sampling(x, f, nsamples):
    '''
    This function performs the inverse transform sampling method to sample from a probability distribution 
    function (PDF) using uniform samples between 0 and 1. It does this by calculating the cumulative distribution 
    function (CDF) of the normalized PDF, inverting it, then sampling along the probability axis from 0 to 1.
    
    Inputs:
    x (array-like) - the x-values used to generate f
    f (array-like) - the PDF
    nsamples (int) - the number of inverse transform sampled values needed, i.e. the size of the output
    
    Outputs:
    samples (array-like) - an inverse transform sampled distribution which resembles the PDF
    '''
    
    dx = x[1]-x[0]
    normalized = f / np.trapz(f, x=x, dx=dx)
    c = cdf(normalized, dx)
    inv = sinterp.interp1d(c, x, fill_value='extrapolate')
    r = np.random.rand(nsamples)
    samples = inv(r)
    samples = samples[np.isfinite(samples)]
    
    return np.array(samples)
    
def uncdist(x, mu, sigma, nsamples):
    '''
    Builds a truncated gaussian uncertainty distribution for a parameter using inverse transform sampling.
    
    Inputs:
    x (array-like) - the x-values used to generate f
    mu (float) - mean value of the parameter
    sigma (float) - uncertainty of the parameter
    nsamples (int) - the number of inverse transform sampled values needed, i.e. the size of the output
    
    Outputs:
    samples (array-like) - an inverse transform sampled distribution which resembles the PDF
    '''
    
    f = gaussian(x, mu, sigma)
    samples = inverse_transform_sampling(x, f, nsamples)
    if nsamples == 1:
        if samples.size != 0:
            return samples[0]
        else:
            return np.nan
    else:
        return samples


def pdf_fit(x, data, kernel='gaussian'):
    '''
    This function fits a probability density function to input data (data) using a non-parametric
    kernel density estimation.
    
    Inputs:
    x (array-like) - range of x-values used to generate the pdf, must cover the range of the data
    data (array-like) - must be in order from least to greatest, and contain no NaN values
    
    Optional:
    kernel (str) - the kernel to use
                {‘gaussian’, ‘tophat’, ‘epanechnikov’, ‘exponential’, ‘linear’, ‘cosine’}, default=’gaussian’
    
    Outputs:
    pdf (array-like) - a fitted probability density function for the histogram of data
    '''
    
    model = KernelDensity(bandwidth=2, kernel=kernel)
    d = data.reshape((len(data), 1))
    model.fit(d)
    
    pdf = model.score_samples(x.reshape((len(x), 1)))
    pdf = np.exp(pdf)
    
    return pdf


    
    
    
################### IMFs #########################

def MillerScalo1979(m, size):
    '''
    This function builds a distribution of masses based on the initial mass function (IMF) derived in Miller-Scalo et al. 1979.
    
    Inputs:
    m (array-like) - range of mass values used to generate the distribution of masses
    size (int) - the number of values that the user wants in the distribution of masses
    
    Outputs:
    inv (array-like) - a distribution of masses of with size values that follows the Miller-Scalo IMF
    '''
    x = np.copy(m)
    a1 = 0   # m < 1
    a2 = 2.3  # m > 1
    m[m<1] = (m[m<1]**(-a1))
    m[(m>1)] = (m[(m>1)]**(-a2))
    inv = inverse_transform_sampling(x, m, size)
    return inv
    
def Romano2005(m, size):
    '''
    This function builds a distribution of masses based on the updated Chabrier initial mass function (IMF) derived in
    Romano et al. 2005.
    
    Inputs:
    m (array-like) - range of mass values used to generate the distribution of masses
    size (int) - the number of values that the user wants in the distribution of masses
    
    Outputs:
    inv (array-like) - a distribution of masses of with size values that follows the Romano/Chabrier IMF
    '''
    x = np.copy(m)
    
    # Chabrier IMF info from Romano+ 2005
    mc = 0.079 #Msun
    sigma = 0.69
    fa = 0.85#1.06
    fb = 0.24#0.30
    c = 1.3#1.7
    m[m<=1] = fa * np.exp(-1 * ((np.log10(m[m<=1]) - np.log10(mc))**2. / (2*sigma**2))) 
    m[m>1] = fb * m[m>1]**(-c)
    inv = inverse_transform_sampling(x, m, size)
    return inv
    
    
    
    