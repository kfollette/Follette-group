##FUNCTIONS:

#1. Line flux to line luminosity
#2. Line luminosity to accretion luminosity
#3. Accretion luminosity to accretion rate
#4. Get accretion rate starting from line flux

#This file additionally imports the scaling relation coefficients.

import astropy.units as u
from astropy.constants import G, M_jup, R_jup, M_earth, R_earth, L_sun, M_sun, R_sun
from numpy import *
import pandas as pd

rel = pd.read_csv('scalingrels_a17.csv')
rel.set_index('Tracer', inplace=True)

def line_lum(line_flux, dist):
    """
    Calculate line luminosity given line flux and distance
    assuming line flux is extinction corrected.
    """    
    line_lum = 4 * pi * (dist*u.pc)**2 * line_flux * u.erg / (u.s * (u.cm)**2)
    line_lum = line_lum.decompose().to(u.W)
    return line_lum/u.W

def accr_lum(L_line, tracer, L_line_err = 0*u.W):
    """
    Translate a line luminosity to accretion luminosity using empirical
    relationships from Alcala et al. 2017.
    
    Included tracers are:
    'Ha'
    'PaB'
    'BrG'
    """
        
    a, a_err, b, b_err = rel['a'][tracer],rel['a_err'][tracer],rel['b'][tracer],rel['b_err'][tracer]
    
    log_L_acc = b + a * log10(L_line*u.W/L_sun)
    
    L_acc = 10**log_L_acc*L_sun/u.W
    
    #error propagation
    
    #c_err = (L_line_err)/(log(10) * L_line)
    #ac_err = a * log10(L_line/L_sun) * ((a_err/a)**2 + (c_err/log10(L_line/L_sun))**2)**0.5
    #log_L_acc_err = (b_err**2 + ac_err**2)**0.5
    #L_acc_err = L_acc * log(10) * log_L_acc_err

    return L_acc

def acc_rate(L_acc, R, M):
    """
    Translate an accretion luminosity and planet mass/radius to accretion rate in Solar masses per year.
    """
    mdot = 1.25*L_acc*u.W*R*u.R_sun/(G*M*u.M_sun)
    mdot = mdot.decompose().to(u.M_sun/u.yr)
    return(mdot/(u.M_sun/u.yr))

def get_rate(line_flux, d, t, R, M):
    """
    Turn a line flux into an accretion rate, 
    given the distance, tracer, object radius, and object mass.
    
    Line flux should be in erg/s/cm^2; distance in pc;
    tracer either Ha, PaB, or BrG; radius in solar radii,
    and mass in solar masses.
    """
    l_line = line_lum(line_flux,d)
    l_acc = accr_lum(l_line, t)
    mdot = acc_rate(l_acc, R, M)
    
    return mdot