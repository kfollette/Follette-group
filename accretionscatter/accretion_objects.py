''''
    
    AUTHOR:
        Joe Palmo

    Last Edited 3/4/2021
'''

################### IMPORT STATEMENTS #########################

import accretion as a
from importlib import reload
reload(a)

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
from matplotlib.animation import FuncAnimation
import scipy.interpolate as sinterp
from sklearn.neighbors import KernelDensity
import pickle

################# OBJECTS #########################


########################
########################
######################## Accretion Object

class Accretion:
    '''
    Create an object for storing all of the information in an accretion MC error propagation simulation
    '''
    
    # class attributes here -- any variable which has the same value for all class instances
    # maybe the accepted accretion rate relationship can go here
    def __init__(self, mass, distance, age, Av, Rin=5):
        self.mass = mass  #Msun
        self.dist = distance   #pc
        self.age = age    #Myr
        self.Av = Av     #mag
        
        #derived variables
        #self.Teff = MtoTfunc(self.mass, self.age)
        #self.radius = MtoRfunc(self.mass, self.age)
        #1d interpolation method done
        #can't incorporate age errors this way
        self.Teff = a.mass_to_Teff(self.mass, self.age)
        self.radius = a.Teff_to_params(self.Teff, self.age)[1]
        self.SpTy = a.Teff_to_SpTy(self.Teff)
        self.Rin = Rin * self.radius
    
        #ideal mdot, Lacc
        self.ideal_mdot = a.empiricalMdot(self.mass)
        self.ideal_Lacc = a.Mdot_to_Lacc(self.ideal_mdot, self.mass, self.radius, self.Rin)
        
        self.mdot = self.ideal_mdot
        #create mdot array for later
        
        
    def MvMdot(self):
        logaccepted_relation = np.log10(self.ideal_mdot)
        logmass = np.log10(np.ones(len(self.mdot))*self.mass)
        logMdot = np.log10(self.mdot)
        fig, ax = plt.subplots()
        ax.plot(logmass, logaccepted_relation, color='cornflowerblue', label='Empirical Relationship')
        ax.scatter(logmass, logMdot, color='darkseagreen', label='Simulated Points')
        ax.set_xlabel('log(Mass) (M$\odot$)')
        ax.set_ylabel('log(Mass Accretion Rate) (M$\odot$/yr)')
        ax.set_title('Monte Carlo Error Propagation')
        ax.legend()
       
    #def getMdot():
        
        
        
    def UVExcessErrorProp(self, SpTyUnc, distUnc, ageUnc, AvUnc, bc, bcUnc, UVExcessUnc, numMC, Rin=5, RinUnc=0, variability=0, age_scatter=False):
            #propagate errors forward and obtain uncertainty distributions
            if age_scatter == True:
                #self.shifted_ideal_mdot = a.empiricalMdot(self.mass, intercept=-a.age_to_intercept_func(self.age))
                self.shifted_ideal_mdot = a.empiricalMdot(self.mass, intercept=-a.age_intercept_exponential(self.age))
                self.ideal_UVExcess = a.Mdot_to_UVExcess(self.shifted_ideal_mdot, bc, self.dist, self.mass, self.radius, self.Av, self.Rin)
            else:
                self.ideal_UVExcess = a.Mdot_to_UVExcess(self.ideal_mdot, bc, self.dist, self.mass, self.radius, self.Av, self.Rin)
            #add in uncertainties to observable itself
            if numMC == 1:
                self.UVExcessUncDist = np.random.normal(self.ideal_UVExcess, UVExcessUnc*self.ideal_UVExcess, size=numMC)[0]
            else:
                self.UVExcessUncDist = np.random.normal(self.ideal_UVExcess, UVExcessUnc*self.ideal_UVExcess, size=numMC)
            
            #add in uncertainties to variables specific to observable
            x = np.linspace(0, 20000, 100000)
            if bcUnc == 0: 
                if numMC == 1:
                    self.bcUncDist = bc
                else:
                    self.bcUncDist = np.ones(numMC)*bc
            else:
                self.bcUncDist = a.uncdist(x, bc, bcUnc, numMC)
            
            #build uncertainty distributions
            if SpTyUnc == 0:
                if numMC == 1:
                    self.SpTyUncDist = self.SpTy
                else:
                    self.SpTyUncDist = np.ones(numMC)*self.SpTy
            else:
                self.SpTyUncDist = a.uncdist(x, self.SpTy, SpTyUnc, numMC)
                
            if distUnc == 0:
                if numMC == 1:
                    self.distUncDist = self.dist
                else:
                    self.distUncDist = np.ones(numMC)*self.dist
            else:
                self.distUncDist = a.uncdist(x, self.dist, distUnc, numMC)
            
            if ageUnc == 0:
                if numMC == 1:
                    self.ageUncDist = self.age
                else:
                    self.ageUncDist = np.ones(numMC)*self.age
            else:
                self.ageUncDist = a.uncdist(x, self.age, ageUnc, numMC)
                
            if AvUnc == 0:
                if numMC == 1: 
                    self.AvUncDist = self.Av
                else:
                    self.AvUncDist = np.ones(numMC)*self.Av
            else:
                self.AvUncDist = a.uncdist(x, self.Av, AvUnc, numMC)
                
            if RinUnc == 0:
                if numMC == 1: 
                    self.RinUncDist = Rin
                else:
                    self.RinUncDist = np.ones(numMC)*Rin
            else:
                self.RinUncDist = a.uncdist(x, Rin, RinUnc, numMC)

            self.TeffUncDist = a.SpTy_to_Teff(self.SpTyUncDist)
            self.massUncDist = a.Teff_to_params(self.TeffUncDist, (self.ageUncDist))[0]
            self.radiusUncDist = a.Teff_to_params(self.TeffUncDist, (self.ageUncDist))[1]
            self.RinUncDist = self.RinUncDist*self.radiusUncDist
            #self.RinUncDist = Rin*self.radiusUncDist
            
            
            #Generate synthetic Mdot distribution
            #Generate synthetic Lacc distribution
            self.mdot = a.UVExcess_to_Mdot(self.UVExcessUncDist, self.bcUncDist, self.distUncDist, self.massUncDist, self.radiusUncDist, self.AvUncDist, self.RinUncDist)
            if variability != 0:
                log_mdot = np.log10(self.mdot)
                if numMC == 1: 
                    draw = np.random.normal(log_mdot, variability, size=numMC)[0]
                else:
                    draw = np.random.normal(log_mdot, variability, size=numMC)
                
                self.mdot = 10**draw
                
                '''
                #specify lower and upper bounds for the truncated gaussian
                lower, upper = 0,1
                mu, sigma = (self.mdot), (self.mdot*variability)
                variability_distribution = st.truncnorm(
                    (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
                if numMC == 1: 
                    self.mdot = variability_distribution.rvs(numMC)[0]
                else:
                    self.mdot = variability_distribution.rvs(numMC)
                '''
                
            self.Lacc =  a.Mdot_to_Lacc(self.mdot, self.massUncDist, self.radiusUncDist, self.RinUncDist)
            
            #if self.mdot < 0:
                #self.mdot = np.nan
            #if self.Lacc < 0:
                #self.Lacc = np.nan

    # not ready            

    def linefluxErrorProp(self, SpTyUnc, distUnc, ageUnc, AvUnc, linefluxUnc, numMC, Rin=5, RinUnc=0, line=None, A=None, AUnc=None, B=None, BUnc=None, variability=0, age_scatter=False):
            #propagate errors forward and obtain uncertainty distributions
            if line == None:
                A = A
                aUnc = AUnc
                B = B
                bUnc = BUnc
            elif line == 'Ha':
                A=1.13 #+/- 0.05
                aUnc = 0.05
                B=1.74  #+/- 0.19
                bUnc = 0.19
            elif line == 'Pab':
                A=1.06 #+/- 0.07
                aUnc = 0.07
                B=2.76  #+/- 0.34
                bUnc = 0.34
            elif line == 'Brg':
                A=1.19 #+/- 0.10
                aUnc = 0.10
                B=4.02  #+/- 0.51
                bUnc = 0.51
            else:
                print('Line not found.')
                return
            
            if age_scatter == True:
                #self.shifted_ideal_mdot = a.empiricalMdot(self.mass, intercept=-a.age_to_intercept_func(self.age))
                self.shifted_ideal_mdot = a.empiricalMdot(self.mass, intercept=-a.age_intercept_exponential(self.age))
                self.ideal_lineflux = a.Mdot_to_lineflux(self.shifted_ideal_mdot, self.dist, self.mass, self.radius, self.Av, Rin=self.Rin, line=line, A=A, B=B)
            else:
                self.ideal_lineflux = a.Mdot_to_lineflux(self.ideal_mdot, self.dist, self.mass, self.radius, self.Av, Rin=self.Rin, line=line, A=A, B=B)
            
            #add in uncertainties to observable itself
            if numMC == 1:
                self.linefluxUncDist = np.random.normal(self.ideal_lineflux, linefluxUnc*self.ideal_lineflux, size=numMC)[0]
            else:
                self.linefluxUncDist = np.random.normal(self.ideal_lineflux, linefluxUnc*self.ideal_lineflux, size=numMC)
               
            x = np.linspace(0, 1000, 100000)
                    
            #add in uncertainties to variables specific to observable
            #self.AUncDist = a.uncdist(x, A, aUnc, numMC) #A
            #self.BUncDist = a.uncdist(x, B, bUnc, numMC) #B
            #Ignore uncertainties of the empirical fit
            self.AUncDist = A #A
            self.BUncDist = B #B
            
            #build uncertainty distributions
            if SpTyUnc == 0:
                if numMC == 1:
                    self.SpTyUncDist = self.SpTy
                else:
                    self.SpTyUncDist = np.ones(numMC)*self.SpTy
            else:
                self.SpTyUncDist = a.uncdist(x, self.SpTy, SpTyUnc, numMC)
                
            if distUnc == 0:
                if numMC == 1:
                    self.distUncDist = self.dist
                else:
                    self.distUncDist = np.ones(numMC)*self.dist
            else:
                self.distUncDist = a.uncdist(x, self.dist, distUnc, numMC)
            
            if ageUnc == 0:
                if numMC == 1:
                    self.ageUncDist = self.age
                else:
                    self.ageUncDist = np.ones(numMC)*self.age
            else:
                self.ageUncDist = a.uncdist(x, self.age, ageUnc, numMC)
                
            if AvUnc == 0:
                if numMC == 1: 
                    self.AvUncDist = self.Av
                else:
                    self.AvUncDist = np.ones(numMC)*self.Av
            else:
                self.AvUncDist = a.uncdist(x, self.Av, AvUnc, numMC)
                
            if RinUnc == 0:
                if numMC == 1: 
                    self.RinUncDist = Rin
                else:
                    self.RinUncDist = np.ones(numMC)*Rin
            else:
                self.RinUncDist = a.uncdist(x, Rin, RinUnc, numMC)
            
            
            self.TeffUncDist = a.SpTy_to_Teff(self.SpTyUncDist)
            self.massUncDist = a.Teff_to_params(self.TeffUncDist, (self.ageUncDist))[0]
            self.radiusUncDist = a.Teff_to_params(self.TeffUncDist, (self.ageUncDist))[1]
            self.RinUncDist = self.RinUncDist*self.radiusUncDist
            
            #Generate synthetic Mdot distribution
            #Generate synthetic Lacc distribution
            self.mdot = a.lineflux_to_Mdot(self.linefluxUncDist, self.distUncDist, self.massUncDist, self.radiusUncDist, self.AvUncDist, Rin=self.RinUncDist, A=self.AUncDist, B=self.BUncDist)
            if variability != 0:
                #specify lower and upper bounds for the truncated gaussian
                lower, upper = 0,1
                mu, sigma = (self.mdot), (self.mdot*variability)
                variability_distribution = st.truncnorm(
                    (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
                if numMC == 1: 
                    self.mdot = variability_distribution.rvs(numMC)[0]
                else:
                    self.mdot = variability_distribution.rvs(numMC)
            
            self.Lacc =  a.Mdot_to_Lacc(self.mdot, self.massUncDist, self.radiusUncDist, self.RinUncDist)





########################
########################
########################Accretion Distribution Object

##bootstrap This object represents a distribution of accreting objects

class AccretionDistribution:
    '''
    Build a distribution of accreting objects
    '''
    
    
    def __init__(self, observed, masses=None, distances=None, ages=None, Avs=None):
        
        # observed is an iteration of Annie's Database
        # it is used to compare simulated values to the observed
        # first we clean it using the lines below
        observed['A_V'] = observed['A_V'].fillna(0)
        self.observed = observed.dropna(axis=0, subset=['Object Mass M_Solar', 'Distance', 'System Age Myr'])
        
        #masses
        self.masses = masses
        
        #distances
        self.distances = distances
        
        #ages
        self.ages = ages
        
        #extinctions, Av
        self.Avs = Avs
        
        
        #save the observed Mdots, Laccs
        self.observed_mdots = observed['Accretion Rate M_solar yr-1'].tolist()
        self.observed_Laccs = (10**observed['Log Accretion Luminosity (solar)']).tolist()
        
        #def build_distribution(self):
         #   if self.masses == None or self.distances == None or self.ages == None or self.Avs == None:
          #      print('Input Masses, Ages, Distances, and Extinctions for your distribution of accreting objects.')
           # else:
            #    #build distribution
             #   self.accretion_distribution = np.array([Accretion(m, self.distances[i], self.ages[i], self.Avs[i]) for i,m in enumerate(self.masses)])
        
        #build_distribution(self)

    def build_distribution(self):
        #build distribution
        self.accretion_distribution = np.array([Accretion(m, self.distances[i], self.ages[i], self.Avs[i]) for i,m in enumerate(self.masses)])

    
    def bootstrap(self):
        
        #masses
        self.masses = self.observed['Object Mass M_Solar'].tolist()
        
        #distances
        self.distances = self.observed['Distance'].tolist()
        
        #ages
        self.ages = self.observed['System Age Myr'].tolist()
        
        #extinctions, Av
        self.Avs = self.observed['A_V'].tolist()
        
        self.accretion_distribution = np.array([Accretion(m, self.distances[i], self.ages[i], self.Avs[i]) for i,m in enumerate(self.masses)])
      
    ##### Mass Input Strategies #####
    
    #From IMF
    def mass_from_IMF(self, size, imf='Chabrier'):
        if imf == 'Chabrier':
            self.masses = a.Romano2005(np.linspace(0,1.4,10000), size)
        elif imf == 'Miller-Scalo':
            self.masses = a.MillerScalo1979(np.linspace(0,1.4,10000), size)
        else:
            print('IMF not found.')
     
    #From prior. PDF fit using scikit-learn -> inverse transform sampling
    def mass_from_prior(self, size, mass_prior, k='gaussian'):
        #masses
        x_mass = np.linspace(0,5,len(mass_prior))
        self.masses = a.inverse_transform_sampling(x_mass, a.pdf_fit(x_mass, mass_prior, k), size)
        
    ##### Age Input Strategies #####
    
    #From prior. PDF fit using scikit-learn -> inverse transform sampling
    def age_from_prior(self, size, age_prior, k='gaussian'):
        #masses
        x_age = np.linspace(0,5,len(age_prior))
        self.ages = a.inverse_transform_sampling(x_age, a.pdf_fit(x_age, age_prior, k), size)
        
    #Gaussian Distribution centered around a value and with a certain standard deviation estimate.
    def age_from_gaussian(self, size, mu, sigma):
        self.ages = np.random.normal(mu, sigma, size=size)
        
    ##### Distance Input Strategies #####
    
    #From prior. PDF fit using scikit-learn -> inverse transform sampling
    def distance_from_prior(self, size, distance_prior, k='exponential'):
        #masses
        x_distance = np.linspace(0,5,len(distance_prior))
        self.distances = a.inverse_transform_sampling(x_distance, a.pdf_fit(x_distance, distance_prior, k), size)
        
    #Gaussian Distribution centered around a value and with a certain standard deviation estimate.
    def distance_from_gaussian(self, size, mu, sigma):
        self.distances = np.random.normal(mu, sigma, size=size)
      
    ##### Extinction Input Strategies #####
    
    #From prior. PDF fit using scikit-learn -> inverse transform sampling
    def Av_from_prior(self, size, Av_prior, k='exponential'):
        #masses
        x_Av = np.linspace(0,5,len(Av_prior))
        self.Avs = a.inverse_transform_sampling(x_Av, a.pdf_fit(x_Av, Av_prior, k), size)
        self.Avs[self.Avs<0] = 0.00
        
    #Gaussian Distribution centered around a value and with a certain standard deviation estimate.
    def Av_from_gaussian(self, size, mu, sigma):
        self.Avs = np.random.normal(mu, sigma, size=size)
        
    def UVExcessErrorProp(self, SpTyUnc, distUnc, ageUnc, AvUnc, bc, bcUnc, UVExcessUnc, numMC, Rin=5, RinUnc=0,  variability=0, age_scatter=False):
        for acc in self.accretion_distribution:
            acc.UVExcessErrorProp(SpTyUnc, distUnc, ageUnc, AvUnc, bc, bcUnc, UVExcessUnc, numMC, Rin=Rin, RinUnc=RinUnc, variability=variability, age_scatter=age_scatter)
    
    def linefluxErrorProp(self, SpTyUnc, distUnc, ageUnc, AvUnc, linefluxUnc, numMC, Rin=5, RinUnc=0, line=None, A=None, AUnc=None, B=None, BUnc=None, variability=0, age_scatter=False):
        for acc in self.accretion_distribution:
            acc.linefluxErrorProp(SpTyUnc, distUnc, ageUnc, AvUnc, linefluxUnc, numMC, Rin=Rin, RinUnc=RinUnc, line=line, A=A, AUnc=AUnc, B=B, BUnc=BUnc, variability=variability,  age_scatter=age_scatter)
            
    def create_df(self):
        idealmdots = [acc.ideal_mdot for acc in self.accretion_distribution]
        idealLaccs = [acc.ideal_Lacc for acc in self.accretion_distribution]
        masses = [acc.mass for acc in self.accretion_distribution]
        radii = [acc.radius for acc in self.accretion_distribution]
        ages = [acc.age for acc in self.accretion_distribution]
        distances = [acc.dist for acc in self.accretion_distribution]
        Avs = [acc.Av for acc in self.accretion_distribution]
        temperatures = [acc.Teff for acc in self.accretion_distribution]
        #SpTys = [a.Num_to_SpTy(acc.SpTy) for acc in self.accretion_distribution]
        Rins = [acc.Rin for acc in self.accretion_distribution]
        mdots = [acc.mdot for acc in self.accretion_distribution]
        Laccs = [acc.Lacc for acc in self.accretion_distribution]
        
        temp = {'Mass (M$_\odot$)':masses, 'Radius (R$_\odot$)':radii, 'Age (Myr)':ages, 'Distance (pc)':distances, 
                'Teff (K)':temperatures, '''''Spectral Type':SpTys,''' 'Rin (R$_\odot$)':Rins, 'Mdot (M$_\odot$)':mdots, 
                'Lacc (L$_\odot$)':Laccs, '"true" Mdot (M$_\odot$)':idealmdots, '"true" Lacc (L$_\odot$)':idealLaccs}
        df = pd.DataFrame(temp)
        
        return df
        
 
        
        
########################
########################
######################## Plotting Functions   
        
        
def MoneyPlot(observed, simulated):
    logmass = np.log10(simulated['Mass (M$_\\odot$)'])
    logMdot = np.log10(simulated['Mdot (M$_\\odot$)'])
    logaccepted_relation = np.log10(simulated['"true" Mdot (M$_\\odot$)'])
    observed_mass = np.log10(observed['Object Mass M_Solar'])
    observed_mdot = np.log10(observed['Accretion Rate M_solar yr-1'])

    fig = plt.figure(figsize=(16, 12))
    frame1 = fig.add_axes((.1,.3,.8,.6))
    ax = plt.gca()
    ax.scatter(logmass, logMdot, color='darkseagreen', s=130, alpha=0.4, label='Simulated')
    ax.scatter(observed_mass, observed_mdot, color='black', s = 50, marker='x', alpha=0.4, label='Observed')
    ax.plot(logmass, logaccepted_relation, color='tomato', label='Empirical Relationship')
    ax.axvline(x=np.log10(0.012), color='black', linestyle='-.', label='Deuterium Burning Limit')
    ax.axvline(x=np.log10(0.1), color='black', linestyle=':', label='Hydrogen Burning Limit')
    ax.set_xlim(-2.6, 0.6)
    ax.set_xlabel('log(Mass) (M$_\odot$)', size=16)
    ax.set_ylabel('log(Mass Accretion Rate) (M$_\odot$/yr)', size=16)
    ax.text(-2.3, -6.6, 'Planets', size = 15, fontstyle='italic')
    ax.text(-1.64, -6.6, 'Brown Dwarfs', size = 15, fontstyle='italic')
    ax.text(-0.29, -6.6, 'Stars', size = 15, fontstyle='italic')
    ax.set_title('Monte Carlo Error Propagation', size=20, fontweight='heavy')
    ax.legend(loc='lower right', prop={'family':'sans-serif', 'style':'normal', 'size': 16}, frameon=False, shadow=True)

    frame2=fig.add_axes((.1,.1,.8,.2))  
    ax = plt.gca()
    observed_difference = np.log10(observed['Accretion Rate M_solar yr-1']) - np.log10(a.empiricalMdot(observed['Object Mass M_Solar']))
    simulated_difference = np.log10(simulated['Mdot (M$_\\odot$)']) - np.log10(a.empiricalMdot(simulated['Mass (M$_\\odot$)']))
    ax.axhline(0, color='tomato')
    ax.set_xlabel('log(Mass) (M$_\odot$)', size=16)
    ax.set_ylabel('Residual (log space)', size=16)
    ax.scatter(np.log10(observed['Object Mass M_Solar']), observed_difference, color='black', marker='x', alpha=0.4, label = 'Observed Residuals')
    ax.scatter(np.log10(simulated['Mass (M$_\\odot$)']), simulated_difference, color='darkseagreen', alpha=0.4, label = 'Simulated Residuals')
    ax.legend(frameon=True)
    
    return fig
    
    
    
def MarginalDistribution(observed, simulated):
    #load in data
    logmass = np.log10(simulated['Mass (M$_\\odot$)'])
    logMdot = np.log10(simulated['Mdot (M$_\\odot$)'])
    logaccepted_relation = np.log10(simulated['"true" Mdot (M$_\\odot$)'])
    observed_mass = np.log10(observed['Object Mass M_Solar'])
    observed_mdot = np.log10(observed['Accretion Rate M_solar yr-1'])

    #Figure Dimensions
    left, width = 0.1, 0.7
    bottom, height = 0.3, 0.55
    spacing = 0.0075

    #create figure
    fig = plt.figure(figsize=(16, 12))

    #Mass vs Mdot
    rect_frame1 = [left, bottom, width, height]

    #Marginal Distributions
    rect_histx = [left, bottom + height + spacing, width, 0.1]
    rect_histy = [left + width + spacing, bottom, 0.1, height]

    #create axes
    frame1 = fig.add_axes(rect_frame1)
    ax_histx = fig.add_axes(rect_histx, sharex=frame1)
    ax_histy = fig.add_axes(rect_histy, sharey=frame1)
    ax = plt.gca()

    #Marginal Distributions
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # determine nice limits by hand:
    binwidth = 0.25
    xymax = max(np.max(np.abs(observed_mass)), np.max(np.abs(observed_mdot)))
    lim = (int(xymax/binwidth) + 1) * binwidth
    bins = np.arange(-lim, lim + binwidth, binwidth)

    #Plot histograms
    ax_histx.hist(observed_mass, bins=bins, color='black', density=True, alpha=0.225, label='Observed')
    ax_histy.hist(observed_mdot, bins=bins, orientation='horizontal', color='black', density=True, alpha=0.225, label='Observed')
    ax_histx.hist(logmass, bins=bins, color='darkseagreen', density=True, alpha=0.3, label='Simulated')
    ax_histy.hist(logMdot, bins=bins, orientation='horizontal', color='darkseagreen', density=True, alpha=0.3, label='Simulated')
    ax_histx.legend()
    ax_histy.legend()

    #Mass vs Mdot
    frame1.scatter(logmass, logMdot, color='darkseagreen', s=130, alpha=0.4, label='Simulated')
    frame1.scatter(observed_mass, observed_mdot, color='black', s = 50, marker='x', alpha=0.4, label='Observed')
    frame1.plot(logmass, logaccepted_relation, color='tomato', label='Empirical Relationship')
    frame1.axvline(x=np.log10(0.012), color='black', linestyle='-.', label='Deuterium Burning Limit')
    frame1.axvline(x=np.log10(0.1), color='black', linestyle=':', label='Hydrogen Burning Limit')
    frame1.set_xlim(-2.6, 0.6)
    frame1.set_ylim(-14, -5.75)
    frame1.set_xlabel('log(Mass) (M$_\odot$)', size=16)
    frame1.set_ylabel('log(Mass Accretion Rate) (M$_\odot$/yr)', size=16)
    frame1.text(-2.3, -6.6, 'Planets', size = 15, fontstyle='italic')
    frame1.text(-1.64, -6.6, 'Brown Dwarfs', size = 15, fontstyle='italic')
    frame1.text(-0.29, -6.6, 'Stars', size = 15, fontstyle='italic')
    ax_histx.set_title('Monte Carlo Error Propagation', size=20, fontweight='heavy')
    frame1.legend(loc='lower right', prop={'family':'sans-serif', 'style':'normal', 'size': 16}, frameon=False, shadow=True)

    #Residual Plot
    rect_frame2 = [left, bottom-0.2, width, bottom-0.1-spacing]
    rect_ax2_hist = [left+width+spacing, bottom-0.2, 0.1, bottom-0.1-spacing]

    #create axes
    frame2 = fig.add_axes(rect_frame2)  
    ax2_hist = fig.add_axes(rect_ax2_hist) 
    ax = plt.gca()

    #Calculate Residuals
    observed_difference = np.log10(observed['Accretion Rate M_solar yr-1']) - np.log10(a.empiricalMdot(observed['Object Mass M_Solar']))
    simulated_difference = np.log10(simulated['Mdot (M$_\\odot$)']) - np.log10(a.empiricalMdot(simulated['Mass (M$_\\odot$)']))

    #Residual Marginal Distribution
    #remove ticks
    ax2_hist.tick_params(axis="y", labelleft=False)

    #bins by hand
    residualbinwidth = 0.25
    residualmax = np.max(np.abs(observed_difference))
    residuallim = (int(residualmax/residualbinwidth) + 1) * residualbinwidth
    residualbins = np.arange(-residuallim, residuallim + residualbinwidth, residualbinwidth)
    
    
    #plot
    ax2_hist.hist(observed_difference, bins=residualbins, orientation='horizontal', color='black', density=True, alpha=0.225, label='Observed')
    ax2_hist.hist(simulated_difference, bins=residualbins, orientation='horizontal', color='darkseagreen', density=True, alpha=0.3, label='Simulated')
    ax2_hist.set_ylim(-3.3, 3.2)
    ax2_hist.axhline(0, color='tomato')
    ax2_hist.legend()

    #Residual Plot
    frame2.axhline(0, color='tomato')
    frame2.set_xlabel('log(Mass) (M$_\odot$)', size=16)
    frame2.set_ylabel('Residual (log space)', size=16)
    frame2.set_ylim(-3.3, 3.2)
    frame2.scatter(np.log10(observed['Object Mass M_Solar']), observed_difference, color='black', marker='x', alpha=0.4, label = 'Observed Residuals')
    frame2.scatter(np.log10(simulated['Mass (M$_\\odot$)']), simulated_difference, color='darkseagreen', alpha=0.4, label = 'Simulated Residuals')
    frame2.legend(frameon=True)
    
    return fig
    
    
    
def residuals(observed, simulated):
    #load in data
    logmass = np.log10(simulated['Mass (M$_\\odot$)'])
    logMdot = np.log10(simulated['Mdot (M$_\\odot$)'])
    logaccepted_relation = np.log10(simulated['"true" Mdot (M$_\\odot$)'])
    observed_mass = np.log10(observed['Object Mass M_Solar'])
    observed_mdot = np.log10(observed['Accretion Rate M_solar yr-1'])

    #Figure Dimensions
    left, width = 0.1, 0.4
    bottom, height = 0.1, 0.7
    spacing = 0.0075
    
    #create figure
    fig = plt.figure(figsize=(16, 12))
    
    #Residual Plot
    rect_frame2 = [left, bottom, width, height]
    rect_ax2_hist = [left+width+spacing, bottom, 0.3, height]

    #create axes
    frame2 = fig.add_axes(rect_frame2)  
    ax2_hist = fig.add_axes(rect_ax2_hist) 
    ax = plt.gca()

    #Calculate Residuals
    observed_difference = np.log10(observed['Accretion Rate M_solar yr-1']) - np.log10(a.empiricalMdot(observed['Object Mass M_Solar']))
    simulated_difference = np.log10(simulated['Mdot (M$_\\odot$)']) - np.log10(a.empiricalMdot(simulated['Mass (M$_\\odot$)']))

    #Residual Marginal Distribution
    #remove ticks
    ax2_hist.tick_params(axis="y", labelleft=False)

    #bins by hand
    residualbinwidth = 0.25
    residualmax = np.max(np.abs(observed_difference))
    residuallim = (int(residualmax/residualbinwidth) + 1) * residualbinwidth
    residualbins = np.arange(-residuallim, residuallim + residualbinwidth, residualbinwidth)
    
    #plot
    ax2_hist.hist(observed_difference, bins=residualbins, orientation='horizontal', color='black', density=True, alpha=0.225, label='Observed')
    ax2_hist.hist(simulated_difference, bins=residualbins, orientation='horizontal', color='darkseagreen', density=True, alpha=0.3, label='Simulated')
    ax2_hist.set_ylim(-3.3, 3.2)
    ax2_hist.axhline(0, color='tomato', label='Empirical Relationship')
    ax2_hist.legend()

    #Residual Plot
    frame2.axhline(0, color='tomato', label='Empirical Relationship')
    frame2.set_xlabel('log(Mass) (M$_\odot$)', size=16)
    frame2.set_ylabel('Residual (log space)', size=16)
    frame2.set_ylim(-3.3, 3.2)
    frame2.scatter(np.log10(observed['Object Mass M_Solar']), observed_difference, color='black', marker='x', alpha=0.4, label = 'Observed Residuals')
    frame2.scatter(np.log10(simulated['Mass (M$_\\odot$)']), simulated_difference, color='darkseagreen', alpha=0.4, label = 'Simulated Residuals')
    frame2.legend(frameon=True)
    
    return fig
    
    
    
# Getting log'd bins
def getLogBins(data, num_bins):
    hist, bins, _ = plt.hist(data, bins=num_bins)
    plt.close()
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    
    return logbins
    
    
    
def histogram(observed, simulated):
    #log bins calculation
    logLaccbins = getLogBins((10**observed['Log Accretion Luminosity (solar)']), 30)
    logMdotbins = getLogBins(observed['Accretion Rate M_solar yr-1'], 30)

    fig, axs = plt.subplots(2,2, figsize=(15,12))
    fig.subplots_adjust(hspace=0.34)
    ax = plt.gca()
    #radius
    radhist, radbins, _ = axs[0,0].hist(observed['Object Radius R_solar'], bins=20, alpha=0.15, density=True, label='Observed')
    axs[0,0].hist(simulated['Radius (R$_\\odot$)'], bins=radbins, color='darkorange', alpha=0.15, density=True, label='Simulated')
    axs[0,0].set_xlabel('Radius (R$_\\odot)$', fontsize=20)
    axs[0,0].set_title('Radius (R$_\\odot)$', fontsize=40)
    axs[0,0].legend(fontsize=20)
    #temperature
    temphist, tempbins, _ = axs[0,1].hist(observed['Effective Temperature K'], bins=20, alpha=0.2, density=True, label='Observed')
    axs[0,1].hist(simulated['Teff (K)'], bins=tempbins, alpha=0.2, density=True, label='Simulated')
    axs[0,1].set_xlabel('Teff (K)', fontsize=20)
    axs[0,1].legend(fontsize=20)
    axs[0,1].set_title('Teff (K)', fontsize=40)
    #Lacc
    axs[1,0].hist(10**observed['Log Accretion Luminosity (solar)'], bins=logLaccbins, alpha=0.2, density=True, label='Observed')
    axs[1,0].hist(simulated['Lacc (L$_\\odot$)'], bins=logLaccbins, color='darkorange', alpha=0.2, density=True, label='Simulated')
    axs[1,0].set_xlabel('Lacc (L$_\\odot$)', fontsize=20)
    axs[1,0].set_xscale('log')
    axs[1,0].set_yscale('log')
    axs[1,0].legend(fontsize=20)
    axs[1,0].set_title('Lacc (L$_\\odot$)', fontsize=40)
    #Mdot
    axs[1,1].hist(observed['Accretion Rate M_solar yr-1'], bins=logMdotbins, alpha=0.2, density=True, label='Observed')
    axs[1,1].hist(simulated['Mdot (M$_\\odot$)'], bins=logMdotbins, color='darkorange', alpha=0.2, density=True, label='Simulated')
    axs[1,1].set_xlabel('Mdot (M$_\\odot$/yr)', fontsize=20)
    axs[1,1].set_title('Mdot (M$_\\odot$/yr)', fontsize=40)
    axs[1,1].set_xscale('log')
    axs[1,1].set_yscale('log')
    axs[1,1].legend(fontsize=20)

    return fig