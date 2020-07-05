#!/usr/bin/env python
# coding: utf-8

# # FolletteLab MagAOX Target Visibility Tool #

# > Goal: Create a prototype web tool that for each inputted object returns plots for airmass and location, as well as visibility statistics for a given night of observation.
# 
# If a module in not currently in your Jupyter Notebook kernel, Python installation, or are just having import issues, execute the following two lines code:
# 
# For example, if you did not have astroplan
# 
# >import sys
# 
# >!{sys.executable} -m pip install astroplan

# ### Imports ###

# In[ ]:

import numpy as np, pandas as pd, scipy.stats as st, matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
matplotlib.use('Agg')



from astroquery.simbad import Simbad

from astropy.io import ascii
from astropy import constants as const
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import *
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
from astropy.table import Table

from astroplan import (AltitudeConstraint, AirmassConstraint, AtNightConstraint)
from astroplan import FixedTarget, Observer
from astroplan import is_observable, is_always_observable, months_observable
from astroplan import observability_table
from astroplan import FixedTarget, Observer




from astropy.coordinates import EarthLocation
from astroplan.plots import plot_airmass
from astroplan.plots import plot_sky
from astropy.time import Time
import scipy.optimize as optimization
from scipy.optimize import curve_fit
from scipy.stats import linregress
from sklearn.neighbors import KernelDensity
from subprocess import *

from flask import Flask
from flask import request, redirect, render_template

def plot(start_Date,end_Date,objects):
    # In[ ]:

    print(objects)
    #This tool is designed for the Magellan Telescope @ Las Camapanas Observatory, in Chile
    las = Observer.at_site('LCO')
    #both las and los are the locations of MagAO, but one is used for the plot and the other for the time
    lco = EarthLocation.of_site('Las Campanas Observatory')


    # In[ ]:


    #start_time = Time('2019-07-18 21:00:00', location = lco)
    #end_time = Time('2019-07-19 13:00:00', location = lco)
    start_time = Time(start_Date, location = lco)
    end_time = Time(end_Date, location = lco)
    date = str(start_time)[:10] + ' to ' + str(end_time)[:10]


    # In[ ]:


    delta_t = end_time - start_time
    observe_time = start_time + delta_t*np.linspace(0, 1, 75)


    # ### Stars ###

    # In[ ]:

    userEntered_list = list(objects.split(","))
    target_list = userEntered_list


    # In[ ]:


    x = 0
    for i in target_list:
        target_list[x] = FixedTarget.from_name(target_list[x])
        x += 1


    # ### Airmass Curve ###

    # In[ ]:


    """
    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['legend.loc'] = "Upper Right"
    plt.rcParams['figure.dpi'] = 300


    x = 0
    for i in target_list:
        plot_airmass(target_list[x], las, observe_time, max_airmass = 2.25, brightness_shading=True, altitude_yaxis=True)
        x += 1


    plt.legend(shadow=True, loc='center left', bbox_to_anchor=(1.1, 1))
    plt.title('Airmass Curve for Night of ' + date)

    plt.savefig('airmass.png', bbox_inches='tight', pad_inches=0.25)
    """

    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['legend.loc'] = "Upper Right"
    plt.rcParams['figure.dpi'] = 300


    x = 0
    for i in target_list:
        plot_airmass(target_list[x], las, observe_time, max_airmass = 2.25, brightness_shading=True, altitude_yaxis=True)
        x += 1


    plt.legend(shadow=True, loc='center left', bbox_to_anchor=(1.1, 1))
    plt.title('Airmass Curve for Night of ' + date)

    #os.remove("static/airmass.png") idk if necesarry
    plt.savefig('static/airmass.png', bbox_inches='tight', pad_inches=0.25)
    plt.clf()
    # In[ ]:

def vis(start_Date,end_Date,objects,obj_tab):
	
	#This tool is designed for the Magellan Telescope @ Las Camapanas Observatory, in Chile
	las = Observer.at_site('LCO')
	#both las and los are the locations of MagAO, but one is used for the plot and the other for the time
	lco = EarthLocation.of_site('Las Campanas Observatory')
	
	userEntered_list = list(objects.split(","))
	target_list = userEntered_list
	
	targets = []
	for i in range(1, len(obj_tab)):
		ra=(obj_tab.iloc[i, 2])[1:] + ' hours'
		dec=(obj_tab.iloc[i, 3])[1:] + ' degrees'
		print(ra + ',' + dec)
		targets.append(FixedTarget(coord=SkyCoord(ra=ra, dec=dec), name=target_list[i-1]))
	
	constraints = [AltitudeConstraint(10*u.deg, 80*u.deg),
               AirmassConstraint(5), AtNightConstraint.twilight_civil()]
	
	start_time = Time(start_Date, location = lco)
	end_time = Time(end_Date, location = lco)
	date = str(start_time)[:10] + ' to ' + str(end_time)[:10]
	
	time_range = Time([start_Date, end_Date])
	
	
    # In[ ]:
	
	
	delta_t = end_time - start_time
	observe_time = start_time + delta_t*np.linspace(0, 1, 75)
    
	# In[ ]:
	
	
	# Are targets *ever* observable in the time range?
	ever_observable = is_observable(constraints, las, targets, time_range=time_range)
	
	# Are targets *always* observable in the time range?
	always_observable = is_always_observable(constraints, las, targets, time_range=time_range)
	
	# During what months are the targets ever observable?
	best_months = months_observable(constraints, las, target_list)
	
	
	# In[ ]:
	
	
	table = observability_table(constraints, las, targets, time_range=time_range)
	print(table)
	
	table = table.to_pandas()
	
	np.savetxt('visibility.txt', table, fmt="%-30s", header = 'Target name                  ever observable                always observable              fraction of time observable')

    #plot_airmass?


    # In[ ]:


    #plt.legend?


    # In[ ]:




