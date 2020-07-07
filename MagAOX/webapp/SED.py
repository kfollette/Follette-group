#!/usr/bin/env python
# coding: utf-8

# # FolletteLab SED/Relevant Photometric Data Tool #

# > Create a tool that will find relevant data and display guide star magnitude and an SED for every object in a Magellan-formatted catalog
# 
# > For all objects in catalog, pull photometry and spectral type from databases. Overplot photometry on a template for that spectral type.
# 
# > If a spectral template does not exist for this object, overplot photometry over a classic blacbody curve
# 
# If a module in not currently in your Jupyter Notebook kernel, execute the following code:
# 
# For example, if you did not have astroplan
# 
# >import sys
# 
# >!{sys.executable} -m pip install astroplan

# In[1]:


#Known Error: Ap (ApSi) spectral type for star HD 167356

#http://classic.sdss.org/dr7/algorithms/spectemplates/
#SDSS Templates
#Warning

#These templates are from the SDSS-I/II spectroscopic pipeline (DR7 and earlier). SDSS-III/IV (DR8 and later) spectroscopic processing pipelines use different templates.
#from astroquery


# ### Imports ###

# In[2]:


import numpy as np, pandas as pd, scipy.stats as st, matplotlib.pyplot as plt

plt.rcParams['font.family'] = "Times New Roman"
plt.rcParams['legend.loc'] = "Upper Right"
plt.rcParams['figure.dpi'] = 900

from astroquery.simbad import Simbad

from specutils import Spectrum1D

from astropy.io import ascii
from astropy import constants as const
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import *
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
from astropy.table import Table

import scipy.optimize as optimization
from scipy.optimize import curve_fit
from scipy.stats import linregress
from sklearn.neighbors import KernelDensity
from subprocess import *


# ### Initializing Simbad and Vizier Queries ###

# In[3]:

#Initializes Simbad queryability
customSimbad = Simbad()

#Fields we wish to query
customSimbad.add_votable_fields('sptype')
customSimbad.add_votable_fields("fluxdata(B)")
customSimbad.add_votable_fields("fluxdata(V)")
customSimbad.add_votable_fields("fluxdata(R)")
customSimbad.add_votable_fields("fluxdata(I)")
customSimbad.add_votable_fields("fluxdata(J)")
customSimbad.add_votable_fields("fluxdata(K)")

#Fields we do not need
customSimbad.remove_votable_fields("coordinates")

#Courtesy of gist.github.com/mfouesneau, an astrophd dude in Germany
try:  # python 3
    from io import BytesIO
    from http.client import HTTPConnection
except ImportError:  # python 2
    from StringIO import StringIO as BytesIO
    from httplib import HTTPConnection

from astropy.table import Table

def query_sed(pos, radius=5):
    """ Query VizieR Photometry 
    The VizieR photometry tool extracts photometry points around a given position 
    or object name from photometry-enabled catalogs in VizieR.
    
    The VizieR photometry tool is developed by Anne-Camille Simon and Thomas Boch
    .. url:: http://vizier.u-strasbg.fr/vizier/sed/doc/
    
    Parameters
    ----------
    pos: tuple or str
        position tuple or object name
    radius: float
        position matching in arseconds.
    
    Returns
    -------
    table: astropy.Table
        VO table returned by the Vizier service.
        
    >>> query_sed((1.286804, 67.840))
    >>> query_sed("HD1")
    """
    try:
        ra, dec = pos
        target = "{0:f},{1:f}".format(ra, dec)
    except:
        target = pos
    
    url = "http:///viz-bin/sed?-c={target:s}&-c.rs={radius:f}"
    host = "vizier.u-strasbg.fr"
    port = 80
    path = "/viz-bin/sed?-c={target:s}&-c.rs={radius:f}".format(target=target, radius=radius)
    connection = HTTPConnection(host, port)
    connection.request("GET", path)
    response = connection.getresponse()
   
    table = Table.read(BytesIO(response.read()), format="votable")
    return table


# ### Targets ###

# In[4]:

def gen(objects):
	#entered_list = ['Vega', 'Sirius', 'HD 79870', 'ieguieiuiug', 'not a real star']
	entered_list = list(objects.split(","))
	
	
	# ### Checking if Names are Resolvable, returns Resolvable Targets ###
	
	# In[5]:
	
	
	#Target names, may be deprecated
	target_list = []
	
	#Target data
	star_list = []
	
	x = 0
	
	for i in entered_list:
	    star = customSimbad.query_object(entered_list[x])
	    if type(star) == type(None):
	        print('\033[1m' +'ERROR: The star ||' 
	              + entered_list[x] + '|| was not resolvable in Simbad. Please check the spelling of the name.')
	        print('\033[0m')
	    else:
	        target_list.append(entered_list[x])
	        star_list.append(star)
	        
	        
	    x += 1

	target_list


	# ### Dedicated DataFrame for Photometric Magnitude Values ###

	# In[12]:


	#Neither camera can look at wavelengths <600nm, so start w/ B-band for estimation if V-band is empty

	Simbad_Mags = Simbad()
	Simbad_Mags.remove_votable_fields("coordinates")
	Simbad_Mags.add_votable_fields("fluxdata(B)","fluxdata(V)","fluxdata(R)","fluxdata(I)",
									"fluxdata(J)","fluxdata(H)","fluxdata(K)")


	# In[13]:


	#Let's collect dataframes for all the targets, then concatenate them to one big one
	magnitudes = []

	x = 0
	for i in target_list:
		star = customSimbad.query_object(target_list[x])
		star = star.to_pandas()
		magnitudes.append(star)
		x += 1
	
	#magnitudes[0]


	# In[14]:


	magnitudes_result = pd.concat(magnitudes, ignore_index=True)
	magnitudes_result


	# In[15]:


	#np.savetxt("mag_sources.txt", magnitudes_result)


	# In[ ]:





	# In[ ]:





	# In[ ]:





	# ### CSVs Containing Spectral Type Data, to later get approx. temp, absolute magnitude, and solar luminosity ###

	# In[16]:


	main_key = pd.read_csv('CSV_Files/mainsequence.csv')
	giantsiii_key = pd.read_csv('CSV_Files/giantsiii.csv')
	supergiantsi_key = pd.read_csv('CSV_Files/supergiantsi.csv')

	# Setting indexer as spectral type
	main_key.set_index('SP_TYPE', inplace=True)
	giantsiii_key.set_index('SP_TYPE', inplace=True)
	supergiantsi_key.set_index('SP_TYPE', inplace=True)


	# ### Spectral Type ###

	# In[17]:


	columns = {'MAIN_ID': [0], 'SP_TYPE': [0]}
	vis = pd.DataFrame(data = columns)

	x = 0
	for i in target_list:
		obj = customSimbad.query_object(target_list[x])
		obj.keep_columns(['MAIN_ID', 'SP_TYPE'])
		obj = obj.to_pandas()
		vis = vis.append(obj)
		x += 1
	
	
	vis = vis.iloc[1:]

	vis = vis.values.tolist()


	# In[18]:


	#Currently in as a bytes type, need to decode it into a string

	x = 0
	for i in target_list:
		vis[x][0] = vis[x][0].decode("utf-8")
		vis[x][1] = vis[x][1].decode("utf-8")
		x += 1


	vis


	# In[19]:




	# ### Temp (k), Absolute Magnitude, Solar Luminosities ###

	# In[20]:


	x = 0
	for i in target_list:
		sp = vis[x][1][0:2]
		app = main_key.loc[sp]
		app = app.tolist()
		print(app)
		vis[x].extend(app)
		x += 1
	

	vis


	# > The list is now [name, spectral type, temp (k), absolute magnitude, solar luminosity]

	# ### Getting Spectral Template ###

	# In[ ]:



	

	# ### SED ###
	
	# In[33]:


	x = 0

	for i in target_list:

		#First, let's get a dataframe of existing photometric data from Vizier
		print(vis[x][0])
		#s = query_sed(vis[x][0].replace(" ", ""))
		#del s['_ID'], s['_tabname'], s['_RAJ2000'], s['_DEJ2000']
		s = {'wavelength' : [0,0.355,0.467,0.616,0.747,0.892,1.031,1.248,1.631,0,0,2.201,0,0], 'sed_freq' : [5.490,0,0,0,0,0,0,2.394,1.802,1.413,1.395,1.364,0.798,0.635], 'sed_flux' : [3.68e-8,3.66e-8,5.41e-8,2.5e-8,1.39e-8,8.32e-9,5.71e-9,2.98e-9,1.16e-9,4.57e-10,4.35e-10,3.95e-10,5.31e-11,2.22e-11]}
		print(str(len(s['wavelength'])) + "," + str(len(s['sed_freq'])) + ',' + str(len(s['sed_flux'])))
		#s = s.to_pandas()
		s = pd.DataFrame(s)

		#Now let us get the spectral template, using first character of spectral type
		spec = Spectrum1D.read("SDSS DR2 Spectral Templates/" + vis[x][1][0] + "(4).fit", format = "SDSS-I/II spSpec")

		#Need to convert the template's x-axis to micrometers, y-axis to janskys, using astropy.units
		template_wavelength = spec.wavelength.to(u.micrometer)


		#spec.new_flux_unit(u.Jy)

		template_flux = spec.flux
		#print(template_flux.unit)

		#Plot template
		plt.plot(template_wavelength, template_flux, label = "Spectral Template for: " + vis[x][1][0] + "Star", color = 'black')
		#plt.xlim(0.2, 1)



		plt.legend(shadow=True, loc='center left', bbox_to_anchor=(1.1, 1))
		plt.xlabel("Wavelength (μm)")
		plt.ylabel("Flux (Jy)")
		plt.title('Spectral Energy Distribution of: ' + target_list[x])


		#Now will plot photometry data from Vizier
		#wavelength micrometer = speed of light / frequency in hz
		s['wavelength'] = (299792458/(s['sed_freq']*(1e+9)))*(1e+6)

		s['wavelength'] = s['wavelength']*u.micrometer
		#s['sed_flux'] = s['sed_flux']*u.Jy


		#print(s["sed_flux"].to_string())


	
		#sed_freq sed_flux sed_eflux sed_filter
		#GHz Jy Jy
		#float64 float32 float32 bytes32
		#1 Jy = 10-23 erg s-1 cm-2 Hz-1 = 10-26 Watts m-2 Hz-1 = FsubV
		#print(vis[x])
		#plt.plot(template_wavelength, template_flux, label = "Spectral Template for: " + vis[x][1][0] + "Star", color = 'black')
		#plt.scatter(s['wavelength'], s['sed_flux']*100, s=5, color = 'red')
		plt.scatter(s['wavelength'], s['sed_flux']*100, s=5, color = 'red')

		#plt.ylim(0, 8)


		plt.xlim(0.3, 1)

		#print(s.to_string())

		#print(s["sed_flux"].to_string())



		plt.savefig('static/'+target_list[x]+'_sed.png', bbox_inches='tight', pad_inches=0.25)


		x += 1	
		
#gen('vega,sirius,gq lupi')
