import numpy as np, pandas as pd, matplotlib.pyplot as plt
from specutils import Spectrum1D
from astropy import units as u
import glob as g
import SED
plt.rcdefaults()

def cam1curves():  #process all camsci2 curves at once
	curves = []
	names = []
	for f in g.glob('/MagAOX_filter_curves/camsci1/*.dat'):
		curves.append(process(f))
		names.append(file.replace('magaox','').replace('.dat','').replace('_',' '))
	return (curves,names)
	
def cam2curves(): #process all camsci2 curves at once
	curves = []
	names = []
	for f in g.glob('/MagAOX_filter_curves/camsci2/*.dat'):
		curves.append(process(f))
		names.append(file.replace('magaox','').replace('.dat','').replace('_',' '))
	return (curves,names)
	
def wfscurves(): #process all wfs curves at once
	curves = []
	names = []
	for f in g.glob('/MagAOX_filter_curves/wfs/*.dat'):
		curves.append(process(f))
		names.append(file.replace('magaox','').replace('.dat','').replace('_',' '))
	return (curves,names)
	
def process(file): 
	"""extracts a MagAOX filter curve from a .dat file
	
	Args:
		file (str): the adress of a .dat file containing a filter curve
	Returns:
		pandas.DataFrame: a DataFrame containing all the points in the curve
	"""
	curve = pd.read_fwf(file, sep=" ", skiprows=skippable())
	curve = curve.drop(['#'], axis=1)
	curve['lambda [m]'] = pd.to_numeric(curve['lambda [m]'])
	curve['transmission'] = pd.to_numeric(curve['transmission'])
	return curve
	
def skippable(): #list of rows to be skipped by pd.read_fwf - first 22 lines + line 24 of the .dat file
	l=[i for i in range(0,23)]
	l.append(24)
	return l

def photons_per_sec(star, curve): 
	curve = process(curve)
	template = SED.get_spectral_template(SED.get_vis([star]), 0)
	wavelength = template.wavelength.to(u.micrometer)
	flux = template.flux
	flux_right_units = flux.decompose()
	flux_right_units = flux_right_units * 1*u.m
	flux_right_units = flux_right_units.to(u.watt/(u.m)**2)
	
	interpolated_flux_density = np.interp(curve["lamda(microns)"], wavelength.value, flux)
	trans_interpolated_flux_density = curve["trans(normalized)"] * interpolated_flux_density
	
	plt.plot(wavelength, flux, label = "Initial SED")
	plt.plot(curve["lamda(microns)"], trans_interpolated_flux_density, label = "Transmission")
	plt.xlabel("Wavelength (microns)")
	plt.ylabel("Flux Density")
	plt.legend()
	plt.savefig('static/data/'+star+'_flux_'+curve.replace('magaox','').replace('_filer_curve.dat','')+'.png', bbox_inches='tight',pad_inches=0.25)
	plt.clf()
	
	integrated_value = np.trapz(trans_interpolated_flux_density, curve["lamda(microns)"])
	photons_sec = integrated_value*(3.8*10**9)
	return photons_sec
	
print(photons_per_sec("Vega", 'wfs/magaox_wfs-open_bs-65-35.dat'))