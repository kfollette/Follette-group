import numpy as np
import astropy.units as u
from astropy.constants import G, M_jup, R_jup, M_earth, R_earth, L_sun, M_sun, R_sun


def contrast_to_Mdot(dist,star_mag,contr,R, M, munits='Mjup', runits='Rjup', law='aoyama',filtwid=0.006,zeropt=2.339e-5):
	L_line = line_lum(dist,star_mag,contr,filtwid=filtwid,zeropt=zeropt)
	L_acc = acc_lum(L_line,law=law)
	Mdot = acc_rate(L_acc, R, M, munits=munits)
	return(Mdot)

def contrast_to_MMdot(dist,star_mag,contr,R,munits='Mjup', runits='Rjup', law='aoyama',filtwid=0.006,zeropt=2.339e-5):
	L_line = line_lum(dist,star_mag,contr,filtwid=filtwid,zeropt=zeropt)
	L_acc = acc_lum(L_line,law=law)
	Mdot = acc_rate_times_mass(L_acc, R, munits=munits)
	return(Mdot)


def line_lum(dist,star_mag,contr,filtwid=0.006,zeropt=2.339e-5):
	"""
	Given distance, flux zeropoint, stellar R-band magnitude, and H-alpha contrast (or contrast limit),
	calculate line luminosity in Watts. 
	
	REQUIRED INPUTS:
	distance = distance in pc
	star_mag =  star's (extinction corrected) R band magnitude 
	contr = measured contrast at H-alpha

	OPTIONAL INPUTS
	fitwid = width of H-alpha filter in microns (default is VisAO width)
	zeropt = flux zeropoint  in erg/cm2/sec/micron (default is VisAO value)

	"""
	#add units to zeropoint and filter width
	zeropt*=u.erg/u.cm**2/u.s/u.um
	filtwid*=u.um

	#translate contrast to delta mag
	delta_mag = -2.5*np.log10(contr)

	#flux to luminosity
	L_line = 4*np.pi*(dist*u.pc)**2*zeropt*filtwid*10**((star_mag+delta_mag)/-2.5)
	
	#get units in order
	L_line = L_line.decompose().to(u.W)
	return(L_line)

#next two are empirical T-Tauri relationships from Rigliaco 2012

def acc_lum(L_line, law='aoyama'):
	"""
	Translate H-alpha line luminosity to accretion luminosity. 
	From model-based planetary accretion estimates in Aoyama et al. 2021
	From empirical T-Tauri relationships in Alcala et al. 2017
	From empirical  T-Tauri relationships in Rigliaco et al. 2102
	"""
	if law=='aoyama':
		b=1.61
		a=0.95
	elif law=='rigliaco':
		b = 2.27
		a = 1.25
	elif law=='alcala':
		b=1.74
		a=1.13
	else:
		print('please specify a valid accretion scaling relation')
	log_acc = b+a*np.log10(L_line/L_sun)
	L_acc=10**log_acc*L_sun
	return(L_acc)

def acc_rate(L_acc, R, M, munits='Mjup', runits='Rjup'):
	"""
	Translate an accretion luminosity and planet mass/radius to accretion rate in jupiter masses 
	per year following standard relation.

	REQUIRED INPUTS:
	L_acc = Accretion liminosity (with units)
	R = object radius in jupiter radii
	M = object mass in jupiter masses

	OPTIONAL INPUTS:
	munits = defines input/output mass units (options are Msun and Mjup (default) )
	runits = defines input radius units  (options are Msun and Mjup (default) )
	"""
	if munits == "Mjup":
		M*=M_jup
	if munits == "Msun":
		M*=M_sun

	if runits == "Rjup":
		R*=R_jup
	if runits == "Rsun":
		R*=R_sun

	mdot = 1.25*L_acc*R/(G*M)
	if munits=='Mjup':
		mdot = mdot.decompose().to(u.Mjup/u.yr)
	elif munits=='Msun':
		mdot = mdot.decompose().to(u.Msun/u.yr)
	else:
		print('I don\'t recognize this unit')
	return(mdot)

def acc_rate_times_mass(L_acc, R, munits='Mjup', runits='Rjup'):
	"""
	Translate an accretion luminosity and planet radius to an accretion rate in jupiter masses 
	per year * M following standard relation.

	REQUIRED INPUTS:
	L_acc = Accretion liminosity (with units)
	R = object radius in solar or jupiter radii
	M = object mass in solar or jupiter masses

	OPTIONAL INPUTS:
	munits = defines output mass units (options are Msun and Mjup (default) )
	runits = defines input radius units  (options are Msun and Mjup (default) )
	"""
	if runits == "Rjup":
		R*=R_jup
	if runits == "Rsun":
		R*=R_sun

	mmdot = 1.25*L_acc*R/(G)
	if munits=='Mjup':
		mmdot = mmdot.decompose().to(u.Mjup**2/u.yr)
	elif munits=='Msun':
		mmdot = mmdot.decompose().to(u.Msun**2/u.yr)
	else:
		print('I don\'t recognize this unit')
	return(mmdot)