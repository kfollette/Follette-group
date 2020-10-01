import numpy as np, pandas as pd, matplotlib.pyplot as plt
plt.rcdefaults()

iprime = pd.read_fwf('VisAO_ip_filter_curve.dat', sep=" ", skiprows = 4)
zprime = pd.read_fwf('VisAO_zp_filter_curve.dat', sep=" ", skiprows = 4)
Ys = pd.read_fwf('VisAO_Ys_filter_curve.dat', sep=" ", skiprows = 4)

curves = [iprime,zprime,Ys]

for x in curves:
	temporary_split = x['RSR       RSR+atm'].str.split(expand=True,)

	x['RSR'] = temporary_split[0]
	x['RSR+atm'] = temporary_split[1]

	#iprime = iprime.drop(['#'], axis=1)
	#iprime = pd.to_numeric(iprime)

	x['RSR'] = pd.to_numeric(x['RSR'])

#plt.plot(iprime['lambda'], iprime['RSR+atm'])
plt.plot(iprime['lambda'], iprime['trans'], label='i\' trans')
plt.plot(iprime['lambda'], iprime['RSR'], label='i\' RSR')
#plt.plot(iprime['lambda'], iprime['RSR+atm'], label='i\' RSR+atm')
plt.plot(zprime['lambda'], zprime['trans'], label='z\' trans')
plt.plot(zprime['lambda'], zprime['RSR'], label='z\' RSR')
#plt.plot(zprime['lambda'], zprime['RSR+atm'], label='z\' RSR+atm')
plt.plot(Ys['lambda'], Ys['trans'], label='YS trans')
plt.plot(Ys['lambda'], Ys['trans'], label='YS RSR')
#plt.plot(Ys['lambda'], Ys['trans'], label='YS RSR+atm')

plt.legend()

plt.show()

#x = np.fromfile("VisAO_ip_filter_curve.dat",dtype=dt)
#x = iprime = pd.read_fwf(‘VisAO_ip_filter_curve.dat’, sep=“\s+“, skiprows = 4)

#temporary_split = iprime['RSR       RSR+atm'].str.split(expand=True,)

#iprime['RSR'] = temporary_split[0]
#iprime['RSR+atm'] = temporary_split[1]

#iprime['RSR'] = pd.to_numeric(iprime['RSR'])
