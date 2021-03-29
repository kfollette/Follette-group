import os
import astropy.io.fits as fits
import pyklip.fitpsf as fitpsf
import matplotlib.pylab as plt
import pyklip.instruments.MagAO as MagAO
import numpy as np
import sys


outputdir = sys.argv[1]
prefix = sys.argv[2]
sep = float(sys.argv[3])
pa = float(sys.argv[4])
length = float(sys.argv[5])
output_prefix = os.path.join(outputdir, prefix)

# get FM frame
fm_frame = fits.getdata(output_prefix + "-fmpsf-KLmodes-all.fits")[0]
fm_header = fits.getheader(output_prefix + "-fmpsf-KLmodes-all.fits")
fm_centx = fm_header['PSFCENTX']
fm_centy = fm_header['PSFCENTY']

# get data_stamp frame
data_frame = fits.getdata(output_prefix + "-klipped-KLmodes-all.fits")[0]
data_header = fits.getheader(output_prefix + "-klipped-KLmodes-all.fits")
data_centx = data_header['PSFCENTX']
data_centy = data_header['PSFCENTY']

# get initial guesses. Should be in the header but aren't?
guesssep = sep
guesspa = pa

# create FM Astrometry object - 13 is fitboxsize
fma = fitpsf.FMAstrometry(guesssep, guesspa, 13)

# generate FM stamp
# padding should be greater than 0 so we don't run into interpolation problems
fma.generate_fm_stamp(fm_frame, [fm_centx, fm_centy], padding=5)

# generate data_stamp stamp
# note that dr=4 means we are using a 4 pixel wide annulus to sample the noise for each pixel
# exclusion_radius excludes all pixels less than that distance from the estimated location of the planet
fma.generate_data_stamp(data_frame, [data_centx, data_centy], dr=4, exclusion_radius=10)

# set kernel, no read noise
corr_len_guess = length
corr_len_label = r"$l$"
fma.set_kernel("matern32", [corr_len_guess], [corr_len_label])

# set bounds
x_range = float(sys.argv[6]) # pixels
y_range = float(sys.argv[7]) # pixels
flux_range = float(sys.argv[8]) # flux can vary by an order of magnitude
corr_len_range = float(sys.argv[9]) # between 0.3 and 30
fma.set_bounds(x_range, y_range, flux_range, [corr_len_range])

import time
t0 = time.time()

# run MCMC fit
fma.fit_astrometry(nwalkers=int(sys.argv[10]), nburn=int(sys.argv[11]), nsteps=int(sys.argv[12]), numthreads=1)

t1 = time.time()
print("time taken: ", str(np.round(t1-t0)), " seconds")

fma.propogate_errs(star_center_err=0.1, platescale=MagAO.MagAOData.lenslet_scale*1000, platescale_err=0.000015, pa_offset=-0.59, pa_uncertainty=0.3)


# show what the raw uncertainites are on the location of the planet
print("\nPlanet Raw RA offset is {0} +/- {1}, Raw Dec offset is {2} +/- {3}".format(fma.raw_RA_offset.bestfit, fma.raw_RA_offset.error,
                                                                                    fma.raw_Dec_offset.bestfit, fma.raw_Dec_offset.error))

# Full error budget included
print("Planet RA offset is at {0} with a 1-sigma uncertainity of {1}".format(fma.RA_offset.bestfit, fma.RA_offset.error))
print("Planet Dec offset is at {0} with a 1-sigma uncertainity of {1}".format(fma.Dec_offset.bestfit, fma.Dec_offset.error))

# Propogate errors into separation and PA space
print("Planet separation is at {0} with a 1-sigma uncertainity of {1}".format(fma.sep.bestfit, fma.sep.error))
print("Planet PA at {0} with a 1-sigma uncertainity of {1}".format(fma.PA.bestfit, fma.PA.error))

print("Flux is {0} with a 1-sigma uncertainty of {1}" .format(fma.raw_flux.bestfit, fma.raw_flux.error))


import pickle

chain_info = pickle.load(open("bka-chain.pkl", "rb"))
fig=plt.figure(figsize=(10,8))
# plot RA offset
ax1 = fig.add_subplot(411)
ax1.plot(chain_info[:,:,0].T, '-', color='k', alpha=0.3)
ax1.set_xlabel("Steps")
ax1.set_ylabel(r"$\Delta$ RA")

# plot Dec offset
ax2 = fig.add_subplot(412)
ax2.plot(chain_info[:,:,1].T, '-', color='k', alpha=0.3)
ax2.set_xlabel("Steps")
ax2.set_ylabel(r"$\Delta$ Dec")

# plot flux scaling
ax3 = fig.add_subplot(413)
ax3.plot(chain_info[:,:,2].T, '-', color='k', alpha=0.3)
ax3.set_xlabel("Steps")
ax3.set_ylabel(r"$\alpha$")

# plot hyperparameters.. we only have one for this example: the correlation length
ax4 = fig.add_subplot(414)
ax4.plot(chain_info[:,:,3].T, '-', color='k', alpha=0.3)
ax4.set_xlabel("Steps")
ax4.set_ylabel(r"$l$")

plt.savefig(outputdir+"/BKA_chain.png")

fig = fma.make_corner_plot()
plt.savefig(outputdir+"/BKA_corner.png")

fig = fma.best_fit_and_residuals()
plt.savefig(outputdir+"/BKA_residuals.png")
