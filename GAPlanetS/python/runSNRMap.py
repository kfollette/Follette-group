"""
used to run SNRMap from the command line, using a fits file input rather than an object input.
COMMAND LINE ARGUMENTS:
    input filename
    True or False, whether you want output fits files of everything to be saved or just the final snr map
    Optional: four integers, corresponding to parameters for a custom wedge mask: starting angle, ending angle, inner radius, outer radius.
        angles are in degrees, counterclockwise with 0 pointing to the right, and the starting angle must be smaller than the ending angle.
"""

import SNRMap as SNR
import SNRMap_helper as SNRhelp
import sys

file = sys.argv[1]
save = sys.argv[2]
if save == 'True':
    save = True
elif save == 'False':
    save = False
if len(sys.argv) > 3:
    theta1 = float(sys.argv[3])
    theta2 = float(sys.argv[4])
    r1 = float(sys.argv[5])
    r2 = float(sys.argv[6])
    mask = SNRhelp.implant_custom_mask(theta1, theta2, r1, r2)
    snr_map, noise_map, multiplied_cube, mask_cube, radial_profile = SNR.create_map(file, save, mask)
else:
    snr_map, noise_map, multiplied_cube, mask_cube, radial_profile = SNR.create_map(file, save)
#med_HD142527_8Apr14short_SDI_a7m3-10KLmodes.fits




