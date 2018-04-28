from astropy.io import fits

def headers():
    dc = fits.getdata('Cont_clip451_flat_circsymreg_nocosmics.fits')
    hc = fits.getheader('Cont_clip451_flat_circsymreg_nocosmics.fits')
    hc.set('WLENGTH', 0.6428)
    fits.writeto('Cont_clip451_flat_circsymreg_nocosmics.fits', dc, hc, clobber=True)
    dl = fits.getdata('Line_clip451_flat_circsymreg_nocosmics.fits')
    hl = fits.getheader('Line_clip451_flat_circsymreg_nocosmics.fits')
    hl.set('WLENGTH', 0.6564)
    fits.writeto('Line_clip451_flat_circsymreg_nocosmics.fits', dl, hl, clobber=True)