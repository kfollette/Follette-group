pro scl_corr

;;;calculate matrix with pearson correlation coefficients between scale factors, airmass, seeing and wavefront error

visao_inventory, /stp

am_sciims=am[sciims]
mag1fwhm_sciims=mag1fwhm[sciims] 
writefits, 'airmass_preproc.fits', am_sciims
writefits, 'mag1seeing_preproc.fits', mag1fwhm_sciims  

scl=readfits('scale_factors.fits')
wfe=readfits('avgwfe_preproc.fits')
am=readfits('airmass_preproc.fits')
mag1=readfits('mag1seeing_preproc.fits')

mat=transpose([[scl],[wfe],[am],[mag1]])
print, correlate(mat)

end