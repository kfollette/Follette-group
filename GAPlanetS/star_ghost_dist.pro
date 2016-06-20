pro star_ghost_dist

x=readfits('./fits/Cont_star_xcen_array.fits', /silent)
y=readfits('./fits/Cont_star_ycen_array.fits', /silent)
x2=readfits('./fits/Cont_ghost_xcen_array.fits', /silent)
y2=readfits('./fits/Cont_ghost_ycen_array.fits', /silent)

r = sqrt((x-x2)^2+(y-y2)^2)
bad=where(r-median(r) gt 2)
r[bad]='NaN'

rmed=median(r)
rstdev=stddev(r, /nan)

print, 'median Continuum star-ghost separation', rmed, 'stdev', rstdev

x3=readfits('./fits/Ha_star_xcen_array.fits', /silent)
y3=readfits('./fits/Ha_star_ycen_array.fits', /silent)
x4=readfits('./fits/Ha_ghost_xcen_array.fits', /silent)
y4=readfits('./fits/Ha_ghost_ycen_array.fits', /silent)

r2 = sqrt((x3-x4)^2+(y3-y4)^2)
bad2=where(r2-median(r2) gt 2)
r2[bad2]='NaN'

rmed2=median(r2)
rstdev2=stddev(r2, /nan)

print, 'median Ha star-ghost separation', rmed2, 'stdev', rstdev2

end