pro unsat_phot

  files = FILE_SEARCH('*.fits',count=nfiles)
  print,'name, exptime, star flux'

  fluxarr = []

  FOR i=0,nfiles-1 DO BEGIN
     
     name = files[i]

     header=HEADFITS(name)
     sxdelpar,header,'NAXIS3'
     MODFITS,name,0,header

     im=MRDFITS(name,0,header,/FSCALE,/SILENT)

     exptime = FXPAR(header, 'EXPTIME')

     starx=FXPAR(header,'CRPIX1A')
     stary=FXPAR(header,'CRPIX2A')

     pad = 5.0 
     sub_im = im[starx-pad:starx+pad, stary-pad:stary+pad]
     fit = MPFIT2DPEAK(sub_im, A, LORENTZIAN, TILT)
     
     objfwhm = A[2] + A[3]

     ; generic: app = 2.5, skya = 3.0, skyb = 6.0
     app = 2.5
     skya = 3.0
     skyb = 6.0

     APER,im,starx,stary,star_flux,0,0,0,0,app*objfwhm,[skya*objfwhm,skyb*objfwhm],[0,0],/EXACT,/FLUX,/SILENT,/NAN

     print, name + ', ' + STRN(exptime) + ', '+ STRN(star_flux)
     
     fluxarr = [fluxarr, star_flux]

  ENDFOR
  
print, 'Median flux for aperture of radius ', STRN(app*objfwhm), " pixels: ", median(fluxarr)  
print, 'Flux per second (in ADU):', median(fluxarr)/exptime
  
END

