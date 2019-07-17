pro cand_click_sat, unsatflux, narrowband, broadband, age_gyr

  print, "Unsaturated flux from star per second:", unsatflux

  files = FILE_SEARCH('*.fits',count=nfiles)

  ; Get the exptime from the first image in the saturated stack:
  name = files[0]
  header = HEADFITS(name)
  sxdelpar, header, 'NAXIS3'
  MODFITS, name, 0, header
  exptime = FXPAR(header, 'EXPTIME')

  
  IF ((narrowband EQ 'Kc2.09') AND (broadband EQ 'Ks')) THEN BEGIN
     star_flux = 16.4 * unsatflux * exptime
     print, "Estimated saturated flux from star in ", exptime, " seconds: ", star_flux
  ENDIF    
     
  deltamag_array = []

  print,'name, dx, dy, sep, position angle, star flux, candidate flux, delta mag'

  FOR i=0,nfiles-1 DO BEGIN
     
     name=files[i]

     header=HEADFITS(name)
     sxdelpar,header,'NAXIS3'
     MODFITS,name,0,header

     im=MRDFITS(name,0,header,/FSCALE,/SILENT)

     starx=FXPAR(header,'CRPIX1A')
     stary=FXPAR(header,'CRPIX2A')

     
     ok=0
     !v->im,im, frame=1
    
     WHILE (ok eq 0) DO BEGIN
        !v->imexam,x,y
        IF ((FINITE(x) EQ 0) OR (FINITE(y) eq 0)) THEN BEGIN
           ok=1
           CONTINUE
        ENDIF ELSE IF (im[x,y] ne im[x,y]) THEN BEGIN
           ok=1
           CONTINUE
        ENDIF
        
        pad=5.0
        
        ; this measures the fwhm of the primary star; we want the companion instead
        sub_im=im[starx-pad:starx+pad,stary-pad:stary+pad]
        fit=MPFIT2DPEAK(sub_im,A,LORENTZIAN,TILT)
        
        candfwhm=A[2]+A[3]

        sub_im=im[x-pad:x+pad,y-pad:y+pad]
        fit=MPFIT2DPEAK(sub_im,A,LORENTZIAN,TILT)
        
        candfwhm = A[2] + A[3]
        
        xcen=A[4]+x-pad
        ycen=A[5]+y-pad
        
        ; generic: app=2.5, skya=3.0, skyb=6.0
        app=2.5
        skya=3.0
        skyb=6.0
        
        !v->circle,xcen-0.5,ycen-0.5,app*candfwhm,COLOR='green',WIDTH=2,/FIXED
        !v->circle,xcen-0.5,ycen-0.5,skya*candfwhm,COLOR='red',WIDTH=2,/FIXED
        !v->circle,xcen-0.5,ycen-0.5,skyb*candfwhm,COLOR='red',WIDTH=2,/FIXED

        ; don't bother measuring the saturated star flux
        ;APER,im,starx,stary,star_flux,0,0,0,0,app*candfwhm,[skya*candfwhm,skyb*candfwhm],[0,0],/EXACT,/FLUX,/SILENT,/NAN

        APER,im,xcen,ycen,cand_flux,0,0,0,0,app*candfwhm,[skya*candfwhm,skyb*candfwhm],[0,0],/EXACT,/FLUX,/SILENT,/NAN

        dx = -(xcen-starx)
        dy = ycen-stary
        sep = sqrt(dx^2 + dy^2)
        dm = 2.5*alog10(star_flux/cand_flux)
        posang = (180/!pi*atan(dy,-dx)+(270)) mod 360


        print,name+', '+STRN(dx)+', '+STRN(dy)+', '+STRN(sep)+', '+STRN(posang)+', '+STRN(star_flux)+', '+STRN(cand_flux)+', '+STRN(dm)
        
        deltamag_array = [deltamag_array, dm]

     ENDWHILE

     ; remove frame for new image
     !v->cmd, "frame delete"

   ENDFOR

print,'Median delta mag:', median(deltamag_array)

; Start the estimation of the mass from the delta magnitude:

IF ((age eq 5) AND (broadbad eq 'Ks')) THEN BEGIN
	READCOL,'5gyr.txt', msun, teff, l, g, r, mv, mr, mi, mj, mh, mk, ml, mm 
	xinterp = findgen(100, start=min(mk), increment=(max(mk)-min(mk))/100.) 
        result = interpol(msun, mk, xinterp)
        ; find nearest match to delta mag in xinterp range, then match to result
ENDIF

END

