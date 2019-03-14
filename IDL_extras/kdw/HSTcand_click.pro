pro hstcand_click,starx,stary

  files = FILE_SEARCH('*.fits',count=nfiles)
  print,'name, dx, dy, sep, position angle, star flux, candidate flux, delta mag'

  FOR i=0,nfiles-1 DO BEGIN
     
     name=files[i]

     im=MRDFITS(name,0,header,/FSCALE,/SILENT)

     ;starx=FXPAR(header,'CRPIX1A')
     ;stary=FXPAR(header,'CRPIX2A')

     ok=0
     !v->im,im
     
    
     WHILE (ok eq 0) DO BEGIN
        !v->imexam,x,y
        IF ((FINITE(x) EQ 0) OR (FINITE(y) eq 0)) THEN BEGIN
           ok=1
           CONTINUE
        ENDIF ELSE IF (im[x,y] ne im[x,y]) THEN BEGIN
           ok=1
           CONTINUE
        ENDIF
        
        pad=10.0
        sub_im=im[starx-pad:starx+pad,stary-pad:stary+pad]
        fit=MPFIT2DPEAK(sub_im,A,/LORENTZIAN,/TILT)
        
        candfwhm=A[2]+A[3]

        sub_im=im[x-pad:x+pad,y-pad:y+pad]
        fit=MPFIT2DPEAK(sub_im,A,/LORENTZIAN,/TILT)
        
        xcen=A[4]+x-pad
        ycen=A[5]+y-pad
        
        app=2.5
        skya=10.0
        skyb=30.0
        
        !v->circle,xcen-0.5,ycen-0.5,app*candfwhm,COLOR='green',WIDTH=2,/FIXED
        !v->circle,xcen-0.5,ycen-0.5,skya*candfwhm,COLOR='red',WIDTH=2,/FIXED
        !v->circle,xcen-0.5,ycen-0.5,skyb*candfwhm,COLOR='red',WIDTH=2,/FIXED

        

        APER,im,starx,stary,star_flux,0,0,0,0,app*candfwhm,[skya*candfwhm,skyb*candfwhm],[0,0],/EXACT,/FLUX,/SILENT,/NAN

        APER,im,xcen,ycen,cand_flux,0,0,0,0,app*candfwhm,[skya*candfwhm,skyb*candfwhm],[0,0],/EXACT,/FLUX,/SILENT,/NAN

        dx = -(xcen-starx)
        dy = ycen-stary
        sep = sqrt(dx^2 + dy^2)
        dm = 2.5*alog10(star_flux/cand_flux)
        posang = (180/!pi*atan(dy,-dx)+(270)) mod 360


        print,name+', '+STRN(dx)+', '+STRN(dy)+', '+STRN(sep)+', '+STRN(posang)+', '+STRN(star_flux)+', '+STRN(cand_flux)+', '+STRN(dm)

        
     ENDWHILE
  ENDFOR
END

