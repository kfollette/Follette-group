pro cfhtarchiveflat,flatlist,darklist,flag
  
  ;Two cases
  ;1 - Dome flats
  ;2 - Sky flats


  ;Simple case

     
  nflat=0
  flatlist=FILELIST(flatlist,nflat)
  nflat_frames=0
  FOR i=0,nflat-1 DO BEGIN
     nflat_frames+=FXPAR(HEADFITS(flatlist[i]),'NAXIS3')
  ENDFOR
  IF nflat_frames eq nflat THEN BEGIN
     flats=FLTARR(1024,1024,nflat_frames)
     FOR i=0,nflat_frames-1 DO BEGIN
        flats[*,*,i]=MRDFITS(flatlist[i],0,header,/FSCALE,/SILENT)
     ENDFOR
  ENDIF ELSE BEGIN
     count=0
     flats=FLTARR(1024,1024,nflat_frames)
     FOR i=0,nflat-1 DO BEGIN
        tmp=MRDFITS(flatlist[i],0,header,/FSCALE,/SILENT)
        naxis=FXPAR(header,'NAXIS3')
        flats[*,*,count:count+naxis-1]=tmp
        count+=naxis
     ENDFOR
  ENDELSE
  IF flag eq 0 THEN BEGIN
     ndark=0
     darklist=FILELIST(darklist,ndark)
     darks=FLTARR(1024,1024,ndark)
     FOR i=0,ndark-1 DO BEGIN
        darks[*,*,i]=MRDFITS(darklist[i],0,header,/FSCALE,/SILENT)
     ENDFOR
     final_dark=MEDIAN(darks,DIMENSION=3)


     FOR i=0,nflat_frames-1 DO BEGIN
        flats[*,*,i]-=final_dark
     ENDFOR
     final_flat=MEDIAN(flats,DIMENSION=3)
     final_flat/=MEDIAN(final_flat)

  ENDIF ELSE IF flag eq 1 THEN BEGIN
     
     ;flat is the gradient of pixel value vs median count
     medvalues=FLTARR(nflat_frames)
     FOR i=0,nflat_frames-1 DO BEGIN
        medvalues[i]=MEDIAN(flats[*,*,i])
     ENDFOR
     final_flat=FLTARR(1024,1024)

     FOR i=0,1023 DO BEGIN
        FOR j=0,1023 DO BEGIN
           A=POLY_FIT(medvalues,flats[i,j,*],1)
           final_flat[i,j]=A[1]
        ENDFOR
     ENDFOR

     final_flat/=MEDIAN(final_flat)

  ENDIF

  bpm=FLTARR(1024,1024)+1.0
  flatmin=0
  flatmax=0
  FOR i=1.0,0.0,-0.01 DO BEGIN
     npixels=N_ELEMENTS(WHERE((final_flat ge i-0.01) and (final_flat lt i)))
     IF (npixels lt 500) and (flatmin eq 0) THEN flatmin=i-0.01
  ENDFOR
  FOR i=1.0,2.0,0.01 DO BEGIN
     npixels=N_ELEMENTS(WHERE((final_flat ge i) and (final_flat lt i+0.01)))
     IF (npixels lt 500) and (flatmax eq 0) THEN flatmax=i
  ENDFOR
  bpm[WHERE((final_flat lt flatmin) or (final_flat gt flatmax))]=0.0
  
  
  
  flat_final=CFHTMASK(final_flat)
  bpm=CFHTMASK(bpm)
  bpm[WHERE(bpm ne bpm)]=1.0
  SXDELPAR,header,'NAXIS3'
  MWRFITS,flat_final,'../Flat.fits',header,/CREATE
  MWRFITS,bpm,'../BPM.fits',header,/CREATE

END
