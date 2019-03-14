pro cfhtarchiveflat,flatlist,flag
;;if flag is zero, then gradient... if flag not zero, then on/off

IF flag eq 0 THEN BEGIN

  
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

ENDIF else begin

 	nflats=0
        flatlist=FILELIST(flatlist,nflats)

        IF nflats ne 2 THEN BEGIN
         print,'Something is going wrong!'
         stop
        ENDIF
        medians=FLTARR(2)
        flata=MRDFITS(flatlist[0],0,header,/FSCALE,/SILENT)
        flatb=MRDFITS(flatlist[1],0,header,/FSCALE,/SILENT)
        print,MEDIAN(flata),MEDIAN(flatb)
        IF MEDIAN(flata) gt MEDIAN(flatb) THEN BEGIN
           flaton=flata
           flatoff=flatb
        ENDIF ELSE BEGIN
           flaton=flatb
           flatoff=flata
        ENDELSE
        flaton=MEDIAN(flaton,DIMENSION=3)
        flatoff=MEDIAN(flatoff,DIMENSION=3)
        flat=flaton-flatoff


        bigflat=flat/MEDIAN(flat)

        ;now create bpm maps
        bpm=FLTARR(1024,1024)+1.0

         flatmin=0
        flatmax=0
        FOR i=1.0,0.0,-0.01 DO BEGIN
           npixels=N_ELEMENTS(WHERE((bigflat ge i-0.01) and (bigflat lt i)))
           IF (npixels lt 500) and (flatmin eq 0) THEN flatmin=i-0.01
        ENDFOR
        FOR i=1.0,2.0,0.01 DO BEGIN
           npixels=N_ELEMENTS(WHERE((bigflat ge i) and (bigflat lt i+0.01)))
           IF (npixels lt 500) and (flatmax eq 0) THEN flatmax=i
        ENDFOR
        bpm[WHERE((bigflat lt flatmin) or (bigflat gt flatmax))]=0.0



        bigflat=CFHTMASK(bigflat)
        bpm=CFHTMASK(bpm)
        bpm[WHERE(bpm ne bpm)]=1.0
        SXDELPAR,header,'NAXIS3'
        MWRFITS,bigflat,'../Flat.fits',header,/CREATE
        MWRFITS,bpm,'../BPM.fits',header,/CREATE

ENDELSE

END
