pro cfhtarchivedark,list

  n=0
  list=FILELIST(list,n)

  ;With CFHT one file per dark combo

  FOR i=0,n-1 DO BEGIN
     dark=MRDFITS(list[i],0,header,/FSCALE,/SILENT)
     SXDELPAR,header,'NAXIS3'
     darkname='../Dark-'+STRN(FXPAR(header,'NAXIS1'))+'-'+STRN(FXPAR(header,'INTTIME'))+'.fits'
     IF FXPAR(header,'NAXIS1') eq 1024 THEN MWRFITS,CFHTMASK(MEDIAN(dark,DIMENSION=3)),darkname,header,/CREATE
     IF FXPAR(header,'NAXIS1') ne 1024 THEN  MWRFITS,MEDIAN(dark,DIMENSION=3),darkname,header,/CREATE
  ENDFOR 
END
