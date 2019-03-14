pro mmtarchivereduction,list,flag,darkflag,alignflag,filter

print,'Have you started ds9?'

  ;if flag is 1, do sky subtraction
  ;if flag is 0, don't
  ;if darkflag eq 1 then do dark subtraction
  

     nims=0                        ;i=number of fits files
     nframes=0                        ;j=number of frames for each star
     imlist=FILELIST(list,nims)
     FOR j=0,nims-1 DO BEGIN
        tmp=MRDFITS(imlist[j],0,header,/SILENT,/FSCALE)
        nframes=nframes+FXPAR(header,'NAXIS3')
     ENDFOR
     s=size(tmp,/DIMENSION)  
     headnum=FLTARR(nframes)
     rawims=FLTARR(s[1],s[1],nframes)
     ims=FLTARR(s[1],s[1],nframes)
     filename=STRARR(nframes)
     name=STRARR(nframes)
     
     jj=0                       ;Use this to count rawframe array
     place=0
     FOR j=0,nims-1 DO BEGIN
        tmp=MRDFITS(imlist[j],0,header,/SILENT,/FSCALE)
        FOR jj=0,FXPAR(header,"NAXIS3")-1 DO BEGIN
           rawims[*,*,place]=tmp[*,*,jj]
           headnum[place]=j
           filename[place]=imlist[j]
           place++
        ENDFOR
     ENDFOR
    
    
     flatname='calib/Flat.fits'
     bpmname='calib/BPM.fits'
     darkname='calib/Dark.fits'
     IF darkflag eq 1 THEN darkframe=MRDFITS(darkname,0,/SILENT,/FSCALE)
     IF darkflag eq 0 THEN darkframe=FLTARR(1024,1024)
     flatframe=MRDFITS(flatname,0,/SILENT,/FSCALE)
     bpmframe=MRDFITS(bpmname,0,/SILENT)
     print,darkname,flatname,bpmname
  

     FOR j=0,nframes-1 DO ims[*,*,j]=(rawims[*,*,j]-darkframe)/flatframe

     IF flag eq 1 THEN BEGIN
        skyframe=MEDIAN(ims,DIMENSION=3)
        MWRFITS,skyframe,'Sky.fits',/CREATE
        ;skyframe=MRDFITS('Sky.fits',0,ssheader,/FSCALE)
        FOR j=0,nframes-1 DO ims[*,*,j]-=skyframe
     ENDIF

     FOR j=0,nframes-1 DO BEGIN
        IF j eq 0 THEN print,nframes,' images'
        FIXPIX,ims[*,*,j],bpmframe,tmp,NPIX=15,/SILENT
        tmp=SIGMA_FILTER(tmp,3,N_SIGMA=3,/ITERATE,/ALL)
        IF s[1] eq 1024 THEN ims[*,*,j]=CFHTMASK(tmp)
        IF s[1] ne 1024 THEN ims[*,*,j]=tmp
     ENDFOR
     

     
;at this point we have large stacks of images, want to reapply the
;header back onto the individual frames, as well as a resonable
;estimate for WCS2 (to be used by diffraction spike algorithm. Need
;human input at this point
     ;skytmp=reducedframe
     FOR j=0,nframes-1 DO BEGIN
        header=HEADFITS(imlist[headnum[j]])
        tmpsmooth=SMOOTH(ims[*,*,j],10,/NAN,/EDGE_TRUNCATE)
        foo=MAX(tmpsmooth,ind,/NAN)
        ind=ARRAY_INDICES(tmpsmooth,ind)
        xpos=ind[0]
        ypos=ind[1]

        IF s[1] lt alignflag*100000. THEN BEGIN
           imedges=[xpos-20,xpos+20,ypos-20,ypos+20]
           IF imedges[0] lt 0 THEN imedges[0]=0
           IF imedges[1] gt s[0]-1 THEN imedges[1]=s[0]-1
           IF imedges[2] lt 0 THEN imedges[2]=0
           IF imedges[3] gt s[1]-1 THEN imedges[3]=s[1]-1
           CNTRD,ims[imedges[0]:imedges[1],imedges[2]:imedges[3],j],20,20,xcen,ycen,3.5,/SILENT ;This is for narrowband
           IF (xcen eq -1) or (ycen eq -1) THEN BEGIN
           ;IF 1 eq 1 THEN BEGIN
              print,xpos,ypos
              !v->im,ims[*,*,j]
              print,'Select star'
              !v->imexam,xcen,ycen
              print,'done!'
              xpos=FLOAT(xcen)
              ypos=FLOAT(ycen)
              CNTRD,ims[*,*,j],xpos,ypos,xcen,ycen,5.5,/SILENT
           ENDIF ELSE BEGIN
              xcen+=imedges[0]
              ycen+=imedges[2]
           ENDELSE
           xpos=xcen
           ypos=ycen
        ENDIF
        

        obj=list
        STRREPLACE,obj,'A',''
        STRREPLACE,obj,'B',''
        STRREPLACE,obj,'C',''

        FXADDPAR,header,'CRPIX1A',xpos
        FXADDPAR,header,'CRPIX2A',ypos
        FXADDPAR,header,'CRVAL1A',0
        FXADDPAR,header,'CRVAL2A',0
        FXADDPAR,header,'CTYPE1A',"Xprime"
        FXADDPAR,header,'CTYPE2A',"Yprime"
        FXADDPAR,header,'CD1_1A',1
        FXADDPAR,header,'CD1_2A',0
        FXADDPAR,header,'CD2_1A',0
        FXADDPAR,header,'CD2_2A',1
        FXADDPAR,header,'BZERO',0
        SXDELPAR,header,'NAXIS3'
        FXADDPAR,header,'OBJECT',obj

        number=STRN(j)
        IF j lt 100 THEN number='0'+STRN(j)
        IF j lt 10 THEN number='00'+STRN(j)

        IF FXPAR(header,'MJDATE') gt 51594.0 THEN BEGIN
		date=FXPAR(header,'DATE-OBS')
	ENDIF ELSE BEGIN
		tmp=FXPAR(header,'DATE-OBS')
		tmp=STRSPLIT(tmp,'/',/EXTRACT)
		date='19'+STRN(tmp[2])+'-'+STRN(tmp[1])+'-'+STRN(tmp[0])
	ENDELSE

        name=obj+'-'+date+'-'+filter+'-CFHT-'+number+'.fits'
       
        MWRFITS,ims[*,*,j],name,header,/create
     ENDFOR
END
