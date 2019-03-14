pro combine,list

  n=0L
  list=FILELIST(list,n)
  names=STRARR(n)
  FOR i=0L,n-1 DO BEGIN
     tmp=STRSPLIT(list[i],'-',/EXTRACT)
     n_tmp=N_ELEMENTS(tmp)
     IF tmp[0] eq 'Radial' THEN BEGIN
        names[i]=tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[n_tmp-3]+tmp[n_tmp-2]
     ENDIF ELSE BEGIN
        names[i]=tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[n_tmp-3]+tmp[n_tmp-2]
     ENDELSE
  ENDFOR
  uname=names[UNIQ(names)]
  nname=N_ELEMENTS(uname)

  FOR i=0,nname-1 DO BEGIN
     print,i,nname,'  ',uname[i]
     newlist=list[WHERE(names eq uname[i])]
     nims=N_ELEMENTS(newlist)
     starcoord=FLTARR(2,nims)
     imsize=FLTARR(2,nims)
     FOR j=0,nims-1 DO BEGIN
        header=HEADFITS(newlist[j])
        starcoord[0,j]=FXPAR(header,'CRPIX1A')
        starcoord[1,j]=FXPAR(header,'CRPIX2A')
        imsize[0,j]=FLOAT(FXPAR(header,'NAXIS1'))
        imsize[1,j]=FLOAT(FXPAR(header,'NAXIS2'))
     ENDFOR


     xpada=CEIL(MAX(starcoord[0,*]))
     xpadb=CEIL(MAX(imsize[0,*]-starcoord[0,*]))
     ypada=CEIL(MAX(starcoord[1,*]))
     ypadb=CEIL(MAX(imsize[1,*]-starcoord[1,*]))
     
     xsize=xpada+xpadb
     ysize=ypada+ypadb
     xcen=xpada+0.5
     ycen=ypada+0.5

     ims=FLTARR(xsize,ysize,nims)
     count=FLTARR(xsize,ysize)
     ims[*,*,*]=!Values.F_NaN
     FOR j=0,nims-1 DO BEGIN
        header=HEADFITS(newlist[j])
        xshift=(xcen-FXPAR(header,'CRPIX1A'))
        yshift=(ycen-FXPAR(header,'CRPIX2A'))
        ims[0+xshift:imsize[0,j]-1+xshift,0+yshift:imsize[1,j]-1+yshift,j]=MRDFITS(newlist[j],0,header,/FSCALE,/SILENT)
        ind=WHERE(ims[*,*,j] eq ims[*,*,j])
        count[ind]+=1
     ENDFOR

     ;total=TOTAL(ims,3,/NAN)
     ;avg=total/count

     avg=MEDIAN(ims,DIMENSION=3)
     
     FXADDPAR,header,'CRPIX1A',xcen
     FXADDPAR,header,'CRPIX2A',ycen
     tmp=newlist[0]
     ;STRREPLACE,tmp,'-000.fits','-Combined.fits'
     tmp = 'Combined-' + tmp
     STRREPLACE,tmp,'.gz',''
     MWRFITS,avg,tmp,header,/CREATE
  ENDFOR
END
