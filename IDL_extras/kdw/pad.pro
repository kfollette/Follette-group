pro pad,list

  n=0L
  list=FILELIST(list,n)

  FOR i=0,n-1L DO BEGIN
     im=MRDFITS(list[i],0,header,/FSCALE,/SILENT)
     x=FXPAR(header,'CRPIX1A')+0.5
     y=FXPAR(header,'CRPIX2A')+0.5
     print,n-i,' '+list[i],x,y
     
     roundx=ROUND(x)
     roundy=ROUND(y)

     IF roundx-x ge 0.0 THEN newx=roundx-0.5
     IF roundx-x lt 0.0 THEN newx=roundx+0.5
     IF roundy-y ge 0.0 THEN newy=roundy-0.5
     IF roundy-y lt 0.0 THEN newy=roundy+0.5

     xshift=newx-x
     yshift=newy-y

     ;newx=511.5
     ;newy=511.5
     
     P=[-xshift,0,1,0]
     Q=[-yshift,1,0,0]
     im=POLY_2D(im,P,Q,2,CUBIC=-0.5,MISSING=!Values.F_NAN)

     FXADDPAR,header,'CRPIX1A',newx
     FXADDPAR,header,'CRPIX2A',newy
     name=list[i]
     STRREPLACE,name,'.gz',''

     MWRFITS,im,'Shifted/'+name,header,/CREATE
  ENDFOR

END
