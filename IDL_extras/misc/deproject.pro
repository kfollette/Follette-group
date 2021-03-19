pro deproject, infile, xcen, ycen, PA, incl, fname=fname, r2=r2

  ;+
  ; NAME: deproject
  ;
  ; PURPOSE:
  ;  deproject a disk image according to major axis PA and inclination
  ;
  ; INPUTS:
  ;  infile : disk image to deproject
  ;  xcen : center of image in x
  ;  ycen:  center of image in y
  ;  PA :    position angle of major axis
  ;  incl:   disk inclination
  ;
  ; INPUT KEYWORDS:
  ;  fname  : output filename
  ;  r2     : scale by r^2
  ;
  ; OUTPUTS:
  ;
  ; OUTPUT KEYWORDS:
  ;    none
  ;
  ; EXAMPLE:
  ;
  ;
  ; HISTORY:
  ;  Written 2014 sometime by Kate Follette, kbf@stanford.edu
  ;  modified 05-05-2016 to include r^2 scaling
  ;
  ;-

  image=infile
  if keyword_set(fname) then filename=fname else fname='test'

  ;rotate so that PA of major axis is directly L/R
  imrot=rot(image, PA-90,1,xcen,ycen,/pivot, cubic=-0.5)

  ndim1=(size(image))[1]
  ndim2=(size(image))[2]

  mag=1./cos(incl*!PI/180)

  ;;stretch image in y by cos(incl) to deproject and divide by that factor to preserve flux
  im_rebin=congrid(imrot, ndim1, ndim2*mag, cubic=-0.5)/mag

  ndim2_after=(size(im_rebin))[2]
  ycen2=ndim2_after/2

  if keyword_set(r2) then begin
    xdim = (size(infile))[1]
    ydim = (size(infile))[2]
    r=dblarr(xdim, ydim)
    for x=0, xdim-1 do begin
      for y=0, ydim-1 do begin
        r[x,y]=sqrt((x-xcen)^2+(y-ycen)^2)
        im_rebin[x,y]=im_rebin[x,y]*(r[x,y])^2
      endfor
    endfor
  endif

  ;;rotate back to original orientation
  im_rebin_rot=rot(im_rebin, -(PA-90), 1, xcen, ycen2, /pivot, cubic=-0.5)


  writefits, string(fname)+'.fits', im_rebin_rot



  stop
end
