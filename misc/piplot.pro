;+
; NAME:
;
; PURPOSE:
;
; INPUTS:
;
; INPUT KEYWORDS:
;
; OUTPUTS:
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
;
; HISTORY:
;-

pro PIplot, objname, expt, type, pix, fmin=fmin, fmax=fmax, nlevels=nlevels, imsz=imsz, mask=mask, fname=fname, ptitle=ptitle, $
  dist=dist, PAs=PAs, ellip=ellip, star=star, points=points, ct=ct, flux=flux, ndivs=ndivs, norm=norm


  image=readfits(string(objname)+'_'+strcompress(string(expt),/remove_all)+'_'+string(type)+'.fits', junkheader, /NO_UNSIGNED)

  if keyword_set(fname) then filename=fname else fname='test'
  if keyword_set(ptitle) then ttl=ptitle else ttl=''
  ndim1=(size(image))[1]
  ndim2=(size(image))[2]
  if keyword_set(ct) then loadct, ct else loadct, 4

  if not keyword_set(imsz) then imsz=ndim1

  ;; change coords to arcsec from center
  xcoord=-(indgen(ndim1)*pix-pix*(ndim1-1)/2.)
  ycoord=indgen(ndim2)*pix-pix*(ndim2-1)/2.
  xmin=-(imsz/2)*pix
  xmax=(imsz/2)*pix
  ymin=xmin
  ymax=xmax

  ;;add a second x axis in AU
  if keyword_set(dist) then begin
    xaxis2=tan(xcoord/206265.)*dist*206265.
    xmin2=tan(xmin/206265.)*dist*206265.
    xmax2=abs(xmin2)
  endif


  ;if not keyword_set(fmin) then fmin=min(image)*1.
  if not keyword_set(fmax) then fmax=max(image)*1.1
  if not keyword_set(nlevels) then nlevels=20
  if not keyword_set(ndivs) then ndivs=10

  if keyword_set(mask) then mask=mask*pix else mask=0

  if keyword_set(bgd) then bgd=bgd else bgd=0

  levels=(fmax-fmin)/(nlevels-1.)*indgen(nlevels)+fmin
  levels[0]=min(image)   ;catch all below floor in first bin so no white
  colors=250/(nlevels-1)*indgen(nlevels)

  if keyword_set(dist) then xsty=9 else xsty=1

  loadct, 4

  set_plot, 'ps'
  device, filename=string(fname)+'.eps', /portrait, bits_per_pixel=8, xsize=17, ysize=15.5, /encapsulated
  !P.MULTI=[0,0,1]
  contour, image, xcoord, ycoord, xrange=[xmax,xmin],yrange=[xmin,xmax], title=ttl, xstyle=xsty, /ystyle,  position=[0.1,0.1,0.8,0.95],levels=levels, $
    c_colors=colors, xtitle='RA offset (in arcsec)', ytitle='Dec offset (in arcsec)', /isotropic, /fill, thick=3, charthick=3, $
    xthick=3, ythick=3

  if keyword_set(dist) then begin
    axis, xaxis=1, xrange=[xmin2,xmax2], xtitle='RA offset (in AU)', xtickinterval=50, charthick=3, xthick=3, subtitle=ttl
  endif

  if keyword_set(PAs) then begin
    print, 'yeah'
    x1=(ndim1-1)/2.
    y1=(ndim2-1)/2.
    x2=dblarr(n_elements(PAs))
    y2=dblarr(n_elements(PAs))
    PAidl=dblarr(n_elements(PAs))
    vec_colors=250/(n_elements(PAs)-1)*indgen(n_elements(PAs))
    for n=0, n_elements(PAs)-1 do begin
      if PAs[n] gt 90 then PAidl[n]=Pas[n]-270
      x2[n]=sin(PAidl[n]*!PI/180)*imsz/2
      y2[n]=cos(PAidl[n]*!PI/180)*imsz/2
      ;arrow,
    endfor
  endif

  loadct, 4
  if keyword_set(flux) then begin
    fsc_colorbar, ncolors=255, /vertical, position=[0.82,0.1,0.87,0.87], divisions=ndivs,$
      maxrange=fmax, minrange=fmin, title='Intensity in mJy/arcsec!U2!N',/right, $
      charthick=3, Format='(f4.1)'
  endif else begin

    if keyword_set(norm) then begin
      fsc_colorbar, ncolors=255, /vertical, position=[0.82,0.1,0.87,0.87], divisions=ndivs, $
        maxrange=fmax, minrange=fmin, title='Normalized Intensity',/right, charthick=3, Format='(f3.1)'
    endif else begin
      fsc_colorbar, ncolors=255, /vertical, position=[0.82,0.1,0.87,0.87], divisions=ndivs, $
        maxrange=fmax, minrange=fmin, title='Intensity in Counts',/right, charthick=3
    endelse
    
  endelse

  tvcircle, mask, 0, 0, /data, /fill   ;software mask
  if keyword_set(star) then begin
    plots, star[0]*pix, star[1]*pix, psym=2, color=cgcolor('white'), thick=4, symsize=2
  endif

  if keyword_set(ellip) then begin
    a=ellip[0]
    b=ellip[1]
    tilt=ellip[2]
    xcen=ellip[3]
    ycen=ellip[4]
    tvellipse, a*pix, b*pix, xcen*pix, ycen*pix, tilt, color=cgcolor('white'), /data, linestyle=2, thick=4
    plots, xcen*pix, ycen*pix, psym=1, color=cgcolor('white'), thick=4, symsize=2
  endif

  if keyword_set(points) then begin
    z=readfits(points)
    oplot, -1.*z[*,0]*pix, z[*,1]*pix, psym=1, color=cgcolor('white'), thick=8, symsize=1.5
  endif

  device, /close_file
  stop
end