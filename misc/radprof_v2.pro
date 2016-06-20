;+
; NAME: RADPROF_V2
;
; PURPOSE: Plot radial profiles at arbitrary PAs around a specified pixel
;
; INPUTS:
;           xcen, ycen is exact center of object in pixel coordinates
;           bin = size of radial bins for profile
;           n_cuts = number of radial bins
;           PAs = array of PA values for cuts
;           pixscale = size of pixels in arcsec
;           lgd = array of labels for each PA cut ex:['major axis', 'minor axis', 'PA=135']
;
; INPUT KEYWORDS:
;           ttl= title of plot
;           fname = name of output file (will be appended with .ps). Default is test.ps.
;           xmin = minimum r in arcsec to begin plot (exclude inner region)
;           vline= ovrplot vertical line at specified number of AU
; OUTPUTS:
;
; OUTPUT KEYWORDS:
;
; EXAMPLE:
;
; HISTORY:
;-

pro radprof_v2, infile, xcen, ycen, bin, n_cuts, PAs, pixscale, lgd, ymax=ymax, dist=dist, $
  xmin=xmin, ttl=ttl, fname=fname, linesty=linesty, irs48=irs48, flux=flux, $
  vline=vline, linelab=linelab, ymin=ymin



  if keyword_set(fname) then filename=fname else fname='test'

  file=readfits(string(infile))
  npix=bin*n_cuts*2
  ;create blank arrays
  cut=dblarr(npix+1, npix+1, n_cuts, n_elements(PAs))
  cut2=dblarr(npix+1, npix+1, n_cuts, n_elements(PAs))
  r=dblarr(npix+1,npix+1)
  pix=dblarr(npix+1,npix+1,n_cuts, n_elements(PAs))
  file_rot=dblarr((size(file))[1],(size(file))[2],n_elements(PAs))

  if keyword_set(xmin) then xmin=xmin else xmin=0
  if keyword_set(ttl) then ttl=ttl else ttl=''
  if keyword_set(linesty) then linesty=linesty else linesty=fltarr(n_elements(PAs))

  ;loop through PAs and rotate image so that that PA value is directly L/R
  for i= 0, n_elements(PAs)-1 do begin
    PAs[i]=PAs[i]-90
    file_rot[*,*,i]=rot(file,PAs[i],1,xcen,ycen,/pivot, cubic=-0.5)
  endfor
  writefits, 'rot.fits', file_rot

  ;create array with distance to each pixel from center
  for x=xcen-(npix)/2, xcen+(npix)/2 do begin
    for y=ycen-(npix)/2, ycen+(npix)/2 do begin
      h=x-xcen+npix/2.
      v=y-ycen+npix/2.
      r[h,v]=sqrt((h-npix/2.)^2+(v-npix/2.)^2)
    endfor
  endfor
  ;writefits, 'rtest.fits', r


  ;make annuli for each radial bin

  for h=0, npix do begin
    for v=0, npix do begin
      for w=0, n_cuts-1 do begin
        for i=0, n_elements(PAs)-1 do begin
          if (r[h,v] lt bin*(w+1)) and (r[h,v] ge bin*w) then begin
            cut[h,v,w,i]=file_rot[h+xcen-npix/2., v+ycen-npix/2.,i]
          endif else begin
            cut[h,v,w,i]=0.
          endelse
          ;;cut annulus to just +/- 5 pixels vertically (to get just one PA, not entire annulus)
          if v le npix/2.+5 and v ge npix/2.-5 and h le npix/2 then begin
            cut2[h,v,w,i]=cut[h,v,w,i]
          endif
          ;count number of pixels in each bin
          if cut2[h,v,w,i] ne 0 then begin
            pix[h,v,w,i]=1
          endif
        endfor
      endfor
    endfor
  endfor

  writefits, 'cut.fits', cut
  writefits, 'cut2.fits', cut2
  writefits, 'pix.fits', pix

  ;;make radial profile
  count=fltarr(n_cuts, n_elements(PAs))
  radprof=fltarr(n_cuts, n_elements(PAs))


  for i=0, n_elements(PAs)-1 do begin
    for w=0, n_cuts-1 do begin
      count[w,i]=total(pix[*,*,w,i],/NAN)
      radprof[w,i]=total(cut2[*,*,w,i],/NAN)/count[w,i]
    endfor
  endfor

  ;make x axis in arcsec
  xaxis=fltarr(101)
  xaxis=(indgen(n_cuts))*pixscale*bin+pixscale*bin/2

  ;;add a second x axis in AU
  if keyword_set(dist) then begin
    xaxis2=tan(xaxis/206265.)*dist*206265.
    if keyword_set(xmin) then xmin2=tan(xmin/206265.)*dist*206265. else xmin2=0
    xmax2=tan(pixscale*n_cuts*bin/206265.)*dist*206265.
  endif

  syms=[6,5,2,1,4,7,9,11,33]
  colors=dblarr(n_elements(PAs))
  for i=0, n_elements(PAs)-1 do begin
    colors[i]=250/(n_elements(PAs)-1)*i
  endfor
  
  if keyword_set(irs48) then begin
    colors= [0, 0, 50, 50, 250, 250]
  endif

  if keyword_set(ymax) then begin
    ymaxm=ymax
  endif else begin
    ymaxm=max(radprof)*1.1
  endelse

  if keyword_set(dist) then xsty=9 else xsty=1

  if keyword_set(flux) then yttl='Intensity in mJy/arcsec!U2!N' else yttl='Normalized Intensity'
  if not keyword_set(ymin) then ymin=0

  loadct, 39
  set_plot, 'ps'
  device, filename=string(fname)+'.eps', /color, bits_per_pixel=8, xsize=17, ysize=17
  !P.MULTI=0
  plot, (xaxis), (radprof[*,0]), psym=6, xtitle='Distance in arcseconds', ytitle=string(yttl), title=string(ttl), $
    xrange=[xmin,pixscale*n_cuts*bin], yrange=[ymin,ymaxm], xstyle=xsty, /ystyle, thick=3, charthick=3, $
    xthick=3, ythick=3, position=[0.18,0.15,0.9,0.85], charsize=2
  oplot, xaxis, radprof[*,0], thick=3, linestyle=linesty[0], color=colors[0]
  for i=1, n_elements(PAs)-1 do begin
    print, 'overplotting', i, colors[i], syms[i]
    oplot, xaxis, radprof[*,i], thick=3, color=colors[i], linestyle=linesty[i]
    oplot, xaxis, radprof[*,i], thick=3, psym=symcat(syms[i]), color=colors[i]
  endfor
  if keyword_set(dist) then begin
    axis, xaxis=1, xrange=[xmin2,xmax2], xtitle='Distance in AU', xtickinterval=50, charthick=3, charsize=2, xthick=3, subtitle=ttl
  endif
  if keyword_set(vline) then begin
    oplot, [vline, vline], [0, ymax], linestyle=2, thick=3
    xyouts, vline*0.95, (ymax)/2, string(linelab), orientation=90
  endif
  legend, lgd, color=colors, psym=syms[0:n_elements(PAs)-1], box=0, /right, /top, thick=3, charthick=3, charsize=1
  device, /CLOSE

  stop

end
