  ;+
  ; NAME:rad_ellipse
  ;
  ; PURPOSE: fit an ellipse to a disk cavity by creating radial profiles in azimuthal bins and selecting those for which a gaussian fit is good
  ;
  ; INPUTS:     xcen, ycen is exact center of object in pixel coordinates
  ;             bin = size of radial bins
  ;             n_cuts = number of radial bins
  ;             pixscale = size of pixels in arcsec
  ;             
  ;
  ; INPUT KEYWORDS:
  ;             angbin = size of azimuthal bins for PA cuts, default is 5 degrees
  ;             mask = radius of pixel mask
  ;             binht = vertical height of each bin from midplane (total bin height = binht*2). Default is 5 pixels (total height =10 pixels)
  ;
  ; OUTPUTS:
  ;             file 'ellipse_fit.fits', which contains the peak locations in [x,y] coordinates and the parameters of the fit in the header
  ;
  ; OUTPUT KEYWORDS:
  ;
  ; EXAMPLE:
  ;
  ; HISTORY:
  ;-
  ;
pro rad_ellipse, infile, xcen, ycen, bin, n_cuts, pixscale, angbin=angbin, mask=mask, fname=fname, binht=binht

  if keyword_set(fname) then filename=fname else fname='ellipsefit_test'
  
  file=readfits(string(infile))
  npix=bin*n_cuts*2
  ;create blank arrays
  if keyword_set(angbin) then angbin=angbin else angbin=5
  PAs=indgen(360./angbin)*angbin
  cut=dblarr(npix+1, npix+1, n_cuts, n_elements(PAs))
  cut2=dblarr(npix+1, npix+1, n_cuts, n_elements(PAs))
  r=dblarr(npix+1,npix+1)
  pix=dblarr(npix+1,npix+1,n_cuts, n_elements(PAs))
  file_rot=dblarr((size(file))[1],(size(file))[2],n_elements(PAs))

  if keyword_set(xmin) then xmin=xmin else xmin=0
  if keyword_set(binht) then binht=binht else binht=5
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
  ;  writefits, 'rtest.fits', r

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
          ;;cut annulus to just +/- binht pixels vertically (to get just one PA, not entire annulus)
          if v le npix/2.+binht and v ge npix/2.-binht and h le npix/2 then begin
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

  ; writefits, 'cut.fits', cut
  ;writefits, 'cut2.fits', cut2
  ;writefits, 'pix.fits', pix

  ;;make radial profile
  count=fltarr(n_cuts, n_elements(PAs))
  radprof=fltarr(n_cuts, n_elements(PAs))

  for i=0, n_elements(PAs)-1 do begin
    for w=0, n_cuts-1 do begin
      count[w,i]=total(pix[*,*,w,i],/NAN)
      radprof[w,i]=total(cut2[*,*,w,i],/NAN)/count[w,i]
    endfor
  endfor

  ;; eliminate bins under the mask
  if keyword_set(mask) then begin
    rmin=ceil(mask/bin)-1
    radprof[0:rmin,*]=-20
  endif

  ; write out radial profiles to file
  writefits, 'ellipse_fit_radprofs.fits', radprof

  ;make x axis in arcsec
  xaxis=fltarr(n_cuts)
  xaxis=(indgen(n_cuts))*pixscale*bin+pixscale*bin/2
  xpix=indgen(n_cuts)*bin+bin/2

  ;;write blank arrays for fitting
  loc=dblarr(2,n_elements(PAs))
  rad=dblarr(n_elements(PAs))
  PA2=PAs+180
  coeff=dblarr(n_elements(PAs), 6)
  peak=dblarr(n_elements(PAs))
  chisq=dblarr(n_elements(PAs))
  fit=dblarr(n_cuts, n_elements(PAs))
  x=''

  ;; guess at gaussian fit parameters ---- REMOVE ME LATER
  guess=[20, 60, 15, 20, 0.1, 0]

  set_plot, 'x'

  for i=0, n_elements(PAs)-1 do begin
    ;; first guess at peak of radial profile
    ;peaks[i]=where(radprof[*,i] eq max(radprof[*,i]))
    ;; fit profile with gaussian
    fit[rmin+1:n_cuts-1,i]=gaussfit(xpix[rmin+1:n_cuts-1],radprof[rmin+1:n_cuts-1,i],a, chisq=b, estimates=guess)
    coeff[i,*]=a
    chisq[i]=b

    ;; show the profile and fit to user
    plot, radprof[*,i], psym=4, xtitle='radial bin number', ytitle= 'Intensity'
    xyouts, 100, 100, 'PA='+string(PAs[i]+90), color=cgcolor('white'), /device
    oplot, fit[*,i], color=250
    oplot, [coeff[i,1]-1,coeff[i,1]-1]/bin, [-100,100], linestyle=2, color=250

    ;; user specifies which are good
    read, x, prompt='is this a good fit? (y/n)   '

    ;; write good fits into array
    if x eq 'n' then coeff[i,*]=0
    if coeff[i,0] gt 0 then rad[i]=coeff[i,1] 
    if coeff[i,0] gt 0 then peak[i]=coeff[i,0]

    ;; compute x,y pixel location for peak of fit at each PA
    loc[*,i]=[rad[i]*cos(PA2[i]*!PI/180), rad[i]*sin(PA2[i]*!PI/180)]
  endfor

  ;;empty zeroes out of array
  xpts=loc[0,*]
  ypts=loc[1,*]
  xkeep=xpts[where(abs(rad) gt 0)]
  ykeep=ypts[where(abs(rad) gt 0)]
  peakkeep=peak[where(abs(rad) gt 0)]
STOP
  ;;fit with ellipse
  fit=mpfitellipse(xkeep, ykeep, /tilt, weights=peakkeep/max(peakkeep), niter=10)

  ;;plot all good points and fit
  plot, xkeep, ykeep, psym=2
  tvellipse, fit[0], fit[1], fit[2], fit[3], fit[4], color=cgcolor('red')
  z=''

  ;;user specifies whether to write into fits array
  read, z, prompt='shall I write this fit to ellipse_fit.fits in current directory (y/n)   '
  if z eq 'y' then begin

    ;;write points into fits file
    pts=[[xkeep],[ykeep]]
    writefits, string(fname)+'.fits', pts

    ;;write ellipse fit parameters into header
    pt_check=readfits(string(fname)+'.fits',head)

    sxaddpar, head, 'major axis=', string(fit[0])
    sxaddpar, head, 'minoraxis=', string(fit[1])
    sxaddpar, head, 'x center=', string(fit[2])
    sxaddpar, head, 'y center=', string(fit[3])
    sxaddpar, head, 'tilt=', string(fit[4]*180/!PI)
  endif

  stop
end
