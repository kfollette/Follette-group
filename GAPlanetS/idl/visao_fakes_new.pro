;+
; NAME: visao_fakes
;
; PURPOSE:
;  By default, will scale down the image of the star and inject it at the specified location(s). For saturated data,
;  injects a Gaussian fake planet with specified fwhm.
;
; INPUTS:
;  imcube  :   the cube of images into which you'd like to inject fake planets
;  rotoffs :   a vector of rotational offsets for these images
;  contrast:   contrast between fake planet peak and stellar peak
;  sep     :   separation, in arcseconds, for the fake planet(s). Can have multiple elements.
;
; INPUT KEYWORDS:
;   saturated  :  two element vector [A,B], where A is the value of stellar peak in raw image counts (should be measured from ghost), and B is the desired FWHM
;   pa         :  position angle (measured east of north) for injected fakes. Can have multiple elements
;   nplanets   :  number of planets to inject at each separation. If not specified, default is 4
;   klipparams : a four element vector of KLIP parameters [A,B,C,D] where A is radial zone size(in pixels), B is azimuthal zone size (in degrees), C is exclusion criterion (in pixels)
;                and D is the minimum radius
;   fixpix     : if feeding in an image with NaNs, interpolate over them before proceeding
; OUTPUTS:
;
; OUTPUT KEYWORDS:
;
;
; EXAMPLE:
;
;
; HISTORY:
;  Written June 2016
;
;-


pro visao_fakes_new, filename, rotoffs, contrast, sep, pa=pa, saturated=saturated, $
  nplanets=nplanets, klipparams=klipparams, fixpix=fixpix

  imcube = readfits(filename+'.fits', head)
  
  dim1=(size(imcube))[1]
  dim2=(size(imcube))[2]

  if keyword_set(fixpix) then begin
    fixpix, imcube, badpix, imcube, /nan, npix=20, /silent
  endif
  ;;; if don't specify a pa or set of pas, then generate a random PA between 0 and 360
  ;;; fill annulus with nplanets (default=4) evenly separated from this random PA
  if not keyword_set(pa) then begin
    seed = systime(1)
    pa1=randomu(seed)*360
    if not keyword_set(nplanets) then begin
      nplanets=4
      PA=dblarr(4)
    endif
    for i=0, nplanets-1 do begin
      PA[i]=pa1+(i+1)*(360/nplanets)
    endfor
  endif

  ;;convert separation to pixels
  platescale=0.00798
  seppix=sep/platescale

  ;;for saturated data, need to inject a gaussian planet with the appropriate brightness
  if keyword_set(saturated) then begin
    fwhm = saturated[1]
    adi_psfmodel, fakes, dim1, dim2, rotoffs+90-0.59, seppix, pa, fwhm=fwhm
    scale = contrast*saturated[0]
  endif else begin
    stop
    adi_psfmodel, fakes, dim1, dim2, rotoffs+90-0.59, seppix, pa, psf0=imcube
    scale = contrast
  endelse

  ;;inject fake planets into image cube
  inim=imcube+scale*fakes

  ;;define KLIP parameters
  if keyword_set(klip) then begin
    rzone=klipparams[0]
    azzone=klipparams[1]
    rotmask=klipparams[2]
    minrad=klipparams[3]
  endif else begin
    ;;apply default klip parameters
    minrad = 10
    rzone=dim1/2-minrad
    azzone=360
    rotmask=1
  endelse

  ;;run KLIP on images with fakes

  fname = filename[:-5] + '_fakes'

  sxaddpar, head_new, 'PAS', pa
  sxaddpar, head_new, 'SEPS', sep
  sxaddpar, head_new, 'CONTRAST', contrast

  writefits, filename+'.fits', inim, head_new
  
  ;pca_regions, finim, inim, rotoffs+50-0.59, rotmask, rzone, azzone, [1,2,3,4,5,10,20,50,100], minrad=minrad, fitsfile=string(fname)+'.fits'

  stop
end
