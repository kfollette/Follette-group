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
;  sep     :   separation, in pixels, for the fake planet(s). Can have multiple elements.
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
;  Written June 2016 by KBF
;  Modified April 2018 by KBF to accommodate multiple contrast values
;
;-


pro visao_fakes_new, filename, rotoff_fname, contrast, sep, pa=pa, saturated=saturated, $
  nplanets=nplanets, klipparams=klipparams, fixpix=fixpix, suffix=suffix

  imcube = readfits(filename+'.fits', head)
  rotoffs = readfits(rotoff_fname+'.fits')

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
  ;platescale=0.00798
  ;seppix=sep/platescale

  ;;for saturated data, need to inject a gaussian planet with the appropriate brightness
  if keyword_set(saturated) then begin
  fwhm = saturated[1]
    if n_elements(contrast) eq 1 then begin
      adi_psfmodel, fakes, dim1, dim2, rotoffs+90-0.59, sep, pa, fwhm=fwhm
      scale = contrast*saturated[0]
      inim=imcube+scale*fakes
    endif else begin
      print, 'multiple contrast loop'
      numplanets = n_elements(contrast)
      fakes=dblarr(dim1,dim2,(size(imcube))[3])
      for i=0, numplanets-1 do begin
        adi_psfmodel, fakepls, dim1, dim2, rotoffs+90-0.59, sep[i], pa[i], fwhm=fwhm
        print, contrast[i]
        fakes=fakes+fakepls*contrast[i]
      endfor
      writefits, 'test_fakes.fits', fakes
      inim=imcube+fakes
    endelse
  endif

  if not keyword_set(saturated) then begin
    if n_elements(contrast) eq 1 then begin
      adi_psfmodel, fakes, dim1, dim2, rotoffs+90-0.59, sep, pa, psf0=imcube
      scale = contrast[0]
      inim=imcube+scale*fakes
      ;writefits, 'test_fakes.fits', fakes
    endif else begin
      numplanets = n_elements(contrast)
      fakes=dblarr(dim1,dim2,(size(imcube))[3])
      for i=0, numplanets-1 do begin
        adi_psfmodel, fakepls, dim1, dim2, rotoffs+90-0.59, sep[i], pa[i], psf0=imcube
        fakes=fakes+fakepls*contrast[i]
      endfor
      ;writefits, 'test_fakes.fits', fakes
      inim=imcube+fakes
    endelse
  endif

  ;;inject fake planets into image cube


  if keyword_set(suffix) then begin
    fname = filename+'_'+suffix
  endif else begin
    fname = filename+'_fakes'
  endelse

  sxaddpar, head, 'PAS', '['+strjoin(string(pa), ',')+']'
  sxaddpar, head, 'SEPS', '['+strjoin(string(sep), ',')+']'
  sxaddpar, head, 'CONTRAST', '['+strjoin(string(contrast), ',')+']'

  writefits, fname+'.fits', inim, head

  ;pca_regions, finim, inim, rotoffs+50-0.59, rotmask, rzone, azzone, [1,2,3,4,5,10,20,50,100], minrad=minrad, fitsfile=string(fname)+'.fits'

end
