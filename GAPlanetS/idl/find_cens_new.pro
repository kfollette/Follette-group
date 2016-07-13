;+
; NAME: find_cens
;
; PURPOSE:
;  fits a 2D gaussian to the brightest object in each image in a sequence of sequentially numbered fits images and writes the x and y values, peak value,
;  and x and y FWHM as keywords into each image header
;
; INPUTS:
;  fname :
;
; INPUT KEYWORDS:
;  ghost : also extracts these values for the ghost and writes to header (coordinates of ghost are hard-coded for MagAO H-alpha, but can be modified)
;  writearrs : writes out measured values as .fits files
; OUTPUTS:
;
; OUTPUT KEYWORDS:
;    none
;
; EXAMPLE:
;
;
; HISTORY:
;  Written 2016-03-29 by Kate Follette, kbf@stanford.edu
;
;-

pro find_cens_new, fname, ghost=ghost, fits=fits

  ;;files have to be named sequentially as fname_0001, etc.

  ;;count files that begin with the string fname
  spawn, 'ls -l ' + fname + '* | wc -l', nfile

  ims=readfits(string(fname)+'.fits')

  dim1=(size(ims))[1]
  dim2=(size(ims))[2]
  nims=(size(ims))[3]
  ;  print, nfiles

  ;;write blank arrays for storing measured values
  xcen=dblarr(nims)
  ycen=dblarr(nims)
  peak=dblarr(nims)
  x_fwhm=dblarr(nims)
  y_fwhm=dblarr(nims)

  if keyword_set(ghost) then begin
    ghost_xcen=dblarr(nims)
    ghost_ycen=dblarr(nims)
    ghost_peak=dblarr(nims)
    ghost_x_fwhm=dblarr(nims)
    ghost_y_fwhm=dblarr(nims)
  endif


  ;;if also measuring ghost center and peak, mask out rest of image (except for r=25pix circle around ghost)
  if keyword_set(ghost) then begin
    ;;measured distance from star to ghost
    ghost_offset=[158.,-6.]
    mkmask, dim1, dim2, mask, 25, cen=[(dim1-1)/2.+ghost_offset[0],(dim2-1)/2.+ghost_offset[1]], /reverse
  endif

  ;; loop over images beginning with fname, measuring stellar (and ghost) centroid and peak and writing these into the header
  for i=0, nims-1 do begin
    print, 'image ', i+1, '   of', nims

    im=ims[*,*,i]

    fixpix, im, badpix, outim, /nan, npix=20, /silent

    starfit = mpfit2dpeak(outim, a)
    xcen[i]=a[4]
    ycen[i]=a[5]
    peak[i]=a[1]
    x_fwhm[i]=a[2]*2
    y_fwhm[i]=a[3]*2

    ;;measure ghost as well (if desired)
    if keyword_set(ghost) then begin

      ;;make sure ghost is in the frame
      if dim1 lt 350 then begin
        print, 'ghost not in image frame'
        break
      endif

      ;;multiply by mask to remove star
      mskim=outim*mask

      ghostfit = mpfit2dpeak(mskim, a2)
      ghost_xcen[i]=a2[4]
      ghost_ycen[i]=a2[5]
      ghost_peak[i]=a2[1]
      ghost_x_fwhm[i]=a2[2]*2
      ghost_y_fwhm[i]=a2[3]*2

    endif

  endfor

;;clean outliers (cosmics) - test is if centroid is >2pix from median
peak[where(xcen - median(xcen) gt 2)]='NaN'
xcen[where(xcen - median(xcen) gt 2)]='NaN'
ycen[where(xcen - median(xcen) gt 2)]='NaN'
x_fwhm[where(xcen - median(xcen) gt 2)]='NaN'
y_fwhm[where(xcen - median(xcen) gt 2)]='NaN'

ghost_peak[where(ghost_xcen - median(ghost_xcen) gt 2)]='NaN'
ghost_xcen[where(ghost_xcen - median(ghost_xcen) gt 2)]='NaN'
ghost_ycen[where(ghost_xcen - median(ghost_xcen) gt 2)]='NaN'
ghost_x_fwhm[where(ghost_xcen - median(ghost_xcen) gt 2)]='NaN'
ghost_y_fwhm[where(ghost_xcen - median(ghost_xcen) gt 2)]='NaN'

  if keyword_set(fits) then begin
      
    writefits, fname+'_star_xcen_array.fits', xcen
    writefits, fname+'_star_ycen_array.fits', ycen
    writefits, fname+'_star_peak_array.fits', peak
    writefits, fname+'_star_xfwhm_array.fits', x_fwhm
    writefits, fname+'_star_yfwhm_array.fits', y_fwhm

    if keyword_set(ghost) then begin
      writefits, fname+'_ghost_xcen_array.fits', ghost_xcen
      writefits, fname+'_ghost_ycen_array.fits', ghost_ycen
      writefits, fname+'_ghost_peak_array.fits', ghost_peak
      writefits, fname+'_ghost_xfwhm_array.fits', ghost_x_fwhm
      writefits, fname+'_ghost_yfwhm_array.fits', ghost_y_fwhm
    endif

  endif

  print, 'median stellar x centroid= ', median(xcen), '+/-', stddev(xcen, /nan)
  print, 'median stellar y centroid= ', median(ycen), '+/-', stddev(ycen, /nan)
  print, 'median stellar peak= ', median(peak), '+/-', stddev(peak, /nan)
  print, 'median stellar x FWHM= ', median(x_fwhm), '+/-', stddev(x_fwhm, /nan)
  print, 'median stellar y FWHM= ', median(y_fwhm), '+/-', stddev(y_fwhm, /nan)
  print, 'ratio stellar xFWHM/yFWHM= ', median(x_fwhm)/median(y_fwhm)

  if keyword_set(ghost) then begin
    print, 'median ghost x centroid= ', median(ghost_xcen), '+/-', stddev(ghost_xcen, /nan)
    print, 'median ghost y centroid= ', median(ghost_ycen), '+/-', stddev(ghost_ycen, /nan)
    print, 'median ghost peak= ', median(ghost_peak), '+/-', stddev(ghost_peak, /nan)
    print, 'median ghost x FWHM= ', median(ghost_x_fwhm), '+/-', stddev(ghost_x_fwhm, /nan)
    print, 'median ghost y FWHM= ', median(ghost_y_fwhm), '+/-', stddev(ghost_y_fwhm, /nan)
    print, 'ratio ghost xFWHM/yFWHM= ', median(ghost_x_fwhm)/median(ghost_y_fwhm)
  endif
  
  ;stop

end

