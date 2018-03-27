;+
; NAME: circsym_reg_generic
;
; PURPOSE:
;  separate, register and align all of the raw science images and create an image cube
;
; INPUTS:
;  sci_imlist :
;
; INPUT KEYWORDS:
;  clip     :  crops the array to specified width (square)
;  scl      :  rescale the continuum images by this factor
;  msk      :  radius of saturated region to mask at center of PSF in pixels
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
;  Written 2017-08-15 based heaviy on the existing visao_reg and visao_circsym routines
;
;-

pro circsym_reg_generic, image, clip=clip, fwhm=fwhm, stp=stp, mask=mask

  
  im = readfits(string(image)+'.fits', head)

  print, size(im)

  dim1=(size(im))[1]
  dim2=(size(im))[2]
  nims=(size(im))[3]
  
  xcen=(dim1-1)/2.
  ycen=(dim2-1)/2.

  ;;if clip keyword not defined (not recommended), image size is full array
  if not keyword_set(clip) then clip=dim1

  ;;create blank arrays
  
  im_smooth=dblarr(dim1, dim2, nims)
  im_reg=dblarr(clip, clip, nims)

  peak_ratio=dblarr(nims)
  expt=dblarr(nims)
  im_max_coords=dblarr(2,nims)
  im_shift_arr=dblarr(2,nims)

  ;; PSF for finding star in each image
  if not keyword_set(fwhm) then fwhm=10
  gauss_cen=psf_gaussian(npixel=[dim1, dim2], centroid=[(dim1-1)/2., (dim2-1)/2.], FWHM=fwhm)

  ;;image loop
  for i=0, nims-1 do begin
    
    ;smooth out cosmic rays
    im_smooth[*,*,i]=gauss_smooth(im[*,*,i], width=fwhm)*1.

    ;; calculate approximate offset by registering against gaussian via Fourier cross-correlation, but don't shift 
    ;; note just finding the max after smoothing doesn't work for faint sources
    subreg, gauss_cen, im_smooth[*,*,i], sft, method='F'
    ;stop
    
    ;; do line image first
    ;make 71x71 grid around max coords
    xr=indgen(51.)-50/2.-sft[0]
    ;; line below is 
    yr=indgen(51.)-50/2.-sft[1]
    
    rmax=25 
    
    ;; calculate center of circular symmetry in this region
    if keyword_set(mask) then begin
      print, 'masked' 
      center_circlesym, im[*,*,i], xr, yr, rmax, im_xc, im_yc, im_grid, mask=mask
    endif else begin
      center_circlesym, im[*,*,i], xr, yr, rmax, im_xc, im_yc, im_grid
    endelse

    ;; record shifts as a check
    im_shift=[(dim1-1)/2.-im_xc,(dim2-1)/2.-im_yc]

    ;; nims=(size(im))[3]
    print, 'center of circular symmetry for image is', im_xc, im_yc
    print, 'shifting image', i+1, ' of', nims, ' by', im_shift
    
    ;Line[*,*,j]=shift_interp(Line[*,*,j], Line_shift-1, spline=-0.5) old ver. - based on bug in Fourier cx correlation
    im[*,*,i]=shift_interp(im[*,*,i], im_shift, spline=-0.5)
    
       

    ;;clip may overshoot edge of frame in some images. Adjust clip to be smaller in these cases

    ;;Line adjustments to clip
    if dim1/2-abs(im_shift[0]) lt clip/2 then begin
      print, 'Line overshoots by', dim1/2-abs(im_shift[0])-clip/2, 'pixels in x'
      im_cropx=floor(dim1/2-abs(im_shift[0]))
    endif else begin
      im_cropx=clip/2
    endelse
    if dim2/2-abs(im_shift[1]) lt clip/2 then begin
      print, 'Line overshoots by', dim2/2-abs(im_shift[1])-clip/2, 'pixels in y'
      im_cropy=floor(dim2/2-abs(im_shift[1]))
    endif else begin
      im_cropy=clip/2
    endelse

    ;stop
    im_reg[clip/2-im_cropx:clip/2+im_cropx-1,clip/2-im_cropy:clip/2+im_cropy-1,i]= $
      im[xcen-im_cropx+1:xcen+im_cropx,ycen-im_cropy+1:ycen+im_cropy,i]

    ;;if clip is odd number, then need to shift by 0.5 to get to center of array. VisAO chip channels
    ;are even numbered, so it's expected that the center of the star will be between two pixels
    if clip/2. mod 1 gt 0 then begin
      im_reg[*,*,i]=shift_interp(im_reg[*,*,i], [0.5, 0.5], spline=-0.5)
    endif

    print, 'processed image', i+1, '        of', nims
    ;stop
    
  endfor

  ;;add parameters to header
  sxaddpar, head, 'CLIP', clip
  if not keyword_set(scl) then scl=1
  sxaddpar, head, 'CONT_SCALE', scl 
  if keyword_set(mask) then begin
    sxaddpar, head, 'MASK', mask
  endif
  sxaddpar, head, 'SMOOTH_FWHM', fwhm

 ;; write files
  if keyword_set(clip) then begin
    writefits, string(image)+'_'+string(clip, format='(i03)')+'_circsymreg.fits', im_reg, head
  endif else begin
    writefits, string(image)+'_circsymreg.fits', im_reg, head
  endelse

  if keyword_set(stp) then  stop

  print, 'registration complete'

end
