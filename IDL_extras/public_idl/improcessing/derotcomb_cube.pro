pro derotcomb_cube, finim, psfsub, derot, xreg=xreg, meancomb=meancomb, sigma=sigma, sigmed=sigmed, relrot=relrot, silent=silent, mask=mask
;+
; NAME: derotcomb_cube
; 
; PURPOSE:
;   De-rotate a cube of images and combine them.
;
; DESCRIPTION:
;   De-rotates the images in psfsub (usually psf subtracted images) by derot.  
;   Then combines them by either the median(default) or mean. NOTE: psfsub is modified, on return it contains
;   the de-rotated images.
;
; INPUTS:
;   psfsub   :  cube of images, format [dim1, dim2, no_images] - modified by this routine.
;   derot   :  vector of angles by which to derotate (specified CCW).
;
; OUTPUTS:
;   finim    : the de-rotated and combined image
;
; KEYWORDS:
;   meancomb     : if set, the images are mean combined.  The default is median.
;   sigma        : if a value is passed here, the sigma clipped mean is used.
;   sigmed       : if set, then the median of the sigma clipped pixels is used
;
; MODIFICATION HISTORY:
;  Written 2013/02/05 by Jared Males (jrmales@email.arizona.edu)
;  2013/07/31: generalized (used to be visao specific).  JRM
;
; BUGS/WISH LIST:
;  None. This is perfect code.
;
;-

get_cubedims, psfsub, dim1, dim2, nims

if(n_elements(sigma) eq 0) then sigma = 0

if(n_elements(derot) eq nims) then begin
   if ~keyword_set(silent) then print, 'derotcomb_cube: De-rotating'
   rot_cube, psfsub, derot, mask=mask
endif


if(keyword_set(xreg)) then begin
   visao_reg_cube, psfsub, /doshift
endif

if(~keyword_set(meancomb) and sigma eq 0) then begin
   if ~keyword_set(silent) then print, 'derotcomb_cube: Forming median image'
   finim = median(psfsub, dim=3) 
endif else begin

   if(sigma gt 0) then begin
      if ~keyword_set(silent) then print, 'derotcomb_cube: Forming sigma-clipped mean image'
      
      finim = fltarr(dim1, dim2)
      for i=0,dim1-1 do begin
         for j=0,dim2-1 do begin
            idx = where(finite(psfsub[i,j,*]))
            if(idx[0] gt -1) then begin
               meanclip, psfsub[i,j,idx], pixmean, clipsig=sigma, subs=subs
               if(keyword_set(sigmed)) then begin
                  pixmean = median((psfsub[i,j,idx])[subs])
               endif
            endif else begin
               pixmean = float('nan')
            endelse
            finim[i,j] = pixmean
         endfor
      endfor
   endif else begin
      if ~keyword_set(silent) then print, 'derotcomb_cube: Forming average image'
      nims = (size(psfsub))[3]
      finim = total(psfsub, 3, /nan)
      finim = finim/nims
   endelse
endelse

end
