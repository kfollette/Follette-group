pro adi_psfmodel, model, dim1, dim2, derot, sep, pa, psf0=psf0, fwhm=fwhm
;+
; NAME: adi_psfmodel
;
; PURPOSE: 
;  Create a cube of images which simulates an ADI observation.
;
; Description:
;  Creates a blank cube of size dim1 X dim2 X n_elements(derot), and injects a PSF 
;  model at the appropriate position based on the derot input.  The psf model is either
;  psf0 or can be a gaussian with fwhm.  If psf0 is a cube of length n_elements(derot)
;  then the individual PSF snapshots are used as the PSF at each step in the cube --
;  use this if you have an unsaturated core and want to model time-variability etc.
;  No scaling is applied.
;
; INPUTS:
;  dim1    :  x dimension of output cube    
;  dim2    :  y  dimension of output cube
;  derot  :  the rotation angle of each step in the cube
;  sep     :  the separation of the model PSF to inject (can be a vector for multiple copies)
;  pa      :  the position angle where to inject the PSF (can be a vector)
;
; OUTPUTS:
;  model   :  the resultant model
;
; INPUT KEYWORDS:
;  psf0    :  the PSF to inject, can be a cube to model a time-variable PSF
;  fwhm    :  the FWHM of a simple Gaussian to use instead of PSF0
;  
; HISTORY:
;  2013-ish:  created by Jared R. Males (jrmales@email.arizona.edu)
;  2014-07-28: documentation updated by JRM
;-

nims = n_elements(derot)
npsfs = n_elements(sep)

model = fltarr(dim1, dim2, nims)


xcen = .5*(dim1-1.)
ycen = .5*(dim2-1.)

if(n_elements(psf0) gt 0) then begin
   get_cubedims, psf0, psfd1, psfd2, psf_nims
   
   padpsf = pad_ims(psf0, 0.5*(dim1-psfd1))
   
endif

for i=0, nims-1 do begin

   ang = (pa) - derot[i]

   for j=0, npsfs-1 do begin
   
      if(n_elements(psf0) gt 0) then begin
         x = xcen - sep[j]*sin(-ang[j]*!dtor)
         y = ycen - sep[j]*cos(-ang[j]*!dtor)

         if(psf_nims gt 1) then begin
            psf = rot(padpsf[*,*,i], 0, 1, x, y, cubic=0.5)
         endif else begin
            psf = rot(padpsf, 0, 1, x, y, cubic=0.5)
         endelse
         
      endif else begin
         x = xcen + sep[j]*sin(-ang[j]*!dtor)
         y = ycen + sep[j]*cos(-ang[j]*!dtor)

         psf = psf_gaussian(npixel=[dim1, dim2], centroid=[x,y], fwhm=fwhm)
      endelse
      
      model[*,*, i]= model[*,*, i] + psf

   endfor
endfor

end

