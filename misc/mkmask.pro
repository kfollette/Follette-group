;+
; NAME: mkmask
;
; PURPOSE:
;  creates image of specified dimensions with a (central unless otherwise specified) mask of a given radius. By default, Pixels under 
;  the mask are 0 and outside are 1
;
; INPUTS:
;  xdim/ydim : x and y dimensions of output image
;  r    : radius of desired mask 
;  
; INPUT KEYWORDS:
;  cen  :  two element array to center mask at location cen[0], cen[1]
;  reverse     :  sets pixels under mask to 1 and outside to 0
;  fits  : writes mask to .fits file called mask.fits
;
; OUTPUTS:
;  mask : a 2D idl mask array
;  
; OUTPUT KEYWORDS:
;    none
;
; EXAMPLE:
;
;
; HISTORY:
;  Written 2015 by Kate Follette, kbf@stanford.edu
;  modified 2016-03-29 to allow off center mask and reversed options
;  modified 2016-06-07 to add nan option
;
;-

pro mkmask, xdim, ydim, r, mask, fits=fits, cen=cen, reverse=reverse, nan=nan


mask=dblarr(xdim, ydim)+1.

if keyword_set(cen) then print, cen[0], cen[1]
for x=0, xdim-1 do begin
  for y=0, ydim-1 do begin
    
    if keyword_set(cen) then begin
       if sqrt((cen[0]-x)^2.+(cen[1]-y)^2.) lt r then mask[x,y]=0.
    endif else begin
          if sqrt(((xdim-1)/2.-x)^2+((ydim-1)/2.-y)^2) lt r then mask[x,y]=0.
    endelse
    
  endfor
endfor


if keyword_set(reverse) then begin
  mask_copy=mask
  mask[where(mask_copy eq 0.)]=1.
  mask[where(mask_copy eq 1.)]=0.
endif

if keyword_set(NaN) then begin
  mask_copy=mask
  mask[where(mask_copy eq 0.)]='NaN'
endif

if keyword_set(fits) then begin
  writefits, 'mask.fits', mask
endif

;return, mask

end
