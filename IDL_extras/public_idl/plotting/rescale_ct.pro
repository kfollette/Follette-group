;+
; NAME: rescale_ct
;
; PURPOSE:
;  Rescales the current color table for a new bias and contrast
;
; DESCRIPTION:
;  Expands or contracts the current color table to have a new bias and contrast.  Same effect as using
;  the right mouse button in ds9.
;
; INPUTS:
;  bias         : the new bias, or midpoint.  0 <= bias <= 1.0
;  contrast     : the new contrast, or relative range.  0 <= contrast
;
; INPUT KEYWORDS:
;  get        : if set, the the new r,g,b values are calculate but not set as the color table
;
; OUTPUTS:
;  newr      : the new red values
;  newg      : the new green values
;  newb      : the new blue values
;
; EXAMPLE:
;  
;  rescale_ct, 0.515, 1.7
; 
;
; HISTORY:
;  Written 2014.09.13 by Jared Males, jrmales@email.arizona.edu
;
;-
pro rescale_ct, bias, contrast, newr, newg, newb, get=get

tvlct, r, g, b, /get

ncol = n_elements(r)

x = findgen(ncol)

centc = (1-bias)*(ncol-1)
print, ':', ncol, centc
minc = centc - .5*(ncol-1)*contrast
maxc = centc + .5*(ncol-1)*contrast

newx = minc + findgen(ncol)*(maxc-minc)/(ncol-1)

linterp, x, r, newx, newr
linterp, x, g, newx, newg
linterp, x, b, newx, newb

;Set bottom of colortable to black
newr[0] = 0
newg[0] = 0
newb[0] = 0

;Now set the new color table
if ~keyword_set(get) then tvlct, newr, newg, newb

end
