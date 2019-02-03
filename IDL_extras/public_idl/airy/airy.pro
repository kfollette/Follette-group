;+
; NAME: airy
; PURPOSE:
;  Calculates the obscured Airy function. 
;
; INPUTS:
;  x  ; location(s) to calculate value of obscured Airy pattern, in lambda/D units
;  eps ; central obscuration fraction, 0<=eps<1
;
; KEYWORDS:
;  none
;
; OUTPUTS:
;  returns the airy pattern
;
; EXAMPLE:
;  airy, myim, .29
;
; HISTORY:
;  Written 2009-09-23 by Jared Males, jrmales@email.arizona.edu
;  2012-10-30. Documentation updated by Jared Males
;
;-
function airy, x, eps

if(n_elements(eps) lt 1) then eps = 0.

y = (1./(1.-eps^2)^2)*(2.*beselj(!pi*x,1)/(!pi*x)-eps*2.*beselj(eps*!pi*x,1)/(!pi*x))^2

;Take care of NaNs
idx = where(finite(y) eq 0)
if(idx[0] ne -1) then y[idx] = 1.

return, y

end
