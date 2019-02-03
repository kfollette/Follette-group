;+
; NAME: angsub
;
; PURPOSE: 
;  Properly subtract angles across 360 degrees
;
; Description:
;  Calculates the difference between two angles even when crossing 0 or 360 degrees.  The resulting difference is always 
;  less than 180 degrees.
;
; INPUTS:
;  q1  :  the first angle, by default in degrees
;  q2  :  the second angle, by default in degrees
;
; INPUT KEYWORDS:
;  radians   :  set if inputs are in radians.  output will be in radians.
;
; OUTPUTS:
;    returns q1-q2, guarenteed to be less than 180 degrees (or pi/2 if in radians)
;
; EXAMPLE:
;  dq = angsub(355, 5)
;       this returns dq = -10 degrees (not 350!)
;
; HISTORY:
;  2012-02-28: Written by Jared Males, jrmales@email.arizona.edu
;  2012-03-16: Jared Males added radians
;-
function angsub, q1, q2, radians=radians

twopi = 360.d
halfpi = (180.d)

;Do this without the slow keyword_set call
if(n_elements(radians) gt 0) then begin
   if(radians) then begin
      twopi = 2.d*!dpi
      halfpi = (0.5d*!dpi)
   endif
endif

dq = (q1-q2)

adq = abs(dq)
idx = where(adq gt halfpi)
if(idx[0] ne -1) then dq[idx] = dq[idx] - twopi*sign(dq[idx])

return, dq

end

