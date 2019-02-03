;+
; NAME: requad_angles
;
; PURPOSE:
;  Adjust a vector of angles such that it is continous if it wraps around 0 or 360.
;
; Description:
;  Adjusts the vector of angles such that it is continous if it wraps around 0 or 360 by moving a portion of it to be lt 0.
;
; INPUTS:
;  angles   :  vector of angles Q such that 0<= Q < 360 (2pi).  Must be monotonic, and have less than 360 (2pi) degree total change.
;
; INPUT KEYWORDS:
;  radians    :  if set, then angles are in radians with range 0<=0<2pi
;
; OUTPUTS:
;  returns the possibly modified vector
;
; HISTORY:
;  Written 2013-01-22 by Jared Males, jrmales@email.arizona.edu
;
;-
function requad_angles, angles, radians=radians

subval = 360.
if(keyword_set(radians)) then subval = 2.*!dpi

;don't modify input
rangles = angles

;find first change
ind = 1
while((ind lt n_elements(angles) - 1) and angles[ind] eq angles[0]) do ind = ind + 1


;Now angles is either increasing, or decreasing.  Choose which part, if any, to make lt 0
if(angles[ind] gt angles[0]) then begin
   
   minang = min(angles)

   
   if(minang lt angles[0]) then begin
   
      idx = where(angles ge angles[0])
            
      if(idx[0] ne -1) then rangles[idx] = rangles[idx] - subval
   endif
   
endif else begin

   maxang = max(angles)
   
   if(maxang gt angles[0]) then begin
      
      idx = where(angles gt angles[0])
      
      if(idx[0] ne -1) then rangles[idx] = rangles[idx] - subval
   endif
   
endelse
         
return, rangles

end

