;+
; NAME: makepupil
; PURPOSE:
;  Creates a pupil image.
;
; INPUTS:
; eps     ; the central obscuration
; pixels  ; linear size of the image (e.g. 64 for 64x64)
;
; KEYWORDS:
;  none       
;
; OUTPUTS:
;  returns the image
;
; EXAMPLE:
;   P = makepupil(6.5, .29, 512.) 
;
; HISTORY:
;  Written 2009-04-03 by Jared Males, jrmales@email.arizona.edu
;  Updated 2012-11-06 documentation updated. (Jared Males)
;-
function makepupil, eps, pixels

D = pixels

Z = dblarr(pixels,pixels)


midpt = 0.5*(pixels-1)
scale = D/pixels

for i=0, pixels-1 do begin
    ;as(i) = (i - midpt)*scale
    x1 = (i - midpt)*scale ;as(i)
    for j=0, pixels-1 do begin
        y1 = (j - midpt)*scale
        r = sqrt((x1)^2 + (y1)^2)
        if(r lt D/2 AND r ge eps*D/2) then Z(j,i) = 1
    endfor
endfor

return, Z

end
