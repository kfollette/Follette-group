function frac_pixel, r, x1, y1, x2, y2, w=w
; (x1,y1) and (x2,y2) are the normalized (0..1) coordinates of the intercepts of the circle on the pixel
; where the origin of the coordinates is the corner closest to the center of the circle.

if(n_elements(w) eq 2) then return, 0. ;a corner

if(n_elements(w) gt 4) then return, 1. ;nearly inscribed

if(n_elements(w) eq 4) then begin
   x1 = w[0]
   y1 = w[1]
   x2 = w[2]
   y2 = w[3]
endif

h = 0.5 * sqrt((x2-x1)^2 + (y2-y1)^2)

sinq = h/r
q = asin(sinq)

Aseg = (q*r^2 - h*sqrt(r^2-h^2))
;print, Aseg
Apix = -Aseg - 1.d

if( (y1 eq 0  and x2 eq 0)  or  (y2 eq 0 and x1 eq 0) ) then begin
   ;Lower corner
   ;single triangle
   h1 = max([y1,y2])
   ;if(y1 eq 0) then h1 = y2
   b1 = max([x1, x2])
   ;if(x1 eq 0) then b1 = x2
   
   Apix = 0.5*h1*b1
   
endif

if ( (x1 eq 1 and y2 eq 1) or (x2 eq 1 and y1 eq 1) ) then begin
   ;Upper corner
   h1 = abs(y2-y1) 
   b1 = abs(x2-x1)
   s1 = min([x2,x1])
   s2 = min([y2,y1])
   ;print, h1, b1, s1, s2
   Apix = 0.5*h1*b1 + b1*s2 + h1*s1 + s1*s2
endif


if ( (x1 eq 0 and x2 eq 1) or (x2 eq 0 and x1 eq 1) ) then begin
   ;crossing-y
   b1 = abs(y2-y1)
   s1 = min([y2,y1])
   Apix = 0.5*b1 + s1
   
endif

if ( (y1 eq 0 and y2 eq 1) or (y2 eq 0 and y1 eq 1) ) then begin
   ;crossing-x
   ;print, 'cross x'
   b1 = abs(x2-x1)
   s1 = min([x2,x1])
   Apix = 0.5*b1 + s1
   
endif

if ( (y1 eq y2 and x1 ne 0. and x1 ne 1.) or (x1 eq x2 and y1 ne 0. and y1 ne 1.) ) then begin
   ;cross-one
   ;print, 'cross one'
   Apix = 0.
endif
   
   
return, Apix +  Aseg

end

