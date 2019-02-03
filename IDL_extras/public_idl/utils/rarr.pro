;+
; NAME: rarr
;
; Description:
;  Create a 2D array, dim1xdim2, with value at each point being the radius of that pixel from the center.  
;  The default center is x=0.5*(dim1-1), y=0.5*(dim2-1), but this can be changed. The radius is specified 
;  relative to dim1/2, that is a pixel a distance dim1/2 from the center will have radius 1.0.  You can
;  have the radius in pixels instead by setting the /pixels keyword.  If desired, will also calculate the 
;  angle at each pixel.
;
; INPUTS:
;  dim1  :  the x dimension
;  dim2  :  the y dimension
;
; INPUT KEYWORDS:
;  xcen   :  the desired x center offset, in units of pixels relative to 0.5*(dim1-1)
;  ycen   :  the desired x center offset, in units of pixels relative to 0.5*(dim1-1)
;  pixels :  output r, x, and y in units of pixels
;
; OUTPUTS:
;    r    :  the return value, the radius array (goes from 0 to 1 along x-axis, unless pixels is set)
;    x    :  the x coordinate of each pixel  (goes from -1 to +1, unless pixels is set)
;    y    :  the y coordinate of each pixel (same as x, unless rectangular)
;
; OUTPUT KEYWORDS:
;   theta :  if passed a variable name, it is set to the angle in radians of each pixel
;
; EXAMPLE:
;  r = rarr(512,255,x,y,xcen=.25, theta=q)
;       this returns a 512x255 array, with center x = 255.75, y = 127.0

; HISTORY:
;  Written by Jared Males, jrmales@email.arizona.edu
;  2012-12-27: updated and documented by Jared Males
;  2013-03-14: Jared Males added pixels keyword
;-
function rarr, dim1, dim2, x, y, xcen=xcen, ycen=ycen, theta=theta, pixels=pixels

;create x and y

x  = (DINDGEN( dim1 ) + 0.5d) / ( double(dim1) ) * 2.0d - 1.0d
x  = x # REPLICATE( 1.0D, dim2)

y = ((DINDGEN( dim2 ) + 0.5d) / ( double(dim2) ) * 2.0d - 1.0d)*(double(dim2)/double(dim1))
y  = REPLICATE( 1.0D, dim1 ) # y

;change center if desired
if(n_elements(xcen) eq 1) then begin
   xc = 2.d*double(xcen)/double(dim1)
   x = x - xc
endif

if(n_elements(ycen) eq 1) then begin
   yc = 2.d*double(ycen)/double(dim2)
   y = y - yc
endif


;Calculate the radius

r = sqrt( x^2 + y^2 )


;If desired calculate the angles
if(arg_present(theta)) then theta = atan( y, x )


;If desired, convert to pixels
if(keyword_set(pixels)) then begin

   r = 0.5*r*dim1
   x = 0.5*x*dim1
   y = 0.5*y*dim1
   
endif

return, r

end
