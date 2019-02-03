pro airy_slit, sz, lodmax, off, sw_as, eps, lam, f, A, x0, y0
; NAME:
;    airy_slitloss
;
; Purpose:
;    Calculates the slitloss for a point source geometrically.
;
; Inputs:
;    sz     - the size of image to use, in pixels (i.e. and array of sz x sz)
;    lodmax - the maximum lambda/D to consider, a vector
;    sw_as  - the slit width, in arc seconds
;    ar_as  - the photometric aperature radius, in arc seconds
;    eps    - the central obscuration fraction
;    lam    - the input wavelengths
;
; Outputs:
;    enc - the fraction of flux enclosed by the slit 
;
; Keywords:
;    
; 
; MODIFICATION HISTORY:
;       Author: Jared R. Males
;       Date: May 29, 2010
;
r0 = rarr(sz, q, x0, y0)
f = dblarr(sz)

r = r0*lodmax
y = (y0+0.*off/sz)*lodmax
x = (x0+2.*off/sz)*lodmax
r = sqrt(x^2+y^2)

A = airy(r, eps)


sw = (.5*sw_as)/(.2063*lam/6.5)



for i =0, sz-1 do begin
   idx = where(abs(y[i,*]) le sw)
   f[i] = total(A[i,idx])

endfor

end

	