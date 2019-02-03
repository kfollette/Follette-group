function radon_integral, t, im, xc, yc, alpha
;+
; Calculate the line integral centered at point (xc,yc) along angle alpha
; Just uses nearest neighbor sampling
;
; INPUTS:
;   t   :   vector of line points, e.g. -256 to +56
;  im   :   the image
;  xc   :   the x-coordinate of the line origin, in pixels
;  yc   :   the y-coordinage of the line origin, in pixels
;  alpha : angle of the line, in degrees
;
; OUTPUS:
;  returns a vector of the pixels at each point on the line
;-

get_cubedims, im, dim1, dim2

res = fltarr(n_elements(t))

cos_alph = cos(alpha*!dtor)
sin_alph = sin(alpha*!dtor)

xco = t*cos_alph + xc
yco = t*sin_alph + yc

idx = where(xco ge 0 and xco le dim1-1 and yco ge 0 and yco le dim2-1)

xco = xco[idx]
yco = yco[idx]

res = im[xco, yco]

; for k=0L, n_elements(t)-1 do begin
; 
;    xco = t[k]*cos_alph + xc
;    yco = t[k]*sin_alph + yc
;    
;    if(xco lt 0 or xco gt dim1-1 or yco lt 0 or yco gt dim2-1) then begin
;       res[k] = 0.
;    endif else begin
;       res[k] =  im[xco, yco]
;    endelse
; endfor

return, res

end

