pro makecircle, rad, x, y, xc=xc, yc=yc, sample=sample, semiminor=semiminor, theta=theta


if(n_elements(sample) ne 1) then sample=100.


xt = (2.*dindgen(sample+1)/double(sample) - 1.)*rad

if(n_elements(semiminor) eq 0) then begin
   sub = rad^2 - xt^2
endif else begin
   sub = semiminor^2*(1.-xt^2/rad^2)
endelse

idx = where(sub lt 0)
if(idx[0] gt -1) then sub[idx] = 0.0
yt = sqrt(sub)

x = [xt, reverse(xt), xt[0]]
y = [yt, -reverse(yt), yt[0]]


if(n_elements(theta) eq 1) then begin
   nx = x*cos(theta*!dtor) - y*sin(theta*!dtor)
   ny = x*sin(theta*!dtor) + y*cos(theta*!dtor)
   
   x = nx
   y = ny
   
endif


if(n_elements(xc) eq 1 and n_elements(yc) eq 1) then begin
   x = x + xc
   y = y + yc
endif




end



