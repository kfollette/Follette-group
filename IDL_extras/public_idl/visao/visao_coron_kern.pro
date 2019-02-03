function visao_coron_kern, dim1, dim2, azsig, pixfwrat, xarr, yarr, x, y, theta, idx, flip=flip

idx = where(abs(xarr-x) le 10*azsig and abs(yarr-y) le 10*azsig)

c  = cos(theta)
s  = sin(theta)
s2 = sin(2*theta)


wx = azsig*pixfwrat
wy = azsig

if(keyword_set(flip)) then begin
   wy = wx
   wx = azsig
endif

kern = dblarr(dim1, dim2)

if(idx[0] eq -1) then return, kern

for k=0, n_elements(idx)-1 do begin
   xp = xarr[idx[k]] - x
   yp = yarr[idx[k]] - y
         
         
   a = .5*(c^2/wx^2 + s^2/wy^2)
   b = .25*(-s2/wx^2 + s2/wy^2)
   d = .5*(s^2/wx^2 + c^2/wy^2)
   
   U = a*xp^2 + 2.*b*xp*yp + d*yp^2
      
   kern[idx[k]] = exp(-U)

endfor

kern = kern/total(kern)

return, kern

end
