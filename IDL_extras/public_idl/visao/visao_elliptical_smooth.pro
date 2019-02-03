function visao_elliptical_smooth, im, azfw, fwrat,  minr=minr, maxr=maxr


b = size(im)

dim1 = b[1]
dim2 = b[2]

r = rarr(dim1, dim2,x, y, theta=theta, /pixels)

if(n_elements(minr) ne 1) then minr = 0
if(n_elements(maxr) ne 1) then maxr = max(r)


outim = dblarr(dim1, dim2)

azsig = (azfw/(2.*sqrt(2.*alog(2))))

for i=0, dim1-1 do begin

   for j=0, dim2-1 do begin
      
      if(r[i,j] lt minr or r[i,j] gt maxr) then continue
      
      pixfwrat = fwrat
      
      kern = visao_coron_kern(dim1, dim2, azsig, pixfwrat, x, y, x[i,j], y[i,j], .5*!pi - theta[i,j], idx, /flip)
      outim[i,j] = total(im[idx]*kern[idx])

   endfor
endfor

return, outim

end
