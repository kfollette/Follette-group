function pad_ims, ims, width, padval=padval

get_cubedims, ims, dim1, dim2, nims


if(n_elements(padval) ne 1) then padval = 0
padims = fltarr(dim1+2*width, dim2+2*width, nims) + padval

for i=0, nims-1 do padims[width:width+dim1-1, width:width+dim2-1, i] = ims[*,*,i]

return, padims

end



