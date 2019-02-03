function imagefilt_catch, im, sm_fwhm=sm_fwhm, usm_fwhm=usm_fwhm

get_cubedims, im, dim1, dim2, nims


fim = im

for i=0,nims-1 do begin

   if(n_elements(usm_fwhm) eq 1) then begin
      if(usm_fwhm gt 0) then  fim[*,*,i] = fim[*,*,i] - filter_image(fim[*,*,i], fwhm = usm_fwhm, /all)
   endif


   if(n_elements(sm_fwhm) eq 1) then begin
      if(sm_fwhm gt 0) then fim[*,*,i] = filter_image(fim[*,*,i], fwhm=sm_fwhm, /all)
   endif
endfor

return, fim 

end 



