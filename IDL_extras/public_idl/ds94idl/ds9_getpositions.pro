pro ds9_getpositions, ims, x, y, smoothim=smoothim, usm=usm, nodoffset=nodoffset

common ds9env, $    ;this is the ds9 environment common block
       ds9obj    ;the ds9 object
       
get_cubedims, ims, dim1, dim2, nims


x = dblarr(nims)
y = dblarr(nims)

for i=0, nims-1 do begin

   im = ims[*,*,i]

   if(keyword_set(nodoffset)) then begin
      ndf = nodoffset
      
      if (i+ndf gt nims-1) then ndf = -1.*ndf
      print, ndf
      
      im = im - ims[*,*, i+ndf]
   endif
   
   if(keyword_set(usm)) then begin
      im = im - smooth(im, usm)
   endif 
   
   if(keyword_set(smoothim)) then begin
      im = smooth(ims[*,*,i], smoothim)
   endif 
   ;else begin
   ;   im = ims[*,*,i]
   ;endelse
   
   ds9, im

   ds9obj->basic_imexam, tkey, xc, yc, zc
   key = tkey
   
   if(key eq 'c') then begin
      x[i] = xc
      y[i] = yc

      print, i, xc, yc
   endif else begin
      if (key eq 'n') then begin
         x[i] = -1;
         y[i] = -1
         print, i, -1, -1
      endif
      
      if (key eq 'q') then return
      if (key eq 'b') then i = i - 2
      if(i lt -1) then i = -1
   endelse
   
endfor


end



