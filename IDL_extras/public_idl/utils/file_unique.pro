function file_unique, basename, extension, ndigit=ndigit, startat=startat


if n_elements(ndigit) lt 1 then ndigit = 4
formstr = strcompress('(I0'+string(ndigit)+')', /rem)


if n_elements(startat) eq 1 then begin
   i = startat
   fname = strcompress(basename+string(i, format=formstr) +  extension, /rem)
endif else begin
   fname = strcompress(basename+extension, /rem)
   i = -1
endelse

if file_test(fname) then begin

   i = i + 1

   fname = strcompress(basename+string(i, format=formstr) +  extension, /rem)
   while file_test(fname)  do begin
      i = i+1
      fname = strcompress(basename+string(i, format=formstr) +  extension, /rem)
   endwhile

   
endif 
   
return, fname

end



