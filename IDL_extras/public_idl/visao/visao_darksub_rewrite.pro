pro visao_darksub_rewrite, fnames, imtypes, outpath, region=region, dark=dark, meancomb=meancomb, stdthresh=stdthresh


ddx = where(imtypes eq 2)

print, 'Found ', n_elements(ddx), ' darks'

dark = mrdfits(fnames[ddx[0]],/silent)

imsz = size(dark)
dim1 = imsz[1]
dim2 = imsz[2]


if(n_elements(region) ne 4) then begin
   reg = [0, dim1-1, 0, dim2-1]
endif else begin
   reg = region
endelse


dims = dblarr(dim1, dim2, n_elements(ddx))
dims[*,*,0] = dark
for i=1l, n_elements(ddx)-1 do dims[*,*,i] = mrdfits(fnames[ddx[i]],/silent)

   
if(n_elements(stdthresh) eq 1) then begin
   ;Find standard deviation of each dark frame for rejection
   std = dblarr(n_elements(ddx))
   for i=0l, n_elements(ddx)-1 do std[i] =  stdev(dims[*,*,i])
   ddxstd = where(std le stdthresh)
   print, strcompress('Rejected ' + n_elements(ddx)-n_elements(ddxstd) + ' darks.') 
endif else begin
   ddxstd = indgen(n_elements(ddx))
endelse

darklist = strarr(n_elements(ddxstd))
for i=0, n_elements(ddxstd)-1 do begin
   darklist[i] = (fnames[ddx])[ddxstd[i]]
endfor

if(keyword_set(meancomb)) then begin
   print, 'Forming mean dark . . .'
   for i=0l, n_elements(ddxstd)-1 do begin
      dark = dark + dims[*,*,ddxstd[i]]
      
   endfor
   dark = double(dark)/double(n_elements(ddxstd))
endif else begin
   print, 'Forming median dark . . .'
   dark = median(dims[*,*,ddxstd], dim=3)
endelse

sdx = where(imtypes eq 0)
nims = n_elements(sdx)

print, 'Subtracting dark and re-writing . . .'

if(keyword_set(meancomb)) then begin
   dsubstr = 'dark subtracted with mean of:'
endif else begin
   dsubstr = 'dark subtracted with median of:'
endelse

for i=0l, nims-1 do begin

   status = strcompress(string(i+1) + '/' + string(nims), /rem)
   statusline, status, 0
   
   im = double(mrdfits(fnames[sdx[i]], 0, header, /silent)) - dark

   outname = strcompress(outpath + '/dsub_' + fnames[sdx[i]], /remove)

   sxaddpar, header,'HISTORY', dsubstr
   for j=0, n_elements(darklist)-1 do sxaddpar, header, 'HISTORY', darklist[j]
   
   mwrfits, im[reg[0]:reg[1], reg[2]:reg[3]], outname, header, /silent
   
endfor

statusline, /clear
print, 'Done'

end





