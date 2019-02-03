;+
; NAME: clio_pca_bgsub
;
; PURPOSE:
;  Performs a PCA b/g subtraction on Clio images
;
; DESCRIPTION:
;   See PURPOSE
;
; INPUTS:
;  fnames_sci   : list of file names of the science images, to be b/g subtracted
;  fnames_sky   : list of file names of the b/g images, to use as the basis set
;  outpath      : the prefix for the output images, which are written to disk
;  nmodes       : vector specifying the number of modes to use 
;
; INPUT KEYWORDS:
;   none
;
; OUTPUTS:
;   klims       :   The K-L images generated from the B/G images
;
; OUTPUT KEYWORDS:
;   none
;
; EXAMPLE:
;  
;
; HISTORY:
;  Written 2016-06-01 by Jared R Males, jaredmales@gmail.com
;          2017-07-10 Doc updated by JRM
;
; BUGS/WISH LIST:
;  Needs much testing
;
;-
pro clio_pca_bgsub, fnames_sci, fnames_sky, outpath, nmodes, klims

clio_getims, fnames_sky, ims=sky, /usef

get_cubedims, sky, dim1, dim2, n_sky

;for i=0,n_sky-1 do begin
;   for j=0,dim1-1 do sky[j,*,i] = sky[j,*,i] - median(sky[j,*,i])
;endfor

sky = reform(sky, dim1*dim2, n_elements(fnames_sky))

for i=0,n_elements(fnames_sky)-1 do sky = sky-median(sky)

pca_covarmat, err, sky ;, /means 


dsubstr = 'clio_pca_bgsub.pro: ' + strcompress(string(nmodes),/rem) + '/' + strcompress(string(n_elements(fnames_sky)),/rem) + ' modes'


pca_klims, klims, err, sky, n_elements(fnames_sky), /silent


bg = fltarr(dim1*dim2)
   
   
for i=0, n_elements(fnames_sci)-1 do begin
   status = 'clio_pca_bgsub: ' + strcompress(string(i+1) + '/' + string(n_elements(fnames_sci)), /rem)
   statusline, status, 0
   

   rim = float(mrdfits(fnames_sci[i], 0, head, /silent, /unsigned))
   rim = rim-median(rim)
   ;for j=0,dim1-1 do rim[j,*] = rim[j,*] - median(rim[j,*])

   rim = reform(rim, dim1*dim2)
 
   cfs = fltarr(nmodes)
   for j=0, nmodes-1 do cfs[j] = klims[*,j]##transpose(rim[*]);

   bg =  cfs[0]*klims[*,0]
   for j=1, nmodes-1 do bg =bg + cfs[j]*klims[*,j]
  
   rim = rim - bg
   rim = rim-median(rim)
 
   rim = reform(rim, dim1, dim2, /over)
   
   
   outname = strcompress(outpath + '/bgsub_' + fnames_sci[i], /remove)

   sxaddpar, head,'HISTORY', dsubstr

   mwrfits, rim, outname, head, /silent
   
endfor

statusline, /clear
print, 'clio_pca_bgsub: done'

klims = reform(klims, dim1, dim2, n_elements(fnames_sky), /over)

end
