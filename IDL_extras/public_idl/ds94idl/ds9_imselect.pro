;+
; NAME: 
;   ds9
;
; PURPOSE:
;   Display images sequentially, allowing user to select good images.
;
; DESCRIPTION:
;   Displays an image cube or list of filenames in ds9 sequentially, and prompts user to select whether image is good or not.  
;   Choices presented to user are y (yes, the default), n (no), and b (back, which moves backwards in the 
;   list without making choices).  Images can be spatially filtered.
;
; INPUTS:
;   ims :  the 3D image cube to display in [dim1, dim2, no_ims] format, or an array of file names.
; 
; INPUT KEYWORDS:
;   usm     :  unsharp mask with a kernel of fwhm=usm
;   sm      :  smooth with a kernel of fwhm=sm
;   iminfo  :  a vector of strings to display for each image, i.e. information about the image. 
;
; OUTPUTS:
;   ims  :  on return, the input array contains only the selected images or file names.
;
; OUTPUT KEYWORDS:
;   index : the original indices of the selected images.
;
; MODIFICATION HISTORY:
;  Written 2013/11/10 by Jared Males (jrmales@email.arizona.edu)
;
; BUGS/WISH LIST:
;
;-
pro ds9_imselect, ims, index=index, usm=usm, sm=sm, iminfo=iminfo


at = size(ims, /type)


fnames = 0
if(at eq 7) then fnames = 1

if fnames then begin
   nims = n_elements(ims)
endif else begin
   get_cubedims, ims, dim1, dim2, nims
endelse


sel = intarr(nims) + 1
   
for i=0,nims-1 do begin

   if fnames then begin
      im = mrdfits(ims[i], /sil)
   endif else begin
      im = ims[*,*,i]
   endelse
   
   if(n_elements(usm) eq 1) then im = im - filter_image(im, fwhm_g=usm)
   if(n_elements(sm) eq 1) then im = filter_image(im, fwhm_g=sm)
   
   ds9, im
        
   prog = 'IMAGE: ' + strcompress(string(i+1) + '/' + string(nims), /rem)
   print, prog
      
   if(n_elements(iminfo) gt 1) then begin
      print, iminfo[i]
   endif

   
   ok_str = 'ok'

   read, ok_str, prompt='Good image? [y/n/b] (y) '
   if(ok_str eq 'n' or ok_str eq 'N' ) then sel[i] = 0
   if(ok_str eq 'b' or ok_str eq 'B' ) then i = i - 2
      
endfor

index = where(sel eq 1)

if fnames then begin
   ims = ims[index]
endif else begin
   ims = ims[*,*,index]
endelse

end


