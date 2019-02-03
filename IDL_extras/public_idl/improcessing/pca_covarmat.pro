;+
; NAME: 
;   pca_covarmat
;
; DESCRIPTION:
;   Calculates the covariance matrix for a cube of images
;
; INPUTS:
;   rims     :  a cube of images with format [dim1, dim2, number_of_images] (possibly with a masked region removed). 
;
; INPUT KEYWORDS:
;   none
;
; KEYWORDS:
;   meansub  : if set, the mean of each image is first subtracted, modifying ims
;              this is the strict interpretation of Soummer et al.
;
; OUTPUTS:
;   err   : the covariance matrix with size [number_of_images, number_of_images] 
;
;
; MODIFICATION HISTORY:
;  Written 2013/01/26 by Jared Males (jrmales@email.arizona.edu)
;
; BUGS/WISH LIST:
;
;-
pro pca_covarmat, err, rims, meansub=meansub, mask=mask


;Do mean subtraction if requested
if(keyword_set(meansub)) then begin

   sz = size(rims)

   nims = sz[2]

   ;Is there a mask?
   domask = 0
   if(n_elements(mask) gt 1) then domask = 1
   
   for i=0, nims-1 do begin
      
      if(~domask) then begin
         mn = mean(rims[*,i])
         rims[*,i] = rims[*,i] - mn
      endif else begin
         idx = where(mask[*,i] ne 0, comp=zdx)
                  
         ;mn = mean(rims[idx,i])
         
         rims[idx,i] = rims[idx,i] - mean(rims[idx,i]); mn

      endelse
      
      
   endfor
   
endif

;Now calculate the covariance matrix
err = rims ## transpose(rims)



end

