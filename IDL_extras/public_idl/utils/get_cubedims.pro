;+
; NAME: get_cubedims
;
; PURPOSE: 
;  Get dimensions of an image cube.
;
; Description:
;  Calculates the dimensions of an image cube, of form [dim1, dim2, nims], where nims is the 
;  number of images.  Really just a wrapper for size().
;  Also works for 2D images, then nims will be 1
;
; INPUTS:
;  ims  :  the image cube
;
; OUTPUTS:
;  dim1 : the size of dimension 1, (size(ims))[1]
;  dim2 : the size of dimension 2, (size(ims))[2]
;  nims : the number of images (size(ims))[3]
;
; HISTORY:
;  2012-03-16: Written by Jared Males, jrmales@email.arizona.edu
;  2012-04-10: Jared Males added handling for 2d 
;-
pro get_cubedims, ims, dim1, dim2, nims

B = size(ims)

dim1 = B[1]
dim2 = B[2]

if(B[0] eq 2) then begin
   nims = 1
endif else begin
   nims = B[3]
endelse

end
