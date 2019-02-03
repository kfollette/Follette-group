;+
; NAME: subradprofim
;
; PURPOSE: 
;  Subtract the median radial profile from each image in im
;
; Description:
;  For each image in im, the median radial profile is calculated, and then subtracted.
;
; INPUTS:
;  im     : the image(s) to subtract, possibly a cube of format im[dim1, dim2, no_ims]
;
; INPUT KEYWORDS:
;  r      : a radius array, as produced by rarr, of the same size as the input image (or one image in the cube, 
;           can be empty in which case is is initialized to a centered array
;
; OUTPUTS:
;  On output im will have its radial profile subtracted
;
; OUTPUT KEYWORDS:
;  none
;
; DEPENDENCIES:
;   This procedure depends on the following code from JRM:
;      get_cubedims.pro
;      radprofim.pro
;
; TODO:
;   
;
; HISTORY:
;  2015-02: Written by Jared Males, jrmales@email.arizona.edu
;-
pro subradprofim, im, r=r

get_cubedims, im, dim1, dim2, nims

for i=0,nims-1 do im[*,*,i] = im[*,*,i] - radprofim(im[*,*,i], r=r)

end


