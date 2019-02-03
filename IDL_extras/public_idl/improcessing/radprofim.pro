;+
; NAME: radprofim
;
; PURPOSE: 
;  Calculate the median radial profile of an image
;
; Description:
;  A median radial profile is calculated by binning and linear interpolation
;
; INPUTS:
;  im     : the image to make a profile from
;
; INPUT KEYWORDS:
;  r      : a radius array, as produced by rarr, of the same size as the input image (or one image in the cube, 
;           can be empty in which case is is initialized to a centered array
;
; OUTPUTS:
;  Returns the radial profile.
;
; OUTPUT KEYWORDS:
;  none
;
; DEPENDENCIES:
;   This procedure depends on the following code from JRM:
;      get_cubedims.pro
;      bin_median.pro
;
; TODO:
;   
;
; HISTORY:
;  2015-02: Written by Jared Males, jrmales@email.arizona.edu
;-
function radprofim, im, r=r

get_cubedims, im, dim1, dim2

if(n_elements(r) ne dim1*dim2) then r=rarr(dim1,dim2, /pix)
mr = max(r)

bin_median, r, im, 0., mr+1, mr+1, binr, binmed 

linterp, binr, binmed, reform(r, dim1*dim2), rprof

rprof=reform(rprof, dim1, dim2, /over)

return, rprof

end


