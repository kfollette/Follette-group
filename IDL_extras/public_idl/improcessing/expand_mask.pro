;+
; NAME: expand_mask
;
; PURPOSE: 
;  expands a binary image mask by replacing each 0 pixel with a circle of 0s with a specified radius
;
; Description:
;  For each 0 pixel in a 1/0 binary mask, this inserts a circle of radius where each pixel is 0.  This
;  then results in an expansion of each masked region.  Will work on a cube of mask images.
;
; INPUTS:
;  maskin : the input mask to expand, possible a cube of format maskin[dim1, dim2, no_ims]
;  rad    : the radius of the replacement circle, the amount by which the mask is expanded
;
; INPUT KEYWORDS:
;  none
;
; OUTPUTS:
;  returns the expanded mask
;
; OUTPUT KEYWORDS:
;  none
;
; DEPENDENCIES:
;   This procedure depends on the following code from JRM:
;      get_cubedims.pro
;      make_apmask.pro
;
; TODO:
;   
;
; HISTORY:
;  2015-04-07: Written by Jared Males, jrmales@email.arizona.edu
;-

function expand_mask, maskin, rad  

mask = maskin

get_cubedims, mask, dim1, dim2, nims

xc = 0.5*(dim1-1)
yc = 0.5*(dim2-1)

xarr = dim1
yarr = dim2

for i=0,nims-1 do begin

   im = mask[*,*,i]

   idx = where(im eq 0)


   for j=0,n_elements(idx)-1 do begin

      pos = array_indices(im, idx[j])

      mdx = make_apmask(pos[0]-xc, pos[1]-yc, rad, xarr, yarr)

      im[mdx] = 0.

   endfor


   mask[*,*,i] = im
endfor

return, mask

end








