;+
; NAME: rot_cube
;
; PURPOSE: 
;  Rotate a cube of images.
;
; Description:
;  A cube of images is rotated by the specified angles, possibly with masking.  Rotation is handled by
;  the rot function.  NOTE: angles here are specified counter-clockwise-positive.
;
; INPUTS:
;  ims     : the images to rotate
;  rotangs : the angles by which to rotate the images, in counter-clockwise degrees.
;
; INPUT KEYWORDS:
;  mask    : an optional 1/0 binary mask which is also rotated.  can be a cube.
;  m0val   : the value to apply to the post-rotation masked pixels, default is nan
;
; OUTPUTS:
;  On output ims contains the rotated images.
;
; OUTPUT KEYWORDS:
;  silent  : shut up.
;
; DEPENDENCIES:
;   This procedure depends on the following code from JRM:
;      get_cubedims.pro
;      
;
; TODO:
;   
;
; HISTORY:
;  2013-ish: Written by Jared Males, jrmales@email.arizona.edu
;  2015-04-07: Docs updated.
;-
pro rot_cube, ims, rotangs, mask=mask, m0val=m0val, silent=silent


nims = (size(ims))[3]

if(n_elements(m0val) eq 0) then m0val = float('nan')

besilent = 0
if(keyword_set(silent)) then besilent = 1


if(n_elements(mask) gt 1) then begin
   get_cubedims, mask, maskd1, maskd2, maskn
endif


for i=0, nims-1 do begin

   if(~besilent) then begin
      status = 'rot_cube: ' + strcompress(string(i+1) + '/' + string(nims), /rem)
      statusline, status, 0
   endif
   
   ims[*,*,i] = rot(ims[*,*,i], -1.d*double(rotangs[i]), cubic=-0.5)

   if(n_elements(mask) gt 1) then begin
   
      if(maskn eq 1) then begin
         rmask = rot(mask, -1.d*double(rotangs[i]), cubic=-0.5)
      endif else begin
         rmask = rot(mask[*,*,i], -1.d*double(rotangs[i]), cubic=-0.5)
      endelse
      
      idx = where(rmask lt 0.2)
      if(idx[0] gt -1) then rmask[idx] = m0val
      
      ims[*,*,i] = ims[*,*,i]*rmask
   endif
      
endfor

if(~besilent) then begin
   statusline, /clear
endif

end

