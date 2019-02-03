pro whereis,image,w,x,y,z,print=print
;+
; NAME: whereis
;
; PURPOSE:
;   Given the 1-d index of a pixel in an array, return the
;   x and y coordinates corresponding to that pixel.
;
;
; NOTES:
;  pro whereis,image,w,x,y
;
; given the index of a pixel in an array return the
; x and y value
;
; jrg/ucb 94/8/16
; 
; if only one pixel, then return scalars rather than vectors
; 4/16/96 MCL
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 1994-2000 by James R. Graham and Michael C. Liu
;   
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;   
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;   
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this
;    software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;   
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;   
; 3. This notice may not be removed or altered from any source
; distribution.
;   
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;   
;###########################################################################




if n_params() lt 4 and ~(keyword_set(print)) then begin
   message,'pro whereis, image, w, x, y, [z]'
endif

sz = size(image)


y = floor(w/sz(1))
x =  w - y*sz(1)

if n_elements(w) eq 1 then begin
    x =  x(0)
    y =  y(0)
 endif

if keyword_set(print) then print," (X,Y) = "+printcoo(x,y)

end
