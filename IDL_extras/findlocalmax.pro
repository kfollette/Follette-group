;+
; NAME:
;   findlocalmax
;
; PURPOSE:
;   Find local maximum near given position in image
;
; CALLING SEQUENCE:
;   maxind = findlocalmax(im, ind, boxsize=)
;
; INPUTS:
;   im      - image 
;   ind     - 1-D index of pixel in image
;
; KEYWORD PARAMETERS:
;   boxsize - size of box in which to explore, in pixels.  (default
;             15)
;
; OUTPUTS:
;   returns maxind, 1-D index of local max near input ind
;
; COMMENTS:
;   Uses 1-D index of a 2-D array, i.e. 0 < ind < N*M-1 for NxM image.
;   This algorithm searches a box of size boxsize x boxsize centered
;   on pixel ind, finding the max within that box.  This iterates
;   until the maximum pixel no longer moves.  Boxsize should be bigger
;   than the typical cosmic ray size, but much smaller than the
;   typical distance between cosmic rays (and other objects). 
;
; MODIFICATION HISTORY:
;   2000 Oct 4 - written by D. Finkbeiner
;   2001 07 25 - keyword silent added. MDP
;   2001 07 26 - modified to accept x,y instead of 1-d index by
;                default. MDP
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2000-2003 by Doug Finkbeiner & Marshall Perrin
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
 
function findlocalmax, im, x,y,ind=ind, boxsize=boxsize,silent=silent


  sx = (size(im))[1]
  sy = (size(im))[2]

if keyword_set(ind) then begin
  x = ind mod sx
  y = ind  /  sx
endif

if x lt 0 or y lt 0 then message,'ERROR: x and y should not be negative!'

if NOT keyword_set(boxsize) then boxsize = 15
  boxrad = (boxsize-1)/2


  repeat begin 
     x0 = (x-boxrad)>0 ; deal with image edges appropriately
     y0 = (y-boxrad)>0
     patch = im[x0:(x+boxrad)<(sx-1),y0:(y+boxrad)<(sy-1)]
     pxsize = (size(patch))[1]
     w = where(patch EQ max(patch))
     dx = (w[0] mod pxsize) - (x-x0)
     dy = (w[0]  /  pxsize) - (y-y0)
     x = ((x+dx) > 0) < (sx-1)
     y = ((y+dy) > 0) < (sy-1)
     if not(keyword_set(silent)) then print,  x, y, dx, dy
  endrep until (dx EQ 0) AND (dy EQ 0)

  maxind = x+y*sx

  return, maxind
end
