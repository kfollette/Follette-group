pro whereismax, image, x, y, maxv, mark=mark, silent=silent,single=single
;+
; NAME:  whereismax
;
; PURPOSE: 
; Given an array, returns location and value of max pixel
;
; INPUTS:
;
; pro Whereismax, image, x, y, maxv, mark=mark,
; silent=silent,single=single
;
;
; note that this will return a long integer if there is
; a unique max pixel, otherwise it will return an long array
; which may be a problem for some routines
;
; 10/31/94 MCL
;
; added mark feature, which will draw a cross at the max
; pixel on the current window
; 8/19/96 MCL
;
; uses a square instead of a circle for /mark
; 06/22/00 MCL
; Added "single" keyword to always return only one value. If multiple
;  points have the same max, return only the first
; 2001-07-06 MDP
;
; 2003-03-04 MDP added /NaN keyword to max
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 1994-2003 by Michael C. Liu & Marshall D. Perrin
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
 

if n_params() lt 1 then begin
   message,'pro whereismax,image,(x),(y),(maxv),[mark],[silent],[single]'
endif


maxv = max(image,/NaN)

whereis,image,where(image eq maxv),x,y

if ((n_elements(x) eq 1) or  keyword_set(single) ) then begin
   if (n_elements(x) gt 1 and keyword_set(single) and ~(keyword_set(silent))) then $
      message,"multiple maxes found; returning only the first",/info
   x = x(0)
   y = y(0)
endif

if keyword_set(silent) eq 0 then $
  print,'max = ',strc(maxv),' at (',strc(fix(x)),',',strc(fix(y)),')'

sz = size(image)
rad = 0.05*avg(sz(1:2))

if keyword_set(mark) then begin
    if (!d.window ne -1) then begin
      oplot, [x], [y], sym = 10, ps = 1, color = 1
      tvbox, 2*rad, x, y, /data, color = 1
      ;tvcircle, 40, x, y, /data, color = 1
   endif else $
      message, 'unable to mark b/c no windows available'
 endif

end
