;+
;  FUNCTION: findmaxstar
;  PURPOSE:
;  given an image, finds the x,y coords of the brightest object
;  present.
;  DESCRIPTION:
;  Uses simple binning as a rough CR rejection routine.
;  Can do more detailed rejection of CRs using qzap
;
;  This is intended to be used for the automatic registration of
;  images for
;  mosaicing.
;  Empirical evidence is that this routine is negliglbly slower than
;  just using whereismax, but is much more robust against hot pixels
;  or CRs which are still in the data. Has not been tested on crowded 
;  fields, where source confusion would probably make it break down.
;
;  It also automatically interpolates values for any bad (NaN) pixels
;  in the image; NaNs will be set to the image median so that the code
; won't choke.
;
;  INPUTS:
;  imgan image
;  /silentwhether to print any messages or not
;  OUTPUTS:
;
;  REQUIREMENTS:
;  needs whereismax.pro and findlocalmax.pro
;  
;  HISTORY:
;  Started July 2001 by Marshall Perrin
;  2002-12-05Added check for NaN pixels.
;
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 2001-2003 by Marshall Perrin
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
    



pro findmaxstar, img, x, y,silent=silent,sigma=sigma,qzap=qzap,reducf=reducf

  if (keyword_set(qzap)) then begin
     message,/info, "Removing CRs..."
     qzap,img,i 
  endif else i=img

  ; definitely need to fix NaNs before rebinning
  wnf = where(finite(i,/NaN),wnfcount)
  if wnfcount gt 0 then begin
     if not(keyword_set(silent)) then message,/info,"Fixing NaNs"
     ;fixpix,i,0,fi,/NaN
     ;i = temporary(fi)
     i[wnf] = median(i)
  endif
  
                                ; reduce the image to something around
                                ; 16x16. Based on correl_images.pro 
  ; from the Goddard IDL Astronomy Library.
  sz=size(i)
  if not(keyword_set(reducf)) then reducf=min(sz[1:2])/16
  
                                ; but don't let reducf be
                                ; greater than 16. That hurts
                                ; our chances too
  ; much of finding the correct peak when we resize
  reducf = reducf < 16
    sz = sz/reducf
    LA = sz * reducf -1   ;may have to drop edges of images.

    imtmp = Rebin( i[ 0:LA[1], 0:LA[2] ], sz[1], sz[2] )

    ; find brightest star in binned down image (doesn't find CRs)
    whereismax,imtmp,x,y,/silent,/single
    x=x*reducf
    y=y*reducf
    ; now find the actual brightest pixel in the full image
    maxpix=findlocalmax(i,x,y,/silent)
    whereis,i,maxpix,x,y

    if keyword_set(silent) eq 0 then $
         print,'max found at (',strc(fix(x)),',',strc(fix(y)),')'



end
