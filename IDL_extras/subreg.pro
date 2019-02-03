;+
;
; FUNCTION: subreg
; PURPOSE:
;subpixel registration of images
;
; NOTES:
;
;Registers multiple images against a reference image based on one of
;several methods for determining the shift between the two images at a
;subpixel
;level. Returns the offsets for each image relative to the first.
;
; PRO subreg,ref,imgs0,shifts,badpix=badpix,silent=silent,$
; goodlist=goodlist0,method=method,boxsize=boxsize
;
; INPUTS:
; ref    Reference image to which all other images are registered.
; imgs0   array of images. Should be bias corrected and flat-fielded.
;               must have same x and y dimensions as ref.
; KEYWORDS:
; badpixbad pixel mask
; goodlistList of which images in the cube are good (=1)
; Only these are registered against the first good one.
; headersan array of FITS headers, used for the "H" alignment
; methodWhich method to use for the subpixel registration?
; "X"Fourier cross-correlation, centroided
; "F"Fourier cross-correlation, fit by gaussian
; "FHP"Fourier cross-correlation, with high pass filter, fit by
;gaussian
; "C"center-of-mass centroid
; "O"Goddard astro library CNTRD routine
; "G"Gaussian fitting iteratively
; "R"Mike Liu's recenter.pro centroid
; "N"no shifts; degenerate case: just return 0
; "M"maximum pixel intensity
; "H"read shifts from FITS headers. Currently shifts
; are relative to the absolute zero point of the
; FITS headers; ref0 is ignored in this case.
; "HF"    Hybrid: use Fourier, but let header values
; override if they're wildly different.
; Default is "F"
;
;Based on empirical evidence, usually "X" or F works best, 
;but which one is best for any specific case
;will depend on the details of your data in that specific instance.
; 
; boxsize box size for fit to cross correlation
; whichref=Which header to use for the reference image? (only used
;forHF)
; OUTPUTS:
; shifts2d array of [x,y] shifts of each image relative to the first
;
; HISTORY:
;Marshall Perrin, starting Aug 2001. Based on pixel-level registration
;code by Mike Liu
;2001-08-22Merged code from pxc.pro
;2001-08-27Split ref into a seperate input, added error checking.
;2002-12-05      Added code to handle if NaNs are present in the
;images
;2003-11-25Added "H" option.
;2004-04-07Renamed Goddard cntrd to "O". Added Gaussian option.
;2005-10-16    Changed all gauss2dfits to mpfit2dpeak.
;and lots of other stuff I haven't kept track of.
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
;software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source
;distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################



PRO subreg,ref0,imgs0,shifts,badpix=badpix,silent=silent,$
           goodlist=goodlist0,method=method,boxsize=boxsize,headers=headers,$
           display=display,whichref=whichref

if n_params() ne 3  then begin
   print,"Usage:  subreg,ref,imgs0, shifts,[badpix=],[silent=]"
   print,"[goodlist=],[centroid=],[boxsize=] "
   stop
endif

sr=size(ref0)
if sr[0] lt 2 then $
   message,"First argument (reference image) must be an array!"
sz=size(imgs0)
if sz[0] gt 2 then ntot=sz[3] else ntot=1
if (sz[1] ne sr[1]) or (sz[2] ne sr[2]) then $
   message,"Reference and images must have same X and Y dimensions!"
; well, they don't really *have* to, but the outputs don't make much
; sense if they don't. There may be some circumstances where you'd
; want
; to do such a thing, but I can't think of any right now...
nx=sz[1]
ny=sz[2]

if not(keyword_set(method)) then method="F" 
if strc(method) eq strc(1) then method="F" ; this catches the case of calling with /method

if not(keyword_set(boxsize)) then begin 
   if (method eq "F") or (method eq "G") then boxsize = 40 else boxsize=20
endif

if keyword_set(goodlist0) then begin
    if n_elements(goodlist0) ne sz(3) then  $
      message, 'goodlist does not have right number of entries!^G'
    goodlist = goodlist0
 endif else  $
  goodlist = bytarr(ntot)+1B
wg = where(goodlist ne 0, ng)

shifts=fltarr(2,ntot)
gshifts=fltarr(2,ntot)
if strc(method) eq "N" then return ; null case - for testing purposes.

;if sz[1] ne sr[1] or sz[2] ne sr[1] then $
;message,"reference image must have same dimensions as images to
;register!"
if sr[0] gt 2 then $
   message,"Reference image must be only 2D, not larger."

if strc(method) eq "H" or strc(method) eq "HF" then begin
   if n_elements(headers) eq 0 then message,"For registration method H, you must supply an array of FITS headers!"
   head0 = headers[*,0]
   if (strc(sxpar(head0,"INSTRUME")) ne "IRCAL") then $
      message,"Method H currently only supported for IRCAL images and headers."
   fitsoffsetx="NOD_EW"
   fitsoffsety="NOD_NS" ; these are in arcsec
   fitsoffsetscale=13.4 ; pixels per arcsec
   fitsoffsetscale=1./0.040 ; pixels per arcsec
   
   ref_offset = [$
                sxpar(headers[*,whichref],fitsoffsetx)*fitsoffsetscale*(-1),$
                sxpar(headers[*,whichref],fitsoffsety)*fitsoffsetscale $
                ] 
endif 

if strc(method) ne "H" then  begin
   ; this is needed for all methods other than H
   ; including HF
   ref = ref0
   wnan = where(finite(ref,/NAN),nanct)
   while(nanct gt 0) do begin ; sometimes this takes more than one pass.
      if not(keyword_set(silent)) then message,/info,"fixing bad pixels before cross-correlation"
      fixpix,ref,0,ref2,/NaN,/quick
      ref=temporary(ref2)
      wnan = where(finite(ref,/NAN),nanct)
   endwhile
   
   findmaxstar,ref,x1,y1,/silent
endif


; Do we need to register the first image?
; If the first image is identical to the reference, then don't bother.
;if max(abs((imgs0[*,*,0] - ref0)/ref0),/nan) lt 0.01 then startreg =
;                           1 else startreg=0
;startreg=0
if array_equal(imgs0[*,*,0],ref0) then startreg=1 else startreg=0



for i=startreg,ng-1 do begin
   n2=wg[i]
   imagei = imgs0[*,*,n2]
   wnan = where(finite(imagei,/NAN),nanct)
   while(nanct gt 0) do begin ; sometimes this takes more than one pass.
      fixpix,imagei,0,imagei2,/NaN,/silent,/quick
      imagei=temporary(imagei2)
      wnan = where(finite(imagei,/NAN),nanct)
   endwhile
   
   if not(keyword_set(silent)) then print,"Registering image "+strc(i)+"/"+strc(ng)+" using method "+method, format = '($,A," ")'
   if keyword_set(display) then display,alogscale(imagei),title="Registering image "+strc(i),/silent

;-----------------------------------------------------------------
   if strc(method) eq "X" then begin ; fourier cross-correlation, fit by centroid.
      
      
      
      cor = convolve(imagei,ref,/correlate,FT_PSF=ftpsf)
      findmaxstar,cor,xi,yi,/silent; get rough center
      mrecenter,cor,xi,yi,x,y,/silent,/nodisp  ; get fine center
      ;print,xi,yi,x,y,nx/2-x,ny/2-y
                                ;s=boxsize
                                ;g=gauss2dfit(cor[0>x-s:x+s<(sz[1]-1),0>y-s:y+s<(sz[2]-1)],a)

      
      shifts[0,n2]=round(nx/2.)-x
      shifts[1,n2]=round(ny/2.)-y
                                ;atv,[[[ref]],[[imagei]],[[shift(imagei,shifts[0,n2],shifts[1,n2])]]],/bl

      ;gshifts[0,n2]=nx/2-(0>x-s)-a[4]
          ;gshifts[1,n2]=ny/2-(0>y-s)-a[5]

      ;print,printcoo(shifts[*,n2]),printcoo(gshifts[*,n2])
      ;stop
      findmaxstar,imagei,x2,y2,/silent
      ;print,shifts[*,n2]
      ;TODO code here to make sure stuff doesn't 'wrap'
      if x2+shifts[0,n2] lt 0 then  shifts[0,n2]+=sz[1]
      if x2+shifts[0,n2] ge sz[1] then shifts[0,n2]-=sz[1]
      if y2+shifts[1,n2] lt 0 then  shifts[1,n2]+=sz[2]
      if y2+shifts[1,n2] ge sz[2] then shifts[1,n2]-=sz[2]
      ;print,shifts[*,n2]
   endif
;-----------------------------------------------------------------
   if strc(method) eq "F" then begin ; old-style fourier cross-correlation, fit by Gaussian
      cor = convolve(imagei,ref,/correlate,FT_PSF=ftpsf)
          findmaxstar,cor,x,y,/silent
          s=boxsize
                                ;g=gauss2dfit(cor[0>(x-s):(x+s)<(sz[1]-1),0>(y-s):(y+s)<(sz[2]-1)],a)
              g=mpfit2dpeak(cor[0>(x-s):(x+s)<(sz[1]-1),0>(y-s):(y+s)<(sz[2]-1)],a)
                                ;if keyword_set(display )then begin
                                ;wset,0 &
                                ;shade_surf,cor[0>x-s:x+s<(sz[1]-1),0>y-s:y+s<(sz[2]-1)],$
                                ;title="Fourier Cross-correlation" &
                                ;endif
                                ;if keyword_set(display )then begin
                                ;wset,1 &  shade_surf,g,title="fit
                                ;gaussian" & endif
              if a[4] lt 0 or a[4] gt ( (x+s<(sz[1]-1)) -(0>x-s) ) then print,"X fit error"
              if a[5] lt 0 or a[5] gt ( (y+s<(sz[1]-1)) -(0>y-s) ) then print,"Y fit error"
              ; TODO add code here to refit when there are errors.
              
                                ; the shift is going to be the
                                ; difference between the image center,
                                ; nx/2
                                ; and the gaussian peak location (plus
                                ; fitting box offset),  (0>(x-s))+a[4]
              ;
                                ; rounding code added 2006-05-25 to
                                ; deal with images with odd-sized axes
                                ; Oddly, IDL does NOT return a
                                ; convolution which is centered on a
              ; half-pixel, for some reason.
              ; Should work the same as always for even sizes.
              shifts[0,n2]=round(nx/2.)-(0>(x-s))-a[4]
                  shifts[1,n2]=round(ny/2.)-(0>(y-s))-a[5]

                  ;findmaxstar,imagei,x2,y2,/silent

; 2005-12-19 The following code has to be wrong.
; Should it be using shifts[0,n2] etc?
; Or is it just entirely unnecessary?
; 2006-07-15  No, sometimes it IS necessary. I'm confused why!
; But it is for, say, PDS 415 on 2006 Jun 08, K band
                  ;TODO code here to make sure stuff doesn't 'wrap'
                  findmaxstar,imagei,x2,y2,/silent
                  if x2+shifts[0] lt 0 then shifts[0]=shifts[0]+sz[1]
                  if x2+shifts[0] ge sz[1] then shifts[0]=shifts[0]-sz[1]
                  if y2+shifts[1] lt 0 then shifts[1]=shifts[1]+sz[2]
                  if y2+shifts[1] ge sz[2] then shifts[1]=shifts[1]-sz[2]
               endif

;-----------------------------------------------------------------
   if strc(method) eq "FHP" then begin ; old-style fourier cross-correlation, fit by Gaussian
      ; HIGH PASS FILTER
      if ~(keyword_set(ftpsf)) then begin
                                ; we can effect the high pass
                                ; filtering in Fourier space just by
         ; multiplying the FFT of the PSF with the filter function.
         ; they both will multiply by the FFT of the image in the
         ; convolution.
         ftpsf = fft(ref,-1)
         highpassfilter = 1 - 1.0 / ( 1.0d + (dist(sz[1],sz[2])/2.0)^2 )
         ftpsf *= highpassfilter
      endif
      cor = convolve(imagei,ref,/correlate,FT_PSF=ftpsf)
          findmaxstar,cor,x,y,/silent
          s=boxsize
              g=mpfit2dpeak(cor[0>(x-s):(x+s)<(sz[1]-1),0>(y-s):(y+s)<(sz[2]-1)],a)
              if a[4] lt 0 or a[4] gt ( (x+s<(sz[1]-1)) -(0>x-s) ) then print,"X fit error"
              if a[5] lt 0 or a[5] gt ( (y+s<(sz[1]-1)) -(0>y-s) ) then print,"Y fit error"
              ; TODO add code here to refit when there are errors.
              
              shifts[0,n2]=round(nx/2.)-(0>(x-s))-a[4]
                  shifts[1,n2]=round(ny/2.)-(0>(y-s))-a[5]
                  ;stop
               endif
;-----------------------------------------------------------------
   if strc(method) eq "HF" then begin ; HYBRID fourier cross-correlation, fit by Gaussian
      cor = convolve(imagei,ref,/correlate,FT_PSF=ftpsf)
          findmaxstar,cor,x,y,/silent
          s=boxsize
                                ;g=gauss2dfit(cor[0>(x-s):(x+s)<(sz[1]-1),0>(y-s):(y+s)<(sz[2]-1)],a)
              g=mpfit2dpeak(cor[0>(x-s):(x+s)<(sz[1]-1),0>(y-s):(y+s)<(sz[2]-1)],a)
                                ;if keyword_set(display )then begin
                                ;wset,0 &
                                ;shade_surf,cor[0>x-s:x+s<(sz[1]-1),0>y-s:y+s<(sz[2]-1)],$
                                ;title="Fourier Cross-correlation" &
                                ;endif
                                ;if keyword_set(display )then begin
                                ;wset,1 &  shade_surf,g,title="fit
                                ;gaussian" & endif
              if a[4] lt 0 or a[4] gt ( (x+s<(sz[1]-1)) -(0>x-s) ) then print,"X fit error"
              if a[5] lt 0 or a[5] gt ( (y+s<(sz[1]-1)) -(0>y-s) ) then print,"Y fit error"
              ; TODO add code here to refit when there are errors.
              
                                ; the shift is going to be the
                                ; difference between the image center,
                                ; nx/2
                                ; and the gaussian peak location (plus
                                ; fitting box offset),  (0>(x-s))+a[4]
              shifts[0,n2]=nx/2-(0>(x-s))-a[4]
                  shifts[1,n2]=ny/2-(0>(y-s))-a[5]

                  findmaxstar,imagei,x2,y2,/silent

                  ; Now the hybrid part:
                  hshifts = [$
                            sxpar(headers[*,i],fitsoffsetx)*fitsoffsetscale*(-1),$
                            sxpar(headers[*,i],fitsoffsety)*fitsoffsetscale $
                            ] - ref_offset
                  ;if max(abs(hshifts -shifts[*,n2])) gt 20 then begin
                                ; changed to 40 on 2006-06-09;
                                ; problems with offset drift during
                                ; LGS
                  ; mode.
                  if max(abs(hshifts -shifts[*,n2])) gt 40 then begin
                      print,strc(shifts[*,n2])
                      print, "    Overriding with FITS header values: ",format = '($,A," ")'
                      shifts[*,n2] = hshifts
                   endif
                  
               endif

;-----------------------------------------------------------------
   if strc(method) eq "C" then begin ; StarFinder method
      s=boxsize
      findmaxstar,imagei,x2,y2,/silent
      
      c1 = centroid(ref[x1-s:x1+s,y1-s:y1+s])
      c2 = centroid(imagei[x2-s:x2+s,y2-s:y2+s])
      dc=(c1+ [x1-s,y1-s] )-(c2 + [x2-s,y2-s]) 
                                ;if not(keyword_set(silent)) then
                                ;print,["dxc= ","dyc= "]+strc(dc)
      shifts[*,n2]=dc

   endif
;-----------------------------------------------------------------
   if strc(method) eq "M" then begin ; 
      s=boxsize
      findmaxstar,imagei,x2,y2,/silent
      
      shifts[*,n2]=[x1-x2,y1-y2]
   endif

;-----------------------------------------------------------------
   if strc(method) eq "G" then begin ; gaussian fitting method
      findmaxstar,imagei,x2,y2,/silent

      s=boxsize < (x1-1) < (y1-1) < (sz[1]-x1) < (sz[2]-y1)
      s=s       < (x2-1) < (y2-1) < (sz[1]-x2) < (sz[2]-y2)
      
      im1 = ref[x1-s:x1+s,y1-s:y1+s]
      im2 = imagei[x2-s:x2+s,y2-s:y2+s]

      roughshift = [x1-x2,y1-y2]
      ;print,"Rough:",roughshift[0],roughshift[1]

      res1 = mpfit2dpeak(im1,g1)
      res2 = mpfit2dpeak(im2,g2)
      dxy = [g1[4]-g2[4],g1[5]-g2[5]]
      ; elements of g are:
      ; const height xwid ywid xcen ycen

          c_iter=0
          ftol=0.005
              repeat begin
                 ;print,"dxy: "+printcoo(dxy)
                 ;print,"iter:",c_iter,"dxy:",dxy
                 if (total(dxy) gt 20) then begin
                    print,"Image shifted by more than 20 pixels from our initial guess position. Check?"
                    print,"try this:   atv,[[[im1]],[[im2s]]],/bl"
                    stop
                 endif
                 im2s=image_shift(im2,dxy[0],dxy[1],interp_type=interp_type,data)
                                ;subreg,im1,im2s,shifts,/silent,method=method,boxsize=boxsize;
                                ;calculate the new relative shifts
                         res2 = mpfit2dpeak(im2s,g2)
                                 s = [g1[4]-g2[4],g1[5]-g2[5]]
                                         dxy=dxy+s          ; new overall shifts for next round
                                                 error=sqrt(total(s^2))
                                                 c_iter=c_iter+1
                                ;print,"Iter "+strc(c_iter)+",
                                ;shifts="+strc(s[0])+", "+strc(s[1])+"
                                ;Error: "+strc(error)
                                              endrep until (error lt ftol)
              ;atv,[[[im1]],[[im2s]],[[im2]]],/bl
              delvarx,data
              
                  ;message,/info, "Image registration complete"
                  dxy = dxy+roughshift

                  shifts[*,n2]=dxy
               endif
   
;-----------------------------------------------------------------
   if strc(method) eq "O" then begin ; goddard method
      s=boxsize
      findmaxstar,imagei,x2,y2,/silent

      cntrd,ref,x1,y1,xc1,yc1,5
      cntrd,imagei,x2,y2,xc2,yc2,5
                                ;if not(keyword_set(silent))
                                ;thenprint,["dxn= ","dyn=
                                ;"]+strc([xc1,yc1]-[xc2,yc2])
      shifts[*,n2]=[xc1,yc1]-[xc2,yc2]
   endif
;-----------------------------------------------------------------
   if strc(method) eq "R" then begin ; Mike Liu's "recenter"
      s=boxsize
      findmaxstar,imagei,x2,y2,/silent

      recenter,ref,x1,y1,xc1,yc1,/nodisp,/silent
      recenter,imagei,x2,y2,xc2,yc2,/nodisp,/silent
      shifts[*,n2]=[xc1,yc1]-[xc2,yc2]
   endif
;-----------------------------------------------------------------
   if strc(method) eq "H" then begin ; use FITS header offsets
      shifts[*,n2] = [$
                     sxpar(headers[*,i],fitsoffsetx)*fitsoffsetscale*(-1),$
                     sxpar(headers[*,i],fitsoffsety)*fitsoffsetscale $
                     ] - ref_offset
      
   endif
   


   if not(keyword_set(silent)) then print,strc(shifts[*,n2])
endfor 

end
