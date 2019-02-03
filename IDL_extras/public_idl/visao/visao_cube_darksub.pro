;+
; NAME: 
;   visao_cube_darksub
;
; PURPOSE:
;   Dark subtract a cube of VisAO images
;
; DESCRIPTION:
;   Given an image cube and a vector of image types corresponding to the VIMTYPE fits header keyword,
;   calculates the dark and subtracts it from each image.  The dark is formed as the median of all
;   darks in the cube.  On completion, the cube contains only the dark subtracted images.
;
; INPUTS:
;   ims : the image cube, in [dim1, dim2, no_ims] format.
;   imtypes : the types of each image in the cube.  see the VisAO documentation.

; INPUT KEYWORDS:
;   medsub   :  mean subtract each image after darksubtraction
;   meandark :  calculate dark using mean, rather than median (the default).
;   masksat  :  if set, then all values >= 16383 (14bits) will be replaced by nan
;
; OUTPUTS:
;   ims   :  on completion, ims will contain only the dark subtracted images.
;
; OUTPUT KEYWORDS:
;   dark       :  the dark image
;   sdx        :  the indices of the science images
;   [KEYWORD]  :  vectors of various standard visao fits header values can be passed.  on exit the vectors are 
;                 shortened so that they correspond to the dark subtracted images.
;
; MODIFICATION HISTORY:
;  Written 2013/04/12 by Jared Males (jrmales@email.arizona.edu)
;
; BUGS/WISH LIST:
;
;-
pro visao_cube_darksub, ims, imtypes, medsub=medsub, meandark=meandark,  masksat=masksat, dark=dark, sdx=sdx, $           
                         AOLOOPST=AOLOOPST, ROTOFF=ROTOFF, EXPTIME=EXPTIME, AVGWFE=AVGWFE, DATEOBS=DATEOBS,$
                         VFOCPOS=VFOCPOS, VGIMXPOS=VGIMXPOS, VGIMYPOS=VGIMYPOS, FNAMES=FNAMES

dim1 = (size(ims[*,*,0]))[1]
dim2 = (size(ims[*,*,0]))[2]



;----Calculate Dark----;
ddx = where( imtypes eq 2 )
   
print, 'Found ', n_elements(ddx), ' darks.'

if(n_elements(ddx) gt 1) then begin
   if(~keyword_set(meandark)) then begin
      dark = median(ims[*,*,ddx], dim=3);
   endif else begin
      dark = total(ims[*,*,ddx], 2)/double(n_elements(ddx))
   endelse
endif else begin
   dark = ims[*,*,ddx]
endelse


sdx = where( imtypes eq 0 )

nims = n_elements(sdx)
   
print, 'visao_cube_darksub: found ', nims, ' science images.'


ims = ims[*,*,sdx]




for i=0, nims-1 do begin
  
      status = strcompress(string(i+1) + '/' + string(nims), /rem)
      statusline, status, 0
      
      if(keyword_set(masksat)) then begin
         mim = ims[*,*,i]
         mdx = where(mim ge 16383)
         mim[mdx] = float("nan")
         ims[*,*,i] = mim
      endif
      
      ims[*,*,i] = ims[*,*,i] - dark
      
      if(keyword_set(medsub)) then ims[*,*,i] = ims[*,*,i] - median(ims[*,*,i])
endfor

statusline, /clear

imtypes=imtypes[sdx]

if(n_elements(AOLOOPST) gt 0) then AOLOOPST=AOLOOPST[sdx]
if(n_elements(ROTOFF) gt 0) then ROTOFF=ROTOFF[sdx]
if(n_elements(EXPTIME) gt 0) then EXPTIME=EXPTIME[sdx]
if(n_elements(AVGWFE) gt 0) then AVGWFE=AVGWFE[sdx]
if(n_elements(DATEOBS) gt 0) then DATEOBS=DATEOBS[sdx]
if(n_elements(VFOCPOS) gt 0) then VFOCPOS=VFOCPOS[sdx]
if(n_elements(VGIMXPOS) gt 0) then VGIMXPOS=VGIMXPOS[sdx]
if(n_elements(VGIMYPOS) gt 0) then VGIMYPOS=VGIMYPOS[sdx]
if(n_elements(FNAMES) gt 0) then FNAMES=FNAMES[sdx]

end







