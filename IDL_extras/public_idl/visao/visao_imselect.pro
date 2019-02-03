;+
; NAME: 
;   visao_imselect
;
; DESCRIPTION:
;   selects images from a subdirectory based on various criteria.  can do so interactively.
;
; INPUTS:
;   none
; 
; INPUT KEYWORDS:
;   region  : 2 or 4 vector defining area of image to extract by either its center(2) or corners(4)
;   regsize : if region is a 2-vector, then this 2-vector defines the size of the image around the center.
;   wfethresh :
;   gimbrad  :
;   vfw3n    :
;   interact : if set, images are shown one at a time and the user is asked if it should be used.
;   prefix   :
;   darksub  :
;   useatv   :  if set, atv is used instead of ds9 for image display when interactive.  
;
; OUTPUTS:
;   fnames  :  the names of the selected files

; OUTPUT KEYWORDS:
;   ims    :  the selected images.  careful - don't use this if you have too many images
;   imtypes
;   aoloopst
;   rotoff
;   exptime
;   dateobs
;   avgwfe
;   vgimxpos
;   vgimypos
;   dark
;
; MODIFICATION HISTORY:
;  Written 2013/01/05-ish by Jared Males (jrmales@email.arizona.edu)
;  2013/04/14 added ds9 as default interactive display (Jared Males)
;
; BUGS/WISH LIST:
;
;-
pro visao_imselect, fnames,  region=region, regsize=regsize, wfethresh=wfethresh, gimbrad=gimbrad, vfw3n=vfw3n,$
                  interact=interact, prefix=prefix, darksub=darksub, $
                  ims = ims, imtypes=imtypes, aoloopst=aoloopst, rotoff=rotoff, exptime=exptime, dateobs=dateobs,$
                  avgwfe=avgwfe, vgimxpos=vgimxpos, vgimypos=vgimypos,  dark=dark, ignorehead=ignorehead,$
                  usefnames=usefnames, subdir=subdir, usm=usm, sm=sm

                   
common atv_state, state
 

if(n_elements(region) eq 2) then begin

   if(n_elements(regsize) eq 2) then begin
   
      nreg = [region[0]-regsize[0]*.5, region[0]+regsize[0]*.5-1, region[1]-regsize[1]*.5, region[1]+regsize[1]*.5-1]
 
   endif else begin
   
      message, 'Must have 2 element regsize for 2 element region'
   endelse
endif else begin

   if(n_elements(region) eq 4) then begin
      nreg = region
   endif
   
endelse


visao_getimtypes, fnames, imtypes, aoloopst=aoloopst, rotoff=rotoff, avgwfe=avgwfe, vgimxpos=gimx, vgimypos=gimy, vfw3posn=vfw3posn, prefix=prefix, exptime=exptime, dateobs=dateobs, usefnames=keyword_set(usefnames), subdir=subdir


print, 'Found ', n_elements(fnames), ' images.'

;--- setup wfe threshold ---
if(n_elements(wfethresh) eq 0) then begin
   wfethresh = max(avgwfe)
endif else begin
   
   ;If wfethresh is 0<=wfethresh<=1 then it is a percentile
   if(wfethresh le 1) then begin
   
      idx = sort(avgwfe)
      
      wfethresh = (avgwfe[idx])[wfethresh*n_elements(idx)]
      ;wfethresh = median(avgwfe)
   endif
   
endelse

;--- setup gimbal threshold ---
medgx = median(gimx)
medgy = median(gimy)

gimr = sqrt((gimx-medgx)^2 + (gimy-medgy)^2)*488.

if(n_elements(gimbrad) eq 0) then gimbrad = 1e6


if(n_elements(vfw3n) gt 0) then begin

   idx = where(imtypes eq 2 or (imtypes eq 0 and aoloopst eq 1 and avgwfe le wfethresh and gimr le gimbrad and vfw3posn eq vfw3n))

endif else begin

   if(~keyword_set(ignorehead)) then begin
      idx = where(imtypes eq 2 or (imtypes eq 0 and aoloopst eq 1 and avgwfe le wfethresh and gimr le gimbrad))
   endif else begin
      idx = where(imtypes eq 2 or imtypes eq 0)
   endelse
endelse
   
   
   
print, 'Selected ', n_elements(idx), ' images'

if(arg_present(ims)) then begin

   print, 'loading images . . .'
   if(n_elements(nreg) eq 4) then begin

      visao_getimtypes, fnames, imtypes, aoloopst=aoloopst, rotoff=rotoff, avgwfe=avgwfe, vgimxpos=gimxpos, vgimypos=gimypos, ims=ims, region=nreg, goodims=idx, prefix=prefix,exptime=exptime, dateobs=dateobs, subdir=subdir

   endif else begin

      visao_getimtypes, fnames, imtypes, aoloopst=aoloopst, rotoff=rotoff, avgwfe=avgwfe, vgimxpos=vgimxpos, vgimypos=vgimypos, ims=ims, goodims=idx, prefix=prefix,exptime=exptime, dateobs=dateobs, subdir=subdir

   endelse

endif else begin

   fnames = fnames[idx]
   imtypes = imtypes[idx]
   rotoff = rotoff[idx]
   aoloopst=aoloopst[idx]
   avgwfe=avgwfe[idx]
   vgimxpos=gimx[idx]
   vgimypos=gimy[idx]

endelse

if(keyword_set(interact)) then begin
   sel = intarr(n_elements(fnames)) + 1
   
   if(keyword_set(useatv)) then atv
   for i=0,n_elements(fnames)-1 do begin

      ;state.zoom_level = 0
      
      
      if(arg_present(ims)) then begin
         tim = ims[*,*,i]
         if(n_elements(usm) eq 1) then tim = tim - filter_image(tim, fwhm_g=usm)
         if(n_elements(sm) eq 1) then tim = filter_image(tim, fwhm_g=sm)
         if(keyword_set(useatv)) then begin
            atv, tim
         endif else begin
            ds9, tim
         endelse
      endif else begin
         tim = mrdfits(fnames[i])
         if(n_elements(usm) eq 1) then tim = tim - filter_image(tim, fwhm_g=usm)
         if(n_elements(sm) eq 1) then tim = filter_image(tim, fwhm_g=sm)
         if(keyword_set(useatv)) then begin
            atv, tim
         endif else begin
            ds9, tim
         endelse
      endelse
      
      prog = strcompress(string(i+1) + '/' + string(n_elements(fnames)), /rem)
      if(imtypes[i] eq 0) then prog =  'SCIENCE ' + prog
      if(imtypes[i] eq 1) then prog =  'ACQ ' + prog
      if(imtypes[i] eq 2) then prog =  'DARK ' + prog
      if(imtypes[i] eq 3) then prog =  'SKY ' + prog
      if(imtypes[i] eq 4) then prog =  'FLAT ' + prog

      prog = prog + ' ' + strcompress(string(avgwfe[i]),/rem) +' nm rms phase.'
      
      print, prog
      
      ok_str = 'ok'

      read, ok_str, prompt='Good image? [y/n] (y) '
      if(ok_str eq 'n' or ok_str eq 'N' ) then sel[i] = 0
      if(ok_str eq 'b' or ok_str eq 'B' ) then i = i - 2
      
   endfor

   idx = where(sel eq 1)
   
   help, idx
   fnames=fnames[idx]
   if(arg_present(imtypes))  then imtypes=imtypes[idx]
   if(arg_present(aoloopst)) then aoloopst = aoloopst[idx]
   if(arg_present(rotoff)) then rotoff = rotoff[idx]
   if(arg_present(exptime)) then exptime = exptime[idx]
   if(arg_present(dateobs)) then dateobs = dateobs[idx]
   if(arg_present(avgwfe)) then avgwfe=avgwfe[idx]
   
   if(arg_present(vgimxpos)) then vgimxpos=vgimxpos[idx]
   if(arg_present(vgimypos)) then vgimypos=vgimypos[idx]
   
   if(arg_present(ims)) then ims = ims[*,*,idx]
   
endif



if (keyword_set(darksub) and arg_present(ims)) then begin
   print, 'Subtracting darks'
   visao_cube_darksub, ims, imtypes, aoloopst=aoloopst, rotoff=rotoff, dark=dark
endif


end
