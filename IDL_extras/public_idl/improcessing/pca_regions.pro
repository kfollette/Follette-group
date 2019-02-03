;+
; NAME: pca_regions
;
; PURPOSE: 
;  Conduct PCA PSF subtraction on a set of images, breaking up the analysis into search regions
;
; Description:
;  Conducts PCA (the KLIP algorithm) on a cube of images, which are assumed to be registered, dark subtracted, etc.  
;  Search regions are specified in terms of radial and angular size.  
;
; INPUTS:
;  ims    :  the fully prepared (registered, dark subtracted, etc.) image cube
;  derot  :  vector of counter-clockwise rotation angles, one for each image in ims
;  mindpx :  the minimum rotation, nominally in pixels, at the inner edge of the search region.
;  regdr  :  region delta-r, the radial size of the search regions in pixels
;  regdq  :  region delta-q, the angular size of the search regions in degrees
;  nmodes :  a scalar or vector specifying the number of modes. if not set then nmodes = 0.3*number_of_images
;
; INPUT KEYWORDS:
;  minrad    :  the minimum radius to analyze (default 0)
;  maxrad    :  the maximum radius to analyze (default is .5*max[dim1, dim2])
;  minang    :  the zero point for the angular search regions (default 0)
;  dpxisdq   :  treat mindpx as an angle, so that minimum angular rotation is constant with radius
;  dqisdn    :  treat mindpx as a number of images to mask, rather than either pixels or angle
;  meancomb  :  mean combine final image, rather than median (the default)
;  sigma     :  sigma clip mean-combine the final image, with sigma specifying how many sigma to clip
;  regmedsub :  subtract the median of each region after PSF subtraction
;  model_ims :  a cube of images to use for forward modeling.  should be identical in size to ims.
;  ref_ims   :  a cube of reference images to use as part of the basis set.
;  refonly   :  if set, then only the reference images are used in determining the PSF
;  refderot  :  derotation angles for the reference images, use if ADI exclusion is desired
;  core_rads :  a 2 element vector, specifying the min and max radius of the core.  if set, then these pixels will be 
;               included in all search regions regardless of their size and location.
;  deltar    :  spacing of search regions.  by default this is regdr.  the amount of each search region kept in the
;            :  final image will be deltar.
;  silent    :  if set then no status messages are printed.
;  mask      :  a binary mask, either single image or a cube, whose pixels are set to excluded from the reduction.  note that these 
;               pixels are assumed to be 0 in the input ims.  If a cube, each image can have a different mask.  Note also that for 
;               best results each image in ims should be high-pass filtered, such as by radial profile subtraction
;  maxdq     :  the maximum rotation to include.  images which have rotated more than this are not used in the psf.
;
; OUTPUTS:
;  finim    : the final derotated, combined image(s)
;
; OUTPUT KEYWORDS:
;  psfsub      :  the individual PSF subtracted images, derotated.
;  model_finim :  the result of forward modeling using model_ims
;  fitsfile    :  base name of the fits file to write output too.  actual name will unique-ified with xxxx. fits file
;                 includes a header with all settings.
;
; DEPENDENCIES:
;   This procedure depends on the following code from JRM:
;      get_cubedims.pro
;      rarr.pro 
;      pca_worker.pro 
;      pca_worker_ref.pro 
;      derotcomb_cube.pro 
;  This procedure depends on the following external code:
;      mwrfits.pro 
;      fxaddpar.pro 
;
; TODO:
;   - investigate the way mean subtraction is handled in the model and reference cases
;   - implement forward modeling with reference images
;   - handle case where mindpx = 0 without recalculating and decomposing covar. mat. each time
;   
;
; HISTORY:
;  2013-03-01: Written by Jared Males, jrmales@email.arizona.edu
;  2014-04-15: Added reference image capabilities. Jared Males
;  2014-04-30: Updated documentation
;  2015-04-06: removed indmean keyword (must always be true), updated masking.  updated docs.
;-
pro pca_regions, finim, ims, derot, mindpx, regdr, regdq, nmodes, minrad=minrad, maxrad=maxrad, minang=minang, $
                       dpxisdq=dpxisdq, dqisdn=dqisdn,$
                       meancomb=meancomb, sigma=sigma,  regmedsub=regmedsub, psfsub=psfsub, $
                       model_ims=model_ims, model_finim=model_finim, $
                       fitsfiles=fitsfile, ref_ims=ref_ims, refonly=refonly, refderot=refderot, $
                       core_rads=core_rads, deltar=deltar, silent=silent, mask=mask, ref_mask=ref_mask, $
                       maxdq=maxdq
                       
;indmean=indmean,

if ~keyword_set(silent) then print, "Start: ", systime()    

dodouble = 0
if(size(ims, /type) eq 5) then dodouble = 1
                       
get_cubedims, ims, dim1, dim2, nims

;Check if a PSF model has been supplied
domodel = 0
if(n_elements(model_ims) gt 1) then begin
   domodel = 1
endif 

;Check if this is an ASDI reduction
doref = 0
if(n_elements(ref_ims) gt 1) then begin
   doref = 1
   get_cubedims, ref_ims, ref_dim1, ref_dim2, ref_nims
endif

if(n_elements(deltar) ne 1) then deltar=regdr

;--------------------------------------------------------------------
;Check if a mask has been supplied
if(n_elements(mask) lt dim1*dim2) then begin
   domask = 0;_mask = dblarr(dim1,dim2) + 1.
endif else begin
   get_cubedims, mask, maskdim1, maskdim2, masknims
   if(masknims eq 1) then begin
      domask = 1
   endif else begin
      domask = 2; changing mask
   endelse
endelse

;and Check if reference mask has been applied
domask_ref = 0
if(doref eq 1) then begin
   if(n_elements(ref_mask) lt dim1*dim2) then begin
      domask_ref = 0;
   endif else begin
      get_cubedims, ref_mask, ref_maskdim1, ref_maskdim2, ref_masknims
      if(ref_masknims eq 1) then begin
         domask_ref = 1
      endif else begin
         domask_ref = 2; changing mask
      endelse
   endelse 
endif

;-------------------------------------------------------------------
;Calculate the average image,and subtract it, unless we're doing individual image means
; if(~keyword_set(indmean)) then begin
;    mnpsf = total(ims, 3)/double(nims)
; 
;    for i=0, nims-1 do ims[*,*,i] = ims[*,*,i] - mnpsf
; 
;    if(domodel) then begin
;       modpsf = total(model_ims, 3)/double(nims)
;       for i=0, nims-1 do model_ims[*,*,i] = model_ims[*,*,i] - modpsf
;    endif
;    
;    if(doref) then begin
;       refpsf = total(ref_ims, 3)/double(ref_nims)
;       for i=0, ref_nims-1 do ref_ims[*,*,i] = ref_ims[*,*,i] - refpsf
;    endif
;    
; endif

;-------------------------------------------------------------------


;allocate the psf subtracted image arrays
if(dodouble) then begin
   psfsub = dblarr(dim1*dim2, nims, n_elements(nmodes))
   if(domodel) then model_psfsub = dblarr(dim1*dim2, nims, n_elements(nmodes))
endif else begin
   psfsub = fltarr(dim1*dim2, nims, n_elements(nmodes))
   if(domodel) then model_psfsub = fltarr(dim1*dim2, nims, n_elements(nmodes))
endelse

;-------------------------------------------------------------------

;Check defaults and establish search region parameters

if(n_elements(maxrad) ne 1) then maxrad = max([dim1, dim2])*0.5
if(n_elements(minrad) lt 1) then minrad = 0.
if(n_elements(minang) lt 1) then minang = 0

if(n_elements(regdr) eq 1) then begin
   drconst = 1
   nregrs = floor(maxrad/regdr)
   nregqs = ceil(360./regdq)
endif else begin
   drconst = 0
   nregrs = n_elements(regdr)-1
endelse

if(n_elements(nmodes) eq 0) then nmodes = floor(nims*.3)

ndiffmodes = n_elements(nmodes)

;-------------------------------------------------------------------

;The radius array for picking regions
r=rarr(dim1, dim2, x, y, theta=q, /pixels)

;Convert to degrees
q = q/!dtor
idx = where(q lt minang)
q[idx] = q[idx] + 360.

;-------------------------------------------------------------------


;reform to vectors
ims = reform( ims, dim1*dim2, nims, /overwrite)
      
if(domodel) then model_ims = reform(model_ims, dim1*dim2, nims, /overwrite)

if(doref) then ref_ims = reform(ref_ims, dim1*dim2, ref_nims, /over)

if(domask eq 2) then mask = reform(mask, dim1*dim2, nims, /overwrite)

if(domask_ref eq 2) then ref_mask = reform(ref_mask, dim1*dim2, nims, /over)

;-------------------------------------------------------------------

;Check if core pixels are included in each region
if n_elements(core_rads) eq 2 then begin
   core_min = core_rads[0]
   core_max = core_rads[1]
endif else begin
   core_min = -1.
   core_max = -1.
endelse

;-------------------------------------------------------------------

;####################################################################
;    ACTUAL WORK STARTS HERE
;####################################################################

;Now loop through the regions      
for i = 0d, nregrs do begin ;radial 

   if(drconst) then begin
      rmin = double(minrad) + i*double(regdr)
      rmax = rmin + double(deltar)
   endif else begin
   
      rmin = minrad[i] 
      rmax = minrad[i] + regdr[i]
   endelse
   
   if(rmin ge maxrad) then break
   
   if(~keyword_set(dpxisdq) and ~keyword_set(dqisdn)) then begin
      mindq = atan(double(mindpx)/rmin)/!dtor
      if(~finite(mindq)) then mindq = 0
   endif else begin
      mindq = mindpx
   endelse
   
   if(drconst eq 0) then begin
      nregqs = floor(360./regdq[i])
   endif
   
   for j=0, nregqs-1 do begin ;angles
   
      if(drconst eq 1) then begin
         qmin = minang + j*regdq
         qmax = qmin + regdq
      endif else begin
         qmin = minang + j*regdq[i]
         qmax = qmin + regdq[i]
      endelse
      
      if ~keyword_set(silent) then begin
         print, 'pca_regions: ' + strcompress(string(rmin),/rem) + ' <= r < ' + strcompress(string(rmax),/rem) + $
                ' ' + strcompress(string(qmin),/rem) + ' <= q < ' + strcompress(string(qmax),/rem)
      endif
      
      if(domask eq 1) then begin
         idx = where( ((r ge rmin and r lt rmax and q ge qmin and q lt qmax) or (r ge core_min and r lt core_max)) and mask ne 0);
      endif else begin
         idx = where( ((r ge rmin and r lt rmax and q ge qmin and q lt qmax) or (r ge core_min and r lt core_max)));
      endelse
      
      if(idx[0] eq -1) then continue
      
      rims = ims[idx, *]
      
      if(domodel) then modelrims = model_ims[idx, *]

      if(doref) then refrims = ref_ims[idx,*]
      
      if(domask eq 2) then begin
         mask_rims = mask[idx, *]
      endif else begin
         mask_rims =0
      endelse
      
      if(domask_ref eq 2) then begin
         ref_mask_rims = ref_mask[idx, *]
      endif else begin
         ref_mask_rims =0
      endelse
      
      if(domodel) then begin
         pca_worker, rpsfsub, rims, derot, nmodes, mindq, $ 
                          dqisdn=keyword_set(dqisdn), modelrims=modelrims, modelsub=modelsub, regmedsub=keyword_set(regmedsub), silent=keyword_set(silent)
      endif else begin
         if(doref) then begin
            pca_worker_ref, rpsfsub, rims, refrims, derot, nmodes, mindq,  $
                       refonly=keyword_set(refonly),dqisdn=keyword_set(dqisdn), regmedsub=keyword_set(regmedsub),$
                             silent=keyword_set(silent), refderot=refderot, adi_mask=mask_rims, ref_mask=ref_mask_rims
                             
         endif else begin
            pca_worker, rpsfsub, rims, derot, nmodes, mindq, $
                     dqisdn=keyword_set(dqisdn), regmedsub=keyword_set(regmedsub), silent=keyword_set(silent), maxdq=maxdq, $
                        mask=mask_rims
         endelse
      endelse
      
      
      psfsub[idx, *,*] = rpsfsub
      
      if(domodel) then model_psfsub[idx, *, *] = modelsub
      
   endfor ;j...nregqs
endfor ;i...nregrs
   
;-------------------------------------------------------------------  
;Post processing   
;Reform back to 2D images
ims = reform(ims, dim1, dim2, nims, /overwrite)

if(doref) then ref_ims = reform(ref_ims, dim1, dim2, ref_nims, /over)

if(domask eq 2) then begin
   mask = reform(mask, dim1, dim2, nims, /overwrite)
endif 

if(domask_ref eq 2) then begin
   ref_mask = reform(ref_mask, dim1, dim2, nims, /overwrite)
endif 

psfsub = reform(psfsub, dim1, dim2, nims, n_elements(nmodes), /overwrite)

;Allocate the final image array
if(dodouble) then begin
   finim = dblarr(dim1, dim2, n_elements(nmodes))
endif else begin
   finim = fltarr(dim1, dim2, n_elements(nmodes))
endelse

;Now derot each individual combination


for k=0, n_elements(nmodes) -1 do begin
   rotpsfsub = psfsub[*,*,*,k]
   
   derotcomb_cube, finimt, rotpsfsub , derot, meancomb=keyword_set(meancomb), sigma=sigma, $
                     silent=keyword_set(silent), mask=mask
   finim[*,*,k] = finimt
   
;    p = -1
;    print, finimt[p]
endfor

;-------------------------------------------------------------------

if(domodel) then begin

   model_ims = reform(model_ims, dim1, dim2, nims, /overwrite)

   model_finim = dblarr(dim1, dim2, n_elements(nmodes))

   model_psfsub = reform(model_psfsub, dim1, dim2, nims, n_elements(nmodes), /overwrite)

   for k=0, n_elements(nmodes) -1 do begin

      derotcomb_cube, finimt, model_psfsub[*,*,*,k], derot, meancomb=keyword_set(meancomb), sigma=sigma
      model_finim[*,*,k] = finimt
   endfor
   
endif


if (n_elements(fitsfile) gt 0) then begin

   mkhdr, header, finim
   fxaddpar, header, 'TIME-UTC', systime(/utc), 'time this file was written' 
   fxaddpar, header, 'MINDPX', string(mindpx)
   fxaddpar, header, 'REGDR', string(regdr)
   fxaddpar, header, 'REGDQ', string(regdq)
   fxaddpar, header, 'NMODES', string(nmodes,format='('+strtrim(n_elements(nmodes))+'(I0,:,","))')
   fxaddpar, header ,'MINRAD', string(minrad)
   fxaddpar, header, 'MAXRAD', string(maxrad)
   fxaddpar, header, 'MINANG', string(minang)
   fxaddpar, header, 'COREMIN', string(core_min)
   fxaddpar, header, 'COREMAX', string(core_max)
   fxaddpar, header, 'DELTAR', string(deltar)
   fxaddpar, header, 'DPXISDQ', string(keyword_set(dpxisdq))
   fxaddpar, header, 'DQISDN', string(keyword_set(dqisdn))
   fxaddpar, header, 'MEANCOMB', string(keyword_set(meancomb))
   fxaddpar, header, 'SIGMA', string(sigma)
   fxaddpar, header, 'INDMEAN', string(keyword_set(indmean))
   fxaddpar, header, 'REGMEDSU', string(keyword_set(regmedsub))
   fxaddpar, header, 'COMMENT', 'pca_regions.pro output',before = 'MINDPX'
   fxaddpar, header, 'COMMENT', 'Command line parameters:',before = 'MINDPX'
   
   fname = file_unique(fitsfile, '.fits')
   mwrfits, finim, fname, header
   

   ;This needs to re-use the number found for the main fits file.
   if(domodel) then begin
      strput, fname, '_mode', strpos(fname, '.fits')
      fname = strcompress(fname+'l.fits', /rem)
      fxaddpar, header, 'COMMENT', 'This is the model'
      mwrfits, model_finim, fname, header
   endif
endif

if ~keyword_set(silent) then begin

   print, "Stop: ", systime()

   print,string([7B]),format='(1A,$)' ;system bell

   print, 'Basic PCA done'
endif

end


