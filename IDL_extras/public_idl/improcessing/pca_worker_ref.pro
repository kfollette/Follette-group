pro pca_worker_ref, psfsub, adi_rims, sdi_rims, derot, nmodes, mindq, err=err, msims=msims, $
                    nocovar=nocovar, silent=silent, indmean=indmean, modelrims=modelrims, $
                    modelsub=modelsub, imno=imno, maxdq=maxdq, refonly=refonly, dqisdn=dqisdn, $
                    regmedsub=regmedsub, refderot=refderot, adi_mask=adi_mask, ref_mask=ref_mask

;-------------------------------------------------------------------  
;Get the number of pixels and number of images
dim1 = (size(adi_rims))[1]
nims = (size(adi_rims))[2]



sdi_dim1 = (size(sdi_rims))[1]
sdi_nims = (size(sdi_rims))[2]

;-------------------------------------------------------------------  
;Combine the two image groups
rims = [[adi_rims], [sdi_rims]]

if(n_elements(adi_mask) gt 1 and n_elements(ref_mask) gt 1) then begin
   mask = [[adi_mask], [ref_mask]]
endif else begin
   mask = 0;
endelse

;-------------------------------------------------------------------  
;Decide if these are doubles or floats
dodouble = 0
if(size(rims, /type) eq 5) then dodouble = 1

;-------------------------------------------------------------------  
;Calculate the covariance matrix if desired
if(~keyword_set(nocovar)) then begin

   pca_covarmat, err, rims , /meansub, mask=mask

   ;for i=0,nims-1 do adi_rims[*,i] = adi_rims[*,i] - mean(adi_rims[*,i])
   adi_rims = rims[*,0:nims-1]
   
   if(arg_present(msims)) then msims = rims

endif

;-------------------------------------------------------------------  
;Allocate the psfsub storage
if(dodouble) then begin
   psfsub = dblarr(dim1, nims, n_elements(nmodes))
endif else begin
   psfsub = fltarr(dim1, nims, n_elements(nmodes))
endelse

;-------------------------------------------------------------------  
;Setup forward modeling if desired
domodel = 0
if(n_elements(modelrims) gt 1 and arg_present(modelsub)) then begin

   domodel = 1

   if(dodouble) then begin
      modelsub = dblarr(dim1, nims, n_elements(nmodes))
   endif else begin
      modelsub = fltarr(dim1, nims, n_elements(nmodes))
   endelse
   
endif

;-------------------------------------------------------------------
;Region median subtraction
doregmedsub = 0
if(keyword_set(regmedsub)) then doregmedsub = 1

;-------------------------------------------------------------------  
;Default maxdq is huge
if(n_elements(maxdq) lt 1) then maxdq = 1e7

;-------------------------------------------------------------------  
;Massage derot so angle subtraction works
roff = requad_angles(derot)

;-------------------------------------------------------------------  
;Setup for image exclusion, rather than rotation exclusion
if(keyword_set(dqisdn)) then imnum=indgen(nims)

;-------------------------------------------------------------------  
;Allocate klims now
actnmodes = max(nmodes)
klims = fltarr(dim1, actnmodes)

;-------------------------------------------------------------------  
;Setup to do a specific image only, or all images by default
if(n_elements(imno) eq 1) then begin
   i0 = imno
   i1 = imno
endif else begin
   i0 = 0
   i1 = nims-1
endelse

;-------------------------------------------------------------------
;Handle no rotation exclusion intelligently
norotmask=0
klims_done = 0
if(mindq eq 0 and maxdq eq 1e5) then begin
   norotmask = 1
endif

;-------------------------------------------------------------------
;Handle reference images only
derot_scale = 1.
if(keyword_set(refonly)) then begin
   derot_scale = 0.
   if(mindq eq 0) then mindq = 0.0001
   
   norotmask = 1 ;this is always true if /refonly is set
   
   print, 'derot_scale  ', derot_scale
endif


;-------------------------------------------------------------------
;Handle derot angles for the reference images
if(n_elements(refderot) eq sdi_nims) then begin
   refdang = refderot ;this means reference images have signal too
endif else begin
   refdang = fltarr(sdi_nims)+1e6 ;add pure reference images with huge dang.
endelse

;-------------------------------------------------------------------

;####################################################################
;    ACTUAL WORK STARTS HERE
;####################################################################

for i=i0, i1 do begin

   ;KL-images are only calculated every time if a rotation mask is applied
   if( norotmask eq 0 or klims_done eq 0 ) then begin
   
      if(keyword_set(dqisdn)) then begin
         dimnum = [abs(imnum - imnum[i])*derot_scale, refdang];fltarr(sdi_nims)+1e6];add sdi images with huge dang.
         idx = where(dimnum ge mindq)
      
      endif else begin
         dang = [abs(angsub(roff,roff[i])) * derot_scale, refdang];fltarr(sdi_nims)+1e6] ;add sdi images with huge dang.
         idx = where(dang ge mindq and dang le maxdq)
      endelse
      
      terr =  err[*, idx]
      terr =  terr[idx, *]
   
      tims = rims[*, idx]

      pca_klims, klims, terr, tims, actnmodes, /silent
      klims_done = 1
   endif
   
   if(~keyword_set(silent)) then begin
      status = 'pca_worker: performing PCA for image ' + strcompress( string(i + 1) + ' / ' + string(nims), /rem) $
      + ' with ' + strcompress(string(n_elements(idx)), /rem) + ' R images.'
      statusline, status
   endif
   
   cfs = dblarr(actnmodes)
   
   if(domodel) then begin
      cfs_mod = dblarr(actnmodes)
      
      for j=0, actnmodes-1 do begin
      
         cfs[j] = klims[*,j]##transpose(adi_rims[*,i])
         cfs_mod[j] = klims[*,j]##transpose(modelrims[*,i])
         
      endfor
      
   endif else begin
   
      for j=0, actnmodes-1 do cfs[j] = klims[*,j]##transpose(adi_rims[*,i])
   
   
   endelse
   
   if(dodouble) then begin
      psf = dblarr(dim1)
      if(domodel) then psfmod = dblarr(dim1)
   endif else begin
      psf = fltarr(dim1)
      if(domodel) then psfmod = fltarr(dim1)
   endelse
   
   j = 0
   for k=0, n_elements(nmodes)-1 do begin
   
      donmodes = nmodes[k]
      
      if donmodes gt actnmodes then begin
         donmodes = actnmodes
      endif
      
      if(domodel) then begin
      
         for j=j, donmodes-1 do begin
            psf = psf + cfs[j]*klims[*,j]
            psfmod = psfmod + cfs_mod[j]*klims[*,j]
         endfor
         
         newim = adi_rims[*,i] - psf
         if(doregmedsub) then psfsub[*,i, k] = newim-median(newim)
         
         newim = modelrims[*,i] - psfmod
         if(doregmedsub) then modelsub[*,i,k] = newim-median(newim)
      
      endif else begin
      
         for j=j, donmodes-1 do psf = psf + cfs[j]*klims[*,j]
        
         if(n_elements(adi_mask) gt 1) then begin
            newim = adi_rims[*,i] - psf*adi_mask[*,i]
         endif else begin
            newim = adi_rims[*,i] - psf
         endelse
         
         psfsub[*,i, k] = newim
         if(doregmedsub) then psfsub[*,i, k] = psfsub[*,i, k] - median(psfsub[*,i, k])
;        
      endelse
      
   endfor
      
   
endfor

statusline, /clear

end