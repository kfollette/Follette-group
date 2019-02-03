pro pca_worker, psfsub, rims, rotoff, nmodes, mindq, err=err, msims=msims, nocovar=nocovar, $
                silent=silent, modelrims=modelrims, modelsub=modelsub, imno=imno, $
                maxdq=maxdq, dqisdn=dqisdn, regmedsub=regmedsub, mask=mask
;
;  psfsubs  :  the psf-subtracted images
;  rims     :  the input images, formatted as row images
;
;  err      :  the covariance matrix (E_RR in Soummer et al notation), can be either input or output
;  nocovar  :  if set, then the covariance matrix is not re-calculated, and err is an input 
;

;-------------------------------------------------------------------  
;Get the number of pixels and number of images
dim1 = (size(rims))[1]
nims = (size(rims))[2]

;-------------------------------------------------------------------  
;Decide if these are doubles or floats
dodouble = 0
if(size(rims, /type) eq 5) then dodouble = 1

;-------------------------------------------------------------------  
;Calculate the covariance matrix if desired
if(~keyword_set(nocovar)) then begin

   pca_covarmat, err, rims, /meansub, mask=mask ;=(keyword_set(indmean))

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

;    if(keyword_set(indmean)) then begin
;    
;       for i=0, nims-1 do modelrims[*,i] = modelrims[*,i] - mean(modelrims[*,i])
;       
;    endif
endif

;-------------------------------------------------------------------  
;Region median subtraction
doregmedsub = 0
if(keyword_set(regmedsub)) then doregmedsub = 1

;-------------------------------------------------------------------  
;Default maxdq is huge
if(n_elements(maxdq) lt 1) then maxdq = 1e5

;-------------------------------------------------------------------  
;Massage rotoff so angle subtraction works
roff = requad_angles(rotoff)

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

;####################################################################
;    ACTUAL WORK STARTS HERE
;####################################################################

for i=i0, i1 do begin


   ;KL-images are only calculated every time if a rotation mask is applied
   if( norotmask eq 0 or klims_done eq 0 ) then begin
   
      if(keyword_set(dqisdn)) then begin
         idx = where(abs(imnum - imnum[i]) ge mindq)
      endif else begin
         dang = abs(angsub(roff,roff[i]))
         idx = where(dang ge mindq and dang le maxdq)
      endelse
   
      
      terr =  err[*, idx]
      terr =  terr[idx, *]

      tims = rims[*, idx]
 
      pca_klims, klims, terr, tims, actnmodes, /silent
      klims_done = 1
      
      ;if(n_elements(mask gt 1)) then klims = klims*mask
      ;print, klims[-1]
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

         cfs[j] = klims[*,j]##transpose(rims[*,i])
         cfs_mod[j] = klims[*,j]##transpose(modelrims[*,i])
         
      endfor
      
   endif else begin
      for j=0, actnmodes-1 do cfs[j] = klims[*,j]##transpose(rims[*,i])
      
  
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
         
         newim = rims[*,i] - psf
         psfsub[*,i, k] = newim
         if(doregmedsub) then psfsub[*,i, k] = psfsub[*,i, k] - median(psfsub[*,i, k])
         
         newim = modelrims[*,i] - psf;psfmod
         modelsub[*,i,k] = newim;-median(newim)
         if(doregmedsub) then modelsub[*,i,k] = modelsub[*,i,k] - median(modelsub[*,i,k])
         
         
      endif else begin
      
         for j=j, donmodes-1 do psf = psf + cfs[j]*klims[*,j]
        
         if(n_elements(mask) gt 1) then begin
            newim = rims[*,i] - psf*mask[*,i]
         endif else begin
            newim = rims[*,i] - psf
         endelse
         psfsub[*,i, k] = newim;-median(newim)
         if(doregmedsub) then psfsub[*,i, k] = psfsub[*,i, k] - median(psfsub[*,i, k])
         
      endelse
   endfor
endfor

statusline, /clear

end