;+
; NAME: bin_median
;
; PURPOSE: 
;  Re-bin data by calculating the median and mean in bins
;
; Description:
;  Calculates the median (and mean) of the data in bins.
;
; INPUTS:
;
; INPUT KEYWORDS:
;  none
;
; OUTPUTS:
;    
;
; EXAMPLE:
;  
;
; HISTORY:
;  2008: Written by Jared Males, jrmales@email.arizona.edu
;  2012-04-19: Jared Males added documentation
;-
pro bin_median, x, y, minx, maxx, nbins, binx, binmed, binmean, stdx, stdy, nostdmean=nostdmean, sigma=sigma


binsz = double(maxx-minx)/double(nbins)
;print, binsz

binx = dblarr(nbins)
binmed = dblarr(nbins)
binmean = dblarr(nbins)
stdx = dblarr(nbins)
stdy = dblarr(nbins)

for i=0l, nbins-1 do begin
   minbin = minx + i*binsz
   maxbin = minbin+binsz
   binx[i] = minbin+.5*binsz

   idx = where(x ge minbin and x lt maxbin)

   if(idx[0] eq -1) then begin
      binmed[i] = -99999999.
      binmean[i] = -99999999.
      stdx[i] = -99999999.
      stdy[i] = -99999999.
   endif else begin
      binmed[i] = median(y[idx])
      
      if(n_elements(sigma) gt 1) then begin
         lsqmean, y[idx], sigma[idx],  lsqm, lsqerr
         binmean[i] = lsqm
         stdy[i] = lsqerr
      endif else begin
         binmean[i] = mean(y[idx])
      
         if(keyword_set(nostdmean)) then begin
            stdy[i] = stddev(y[idx])
         endif else begin
            stdy[i] = stddev(y[idx])/double(n_elements(idx))
         endelse
      endelse
      stdx[i] = stddev(x[idx])

   endelse
endfor


end





