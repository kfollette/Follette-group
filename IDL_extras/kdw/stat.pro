pro Stat, data, out, lower=lower, upper=upper, nolabel=nolabel, $
          nobadval = nobadval, $
          silent=silent 

;+
; PURPOSE:
; gives median, mean, min, max, and std deviation for a given image
;
; NOTES:
; able to specify lower and upper boundaries to image
; program uses IDL's median function so if there are an even number of pixels,
;   it will take the value just above the median location
;
; checked against IRAF "imstat" and everything looks good
;   but runs somewhat slower than IRAF imstat
;
; Now this ignores any NaN values.
;
; INPUTS
;   data     the data (any IDL array)
;
; OUTPUTS
;   out(0) = number of pixels
;   out(1) = mean
;   out(2) = median
;   out(3) = standard deviation
;   out(4) = minimum pixel
;   out(5) = maximum pixel
;
; KEYWORD PARAMETERS
;   lower    lower limit of data to find stats
;   upper    lower limit of data to find stats
;   nolabel  output results w/o the column headings
;   nobadval exclude BADVAL pixels from calculation of min pixel value
;   silent
;
; HISTORY 
; Written by M. Liu (UCB): 06/22/94
; 05/30/97 (MCL): added /nobadval
; 2002-12-05: Added check to ignore NaNs  - M. Perrin
;
; Please send comments/questions to <mperrin@astro.berkeley.edu>
;-


; return to program calling this one if it bombs
on_error,2
BADVAL = -1e6

if n_params() eq 0 or n_elements(data) eq 0 then begin
    print, 'pro Stat,data,(out),[lower=],[upper=],[nolabel],[nobadval],[silent]'
   return
endif


if keyword_set(lower) eq 0 then begin
    if keyword_set(nobadval) then begin
        llim = min(data(where(data ne BADVAL))) 
        lower = 1
    endif else $
      llim = min(data)
endif else $
  llim=lower
if keyword_set(upper) eq 0 then ulim = max(data) else ulim=upper
; if keyword_set(log) eq 1 then openw,1,log

if keyword_set(silent) eq 0 and keyword_set(nolabel) eq 0 then $
	print, format = '(6(A10,"  "))', $
	"NPIX","MEAN","MEDIAN","STDDEV","MIN","MAX"

; if keyword_set(log) eq 1 then begin
; 	openw,unit0,log,/get_lun
; 	printf,unit0,format = '(5(A8,"  "))', $
; 	   "MEAN","MEDIAN","STDDEV","MIN","MAX"
; endif

; if upper/lower limit(s) set
if ((keyword_set(lower) eq 1) or (keyword_set(upper) eq 1)) then begin
;   print,'limits:  lower=',llim,' upper=',ulim
   w = where(data ge llim and data le ulim and finite(data))
   if (w(0) eq (-1)) then begin
	print,'* limits outside domain of image! *'
	return
   endif
   image  = data(w)
endif else $
   image = data(where(finite(data)))  ; ignore NaN's

; calculate stats
npix = n_elements(image)
mean = total(image)/npix
if (npix gt 1) then $
  sd = stdev(image) $
else begin
    message, 'only one element - not calculating std. deviation', /info
    sd = 0.0
end

med = median(image)
min = min(image)
max = max(image)

; print the results
if keyword_set(silent) eq 0 then $
;  print, format = '(6(G10.3,"  "))', npix, mean, med, sd, min, max
	print,format = '(6(E10.3,"  "))',npix,mean,med,sd,min,max

; dump results to output array
if (n_params() gt 1) then begin
	out = fltarr(6)
	out(0) = npix  &  out(1) = mean  &  out(2) = med  
	out(3) = sd    &  out(4) = min   &  out(5) = max
endif

; if keyword_set(log) eq 1 then begin
;	printf,1,format = '6(F8.0,"  ")', $
;	   npix,mean,med,sd,min,max
; endif

return
end


