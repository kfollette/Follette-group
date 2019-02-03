;+
; NAME: pupilfft
; PURPOSE:
;  Calculates a PSF given a pupil image, by FFT.  Calculates pixel scale.
;
; INPUTS:
; P       ; the pupil image (2D array)
; Pw      ; the pupil image width (m)
; lam     ; the central wavelength (m)
; Nfft    ; linear size of fft, e.g. 2048
; Imw     ; desired image width (pixels)
;
; INPUT KEYWORDS:
; Phase   ; a 2D array, of the same size as P, which is the wavefront phase in meters      
;
; OUTPUTS:
;  returns the PSF   
;  as     ;  1D array, the pixel scale of the image (arcsecs).  Platescale is as[1]-as[0].
;
; EXAMPLE:
;  I = pupilfft(P, 6.5, .765d-6, 2048, 512, as)
;
; HISTORY:
;  Written 2009-02-20 as pupilimage by Jared Males, jrmales@email.arizona.edu
;  Updated 2012-11-06 renamed to pupilfft, documentation updated. (Jared Males)
;-
function pupilfft, P, Pw, lam, Nfft, Imw, as, phase=phase, subcoron=subcoron, complexamp=complexamp, radians=radians

N = size(P)

N=N(1)

;First pad out to the desired size (2^n is best)
pad = complexarr(Nfft, Nfft)

_P  = P

if(n_elements(phase) eq n_elements(_P)) then begin
   ;pad = complex(pad, pad)
   arg = 1.
   if(~keyword_set(radians)) then begin
      arg = (2.d*!dpi/lam)
   endif
   _P = P*exp(arg*phase*complex(0,1))
endif

if(keyword_set(subcoron)) then begin
;_P = _P-P
;exitpup = primary;
;     exitpup.E = primary.E - sum(sum(primary.E.*primary.A))/sum(sum(primary.A.^2))*primary.A;

_P = _P - total(_P*P)/total(P*P)*P

endif

pad[0:N-1, 0:N-1] = _P

if ~keyword_set(complexamp) then begin
   ;Compute fft, and square the abs value
   J = abs( fft(pad,-1, /center) )^2
endif else begin
   J = fft(pad,-1, /center)
endelse

J = shift(J, -1, -1)

;Cut out the desired portion
I = J[Nfft/2-Imw/2:Nfft/2+Imw/2-1,Nfft/2-Imw/2:Nfft/2+Imw/2-1]




as = (double(N)/Nfft) * 206265. * (lam/Pw)


return, I

end
