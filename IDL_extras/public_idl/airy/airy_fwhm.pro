function airy_fwhm, eps1
;Calculates the FWHM of an obscured Airy peak.
;Units are lambda/D
;
; Author: Jared R. Males, University of Arizona
; Date: 23 Sep 2009

common airysolve, eps
eps = eps1

return, 2.*zbrent(0.1, 1.22, func_name="airy_solver_fwhm")

end
