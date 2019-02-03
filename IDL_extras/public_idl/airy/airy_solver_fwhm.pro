function airy_solver_fwhm, x
;Solver working function for airy_fwhm
;
; Author: Jared R. Males, University of Arizona
; Date: 23 Sep 2009

common airysolve, eps

return, airy(x, eps) - .5

end
