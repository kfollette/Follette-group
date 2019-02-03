function airy_solver_min, x
;Solver working function for airy_1stmin
;
; Author: Jared R. Males, University of Arizona
; Date: 23 Sep 2009

common airysolve, eps

return, beselj(x, 1) - eps*beselj(eps*x,1)

end
