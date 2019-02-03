function airy_solver_max, x
;Solver working function for airy_1stmax
;
; Author: Jared R. Males, University of Arizona
; Date: 23 Sep 2009

common airysolve, eps

return, beselj(x, 2) - eps*eps*beselj(eps*x,2)

end
