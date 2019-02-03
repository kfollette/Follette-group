function airy_1stmax, eps1
;Calculates the location of the 1st maximum, in lambda/D units.
;
; Author: Jared R. Males, University of Arizona
; Date: 23 Sep 2009

common airysolve, eps
eps = eps1

return, zbrent(3.14, 7., func_name="airy_solver_max")/ !pi

end
