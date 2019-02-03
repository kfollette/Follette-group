function airy_1stmin, eps1
;Calculates the location of the 1st Airy minimum, in lambda/D units
;
; Author: Jared R. Males, University of Arizona
; Date: 23 Sep 2009

common airysolve, eps
eps = eps1

return, zbrent(0.1, 3.9, func_name="airy_solver_min")/ !pi

end
