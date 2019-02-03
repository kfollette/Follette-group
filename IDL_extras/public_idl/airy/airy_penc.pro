function airy_penc, x, eps
;This calculates the fraction of enclosed power.
;This is different from the formula given by wikipedia, which isn't normalized.
;See: Mahajan, V. N., JOSA 3,4,470 (1986)
;Units of x are lambda/D
;
; Author: Jared R. Males, University of Arizona
; Date: 23 Sep 2009

x1= double(x)*!dpi

JINT = qpint1d('beselj(x,1)*beselj(P*x,1)/x', 0.d, x1, double(eps), /EXPRESSION)

encp = 1 - (beselj(x1,0))^2 - (beselj(x1,1))^2 + (eps^2)*(1 - (beselj(x1,0))^2 - (beselj(x1,1))^2)
encp = encp - 4*eps*JINT
encp = encp/(1-eps^2)

return, encp

end

