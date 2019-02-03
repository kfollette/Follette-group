;+
; NAME: visao_getderot
;
; PURPOSE: 
;  Calculate the value of VISAO_DEROT
;
; Description:
;  Calculates VISAO_DEROT, the angle by which to rotate a VisAO image to get North up, East left.
;
; INPUTS:
;  rotoff   :   the value of the rotator offset, from the ROTOFF fits keyword
;
; INPUT KEYWORDS:
;  None.
;
; OUTPUTS:
;  VISAO_DEROT   :  the value by which to rotate CCW.
;
; OUTPUT KEYWORDS:
;  None.
;
; HISTORY:
;  2013-07-13: Written by Jared Males, jrmales@as.arizona.edu
;-
function visao_getderot, rotoff

visao_getcal, NORTH_VISAO

DEROT_VISAO = rotoff + 90. + NORTH_VISAO

return, derot_visao   

end


