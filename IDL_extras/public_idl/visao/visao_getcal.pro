;+
; NAME: visao_getcal
;
; PURPOSE: 
;  Get the value of various VisAO camera calibrations
;
; Description:
;  Returns hard coded values.
;
; INPUTS:
;  None.
;
; INPUT KEYWORDS:
;  None.
;
; OUTPUTS:
;  NORTH_VISAO   :   the value of the NORTH_VISAO parameter
;  visao_gx      :   the x position of the ghost, relative to the star
;  visao_gy      :   the y position of the ghost, relative to the star
;
; OUTPUT KEYWORDS:
;  None.
;
; HISTORY:
;  2013-07-13: Written by Jared Males, jrmales@as.arizona.edu
;-
pro visao_getcal, NORTH_VISAO, visao_gx, visao_gy

NORTH_VISAO = -0.591

visao_gx = 159.8
visao_gy = -9.28

end
