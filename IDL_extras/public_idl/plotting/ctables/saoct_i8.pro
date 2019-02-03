;+
; NAME: saoct_i8
;
; PURPOSE:
;  Get the relative color values of the SAO "i8" color table.
;
; DESCRIPTION:
;  Returns the relative color values of the SAO "i8" color table, taken from default.C in the ds9 source.
;  Don't call directly, rather use saoct.pro to actually interpolate and set the color table
;
; INPUTS:
;   NONE
;
; INPUT KEYWORDS:
;   NONE
;
; OUTPUTS:
;  red_x   :  the x values (0<= x <=1) of the red points
;  red_c   :  the color values (0<= c<= 1) of the red points
;  green_x   :  the x values (0<= x <=1) of the green points
;  green_c   :  the color values (0<= c<= 1) of the green points
;  blue_x   :  the x values (0<= x <=1) of the blue points
;  blue_c   :  the color values (0<= c<= 1) of the blue points
;
; EXAMPLE:
;  
;  saoct_i8, red_x, red_c, green_x, green_c, blue_x, blue_c
; 
;
; HISTORY:
;  Written 2014.09.13 by Jared Males, jrmales@email.arizona.edu
;
;-
pro saoct_i8, red_x, red_c, green_x, green_c, blue_x, blue_c

red_c = [0., 0., 0.,0.,1.,1.,1.,1.]

green_c = [0., 1., 0., 1.,0., 1., 0., 1.]

blue_c = [0., 0., 1., 1., 0., 0., 1., 1.]

red_x = findgen(n_elements(red_c))/(1.*n_elements(red_c) - 1.)
green_x = red_x
blue_x = red_x
end

