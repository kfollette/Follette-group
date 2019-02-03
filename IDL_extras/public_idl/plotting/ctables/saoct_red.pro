;+
; NAME: saoct_red
;
; PURPOSE:
;  Get the relative color values of the SAO "red" color table.
;
; DESCRIPTION:
;  Returns the relative color values of the SAO "red" color table, taken from default.C in the ds9 source.
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
;  saoct_red, red_x, red_c, green_x, green_c, blue_x, blue_c
; 
;
; HISTORY:
;  Written 2014.09.13 by Jared Males, jrmales@email.arizona.edu
;
;-
pro saoct_red, red_x, red_c, green_x, green_c, blue_x, blue_c

red_x = [0, .5, 1.0]
red_c = [0, 0.5, 1.]

green_x = [0, 0.5, 1.]
green_c = [0., 0., 0.]

blue_x = [0., 0.5, 1.]
blue_c = [0., 0., 0.]

end



