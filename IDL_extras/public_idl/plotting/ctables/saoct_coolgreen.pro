;+
; NAME: saoct_coolgreen
;
; PURPOSE:
;  Get the relative color values of the SAO "coolgreen" color table.
;
; DESCRIPTION:
;  Returns the relative color values of the SAO "coolgreen" color table.  This is not actually in the SAO toolkit
;  (or in ds9), but rather is a custom table based on "cool", where blue and green have been swapped.
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
;  saoct_coolgreen, red_x, red_c, green_x, green_c, blue_x, blue_c
; 
;
; HISTORY:
;  Written 2014.09.13 by Jared Males, jrmales@email.arizona.edu
;
;-
pro saoct_coolgreen, red_x, red_c, green_x, green_c, blue_x, blue_c

red_x = [0, .29, .76, 1]
red_c = [0, 0, .1, 1]

blue_x = [0, .22, .96, 1.0]
blue_c = [0., 0., 1., 1.]

green_x = [0., .53, 1.0]
green_c = [0., 1., 1.]

end

