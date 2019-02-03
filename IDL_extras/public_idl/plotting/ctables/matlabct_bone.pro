;+
; NAME: matlabct_cool
;
; PURPOSE:
;  Get the relative color values of the Matlab "bone" color table.
;
; DESCRIPTION:
;  Returns the relative color values of the Matlab "bone" color table.
;  Don't call directly, rather use matlabct.pro to actually interpolate and set the color table
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
;  matlab_bone, red_x, red_c, green_x, green_c, blue_x, blue_c
; 
;
; HISTORY:
;  Written 2014.09.13 by Jared Males, jrmales@email.arizona.edu
;
;-
pro matlabct_bone, red_x, red_c, green_x, green_c, blue_x, blue_c

red_c = [0., 4, 7, 11, 14, 18,21, 25,28, 32, 35, 39,  43,  46, 50, 53, 57,  60, 64, 67,  71,$
        74, 78, 81, 85, 89, 92, 96, 99, 103, 106, 110, 113,  117, 120, 124, 128, 131, 135, 138, 142, $
        145, 149, 152, 156, 159, 163, 166, 172, 178, 183, 189, 194, 200, 205, 211, 216, 222, 227, 233, $
        238, 244, 249, 255] /255. 
  
green_c = [0., 4,7,11 ,14,18,21,25,28 ,32,35,39,43,46,50,53,57,60,64,67,71,74,78,81,86,91,96,101,106,111,116,$
         120, 125,130,135,140,145,150,155,159,164,169,174,179,184,189,193,198,202,205,209,213,216,220,223,$ 
         227,230,234,237,241,244,248,251,255] /255.


blue_c =  [1.,6, 11, 16, 21, 26, 31, 35, 40, 45, 50, 55, 60, 65, 70, 74, 79, 84, 89, 94, 99,104,108,113,117,120,$
           124,128,131,135, 138,142,145,149  ,152  ,156  ,159,163,166,170,174,177,181,184,188,$
           191,195,198,202,205,209,213,216,220,223,227,230,234,237,241,244,248,251,255]/255.
           
red_x = findgen(n_elements(red_c)) / (1.*n_elements(red_c)  - 1.)
green_x = findgen(n_elements(green_c)) / (1.*n_elements(green_c)  - 1.)
blue_x = findgen(n_elements(blue_c)) / (1.*n_elements(blue_c)  - 1.)   

end