;+
; NAME: 
;   ds9_shutdown
;
; PURPOSE:
;   This procedure deallocates the ds9 common block.
;
; DESCRIPTION:
;   Destroys the ds9 object, and sets the stored image array to 0.  
;
; INPUTS:
;   none
; 
; INPUT KEYWORDS:
;   none
;
; OUTPUTS:
;   none
;
; OUTPUT KEYWORDS:
;   none
;
; MODIFICATION HISTORY:
;  Written 2013/04/12 by Jared Males (jrmales@email.arizona.edu)
;
; BUGS/WISH LIST:
;
;-
pro ds9_shutdown

common ds9env, $    ;this is the ds9 environment common block
       ds9obj    ;the ds9 object
      

obj_destroy, ds9obj
undefine, ds9obj


end




