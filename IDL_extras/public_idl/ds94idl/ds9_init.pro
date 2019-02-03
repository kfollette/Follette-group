;+
; NAME: 
;   ds9_init
;
; PURPOSE:
;   Initialize the ds9 facility.
;
; DESCRIPTION:
;   Creates the underlying ds9 object, which will spawn a ds9 window if not already open.  Also initializes the 
;   ds9env common block, setting the imexam parameters to their defaults.
;
; INPUTS:
;   none
; 
; INPUT KEYWORDS:
;   force : if set, initializes even if ds9obj is not empty.  use with caution, as this will not 
;           deallocate an existing ds9obj.
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
pro ds9_init, force=force

common ds9env, $    ;this is the ds9 environment common block
       ds9obj    ;the ds9 object
       

if (n_elements(ds9obj) ne 1 or keyword_set(force)) then begin
   ds9obj = obj_new('ds9')
endif

end




