;+
; NAME: 
;   ds9_coord
;
; PURPOSE:
;   This procedure gets mouse click coordinates from the ds9 image display.
;
; DESCRIPTION:
;   Interacts with ds9, changing the cursor to notify you to click.  Returns once you click somewhere
;   in the displayed image, with the coordinates.  Works with data cubes, i.e. 3D images.  Note that ds9
;   gives coordinates starting with 1, but the ds94idl routines convert to 0 counting.
;
; INPUTS:
;   none
; 
; INPUT KEYWORDS:
;   image : the image to display, optional.
;
; OUTPUTS:
;   x  :  the x coordinate of the mouse click, converted to 0 counting.
;   y  :  the y coordinate of the mouse click, converted to 0 counting.
;   z  :  the z coordinate of the mouse click, converted to 0 counting.
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
pro ds9_coord, x, y, z, image=image

common ds9env, $    ;this is the ds9 environment common block
       ds9obj    ;the ds9 object
       

;check if ds9 is initialized
ds9_init

;display image if set
if(n_elements(image) gt 0) then ds9, image

;only ask for z coordinate if needed
if(arg_present(z)) then begin
   ds9obj->coordinate, tx, ty, tz
   z=tz
endif else begin
   ds9obj->coordinate, tx, ty
endelse

x=tx
y=ty


end

