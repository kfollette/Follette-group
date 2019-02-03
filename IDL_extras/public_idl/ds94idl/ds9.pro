;+
; NAME: 
;   ds9
;
; PURPOSE:
;   Display an image or image cube in ds9.
;
; DESCRIPTION:
;   Displays an image or image cube in ds9. Based on the ds9 class (ds9__define.pro) by Nick Konidaris, 
;   with modifications by Jared Males.  Requires that ds9 is installed, and that the XPA binaries are 
;   installed (xpaccess, xpaget, xpaset, etc.).
;
; INPUTS:
;   image :  the 2D image or 3D image cube to display, either an array of data or a filename
; 
; INPUT KEYWORDS:
;   preserve : if set, then ds9 will preserve the image scaling
;   frame    : set to a value >=1 to specify which frame to display the image in
;   newframe : if set, then a new frame will be added and the image displayed in it
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
pro ds9, image, preserve=preserve, newframe=newframe, frame=frame

common ds9env, $    ;this is the ds9 environment common block
       ds9obj   ;the ds9 object
       

;check if ds9 is initialized
ds9_init

if (size(image, /type) eq 7) then begin

   ds9obj->dispimage, mrdfits(image), preserve=keyword_set(preserve), newframe=keyword_set(newframe), frame=frame

endif else begin

   ds9obj->dispimage, image, preserve=keyword_set(preserve), newframe=keyword_set(newframe), frame=frame

endelse

end

