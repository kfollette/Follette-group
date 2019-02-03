;+
; NAME: center_cube
; 
; DESCRIPTION:
;   Centers a cube of images by either centroiding on a given star or using rotational centration.
;
; INPUTS:
;   ims  :  the cube of 2D input images
;   xc   : approximate x coord of the star to center, either a single value or a vector of length nims
;   yc   : approximate y coord of the star to center, either a single value or a vector of length nims  
;  
; INPUT KEYWORDS:
;   noshift   : if set, the images are not shifted
;   fwhm      : the fwhm parameter for gcntrd, default is 4
;   rotcen    : if set, rotational centration is used, otherwise gcntrd is used
;   precenter : if set and using rotational centration, a single 180 rotation is used first to get close.
;   box  :  the box size to use for rotational centration.  default is 32.
; OUTPUTS:
;   ims     : on completion, the input array contains the centered images, unless noshift is set.
;   xstar   : the shifts in x
;   ystar   : the shifts in y
;
; MODIFICATION HISTORY:
;  Written 2013/06/16 by Jared Males (jrmales@email.arizona.edu)
;
; BUGS/WISH LIST:
;  Should pass interp parameters via keywords
;
;-
pro center_cube, ims, xc, yc, xstar, ystar, noshift=noshift, fwhm=fwhm, $
                              rotcen=rotcen, precenter=precenter, noall=noall, box=box, shiftonly=shiftonly, padsz=padsz, mask=mask, n_angles=n_angles, padval=padval, padran=padran, good=good, silent=silent

;--------------- Setup ------------------------------

get_cubedims, ims, dim1, dim2, nims

if n_elements(padsz) eq 1 then begin

   if (n_elements(padval) ne 1) then padval = 0

   if(n_elements(padran) eq 1) then begin
   
      padims = randomn(s1, dim1+2*padsz, dim2+2*padsz, nims) * padran
      
   endif else begin
   
      padims = fltarr(dim1+2*padsz, dim2+2*padsz, nims) + padval
      
   endelse
   
   for i=0, nims-1 do padims[padsz:padsz+dim1-1,padsz:padsz+dim2-1,i] = ims[*,*,i]

   ims = padims
   
   padims = 0
   
   get_cubedims, ims, dim1, dim2, nims
   
endif else begin
   padsz = 0
endelse


x0=.5*(dim1-1)
y0=.5*(dim2-1)

if(n_elements(fwhm) ne 1) then fwhm = 4.

xstar = dblarr(nims)
ystar = dblarr(nims)

if keyword_set(shiftonly) then begin

   for i=0, nims-1 do begin
      ims[*,*,i] = rot(ims[*,*,i], 0., 1., xc[i] + padsz, yc[i] + padsz, cubic=-0.5)
   endfor
   
   return
endif
      

if keyword_set(rotcen) then begin

   
   if(n_elements(box) lt 1) then box = 32.
   
   for i=0, nims-1 do begin
   
      status = 'center_cube: rot. centering: ' + strcompress(string(i+1) + '/' + string(nims), /rem)
      statusline, status, 0
      
      im = ims[*,*,i]

      xstar[i] = 0.
      ystar[i] = 0.
      if (keyword_set(precenter)) then begin
         im = subpix_centration( im, 180., dx, dy, box=box)
         xstar[i] = xstar[i] + dx
         ystar[i] = ystar[i] + dy
      endif
   
      if(~keyword_set(noall)) then begin   
         im = subpix_centration_allangles(im, dx, dy, box=box, n_angles=n_angles)
         xstar[i] = xstar[i] + dx
         ystar[i] = ystar[i] + dy
      endif
   
      if (~keyword_set(noshift)) then ims[*,*,i] = im
      
   endfor
   
endif else begin ;use gcntrd

   good = intarr(nims)+1.
   txc = xc + padsz
   if n_elements(xc) eq 1 then txc = dblarr(nims)+xc + padsz

   tyc = yc + padsz
   if n_elements(yc) eq 1 then tyc = dblarr(nims)+yc + padsz
   
   for i=0, nims-1 do begin
   
      status = 'center_cube: centering ' + strcompress(string(i+1) + '/' + string(nims), /rem)
      statusline, status, 0
   
      gcntrd, ims[*,*,i], txc[i], tyc[i], xst, yst, fwhm
      
      if(xst eq -1 or yst eq -1) then begin
         good[i] = 0.
         if(~keyword_set(silent)) then print, 'bad center #', i
         xst = txc[i]
         yst = tyc[i]
      endif
      
      xstar[i] = xst
      ystar[i] = yst
      
      if(~keyword_set(noshift)) then begin
         ims[*,*,i] = rot(ims[*,*,i], 0., 1., xstar[i], ystar[i], cubic=-0.5)
      endif
   endfor
endelse
      

        
statusline, /clear

q = -1
if keyword_set(mask) then begin
   
   maskno = float('NaN')
   r = rarr(dim1, dim2, xarr, yarr, /pix)
   xarr = xarr - min(xarr)
   yarr = yarr - min(yarr)
   
   for i=0, nims-1 do begin
   
      x0 = dim1 - xstar[i]; + .5*(dim1-1.)
      y0 = dim2 - ystar[i]; + .5*(dim2-1.)

      idx = where( xarr lt x0-.5*dim1+padsz+1 or xarr gt x0+.5*dim1-padsz-1 or yarr lt y0-.5*dim2+padsz+1 or yarr gt y0+.5*dim2-padsz-1 )
   
      im = ims[*,*,i]
      im[idx] = maskno
      ims[*,*,i] = im
      
   endfor
endif

   
   
end

