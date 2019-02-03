;+
; NAME: 
;   ds9_centroid
;
; PURPOSE:
;   Calculates centroid with given coords as initial guess.
;
; DESCRIPTION:
;   Using gcntrd, this procedure centroids in ds9image (in the ds9env common block).  The x, y,z inputs are 
;   the initial guesses and ds9obj->fwhm() is the full-width at half-maximum supplied to gcntrd.
;   If no output variables are passed, or the printout keyword is set, then results are printed.
;
; INPUTS:
;   xin  :  the input x coordinate.  If none of xin,yin,zin are supplied then ds9 cursor is used.
;   yin  :  the input y coordinate.  If none of xin,yin,zin are supplied then ds9 cursor is used.
;   zin  :  the input z coordinate.  If none of xin,yin,zin are supplied then ds9 cursor is used.
;
; OUTPUTS:
;  xout  :  the output x coordinate of the centroid.
;  yout  :  the output y coordinate of the centrid.
;  sep   :  the separation of the centroid position from the image center
;  pa    :  the position angle of the centroid, E of N, from the image center.
;
; OUTPUT KEYWORDS:
;  subim    : the subimage extracted for centroiding.
;  zout     : the out put z coordinate, in case it's a cube
;  printout : prints the output if set
;
; MODIFICATION HISTORY:
;  Written 2013/05/26 by Jared Males (jrmales@email.arizona.edu)
;
; BUGS/WISH LIST:
;
;-
pro ds9_centroid, xout, yout, sep, pa, xin, yin, zin, subim=subim, zout=zout, printout=printout


common ds9env, $    ;this is the ds9 environment common block
       ds9obj   ;the ds9 object

if(n_elements(xin) lt 1 or n_elements(yin) lt 1 or n_elements(zin) lt 1 ) then begin
   ds9_coord, xin, yin, zin
endif
       
if(xin eq -1 or yin eq -1 or zin eq -1) then begin
   ds9_coord, xin, yin, zin
endif

ds9_subim, xin, yin, zin, ds9obj->plotw(), subim

ds9obj->getimage, ds9im

dim1 = (size(ds9im))[1]
dim2 = (size(ds9im))[2]

gcntrd, subim, 0.5*ds9obj->plotw(), 0.5*ds9obj->plotw(), xcen, ycen, ds9obj->fwhm()
   
xout = xcen + floor(xin-0.5*ds9obj->plotw())
yout = ycen + floor(yin-0.5*ds9obj->plotw())

if (arg_present(zout)) then zout = zin

imxc = ds9obj->xcen()
if(imxc eq -1) then imxc = 0.5*(dim1-1.)
imyc = ds9obj->ycen()
if(imyc eq -1) then imyc = 0.5*(dim2-1.)

print, imxc, imyc
print, xout, yout

sep = sqrt((yout-imyc)^2 + (xout-imxc)^2)
pa =  (atan(yout-imyc, xout-imxc)/!dtor + 3600. - 90.) mod 360.

if(~arg_present(xout) or keyword_set(printout)) then begin

   print, 'Centroid:'
   print, format='(A8,A8,A8,A8)', 'x','y', 'sep', 'PA'
      
   print, xout, yout, sep, pa, format='(F8.2,F8.2,F8.2,F8.2)'
endif

end
