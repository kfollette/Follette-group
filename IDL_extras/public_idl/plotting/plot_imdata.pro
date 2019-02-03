;+
; NAME: plot_imdata
;
; PURPOSE:
;  Make a 2D image plot, using the plotimage procedure, using platescale to size the plot. 
;
; DESCRIPTION:
;  Uses the provided platescale and desired plot dimensions to scale the image.  Provides axis labels
;  and tic marks, etc.
;
; INPUTS:
;  im         : the image data to plot   
;  platescale : the platescale of the image
;  plot_dim   : the desired size of the image, in the units of platescale
;  minc       : the minimum value for the colortable
;  maxc       : the maximum value for the colortable
;
; INPUT KEYWORDS:
;  x0        : the center of the plot, default is 0.5*(dim1-1), where dim1 is the number of columns
;  y0        : the center of the plot, default is 0.5*(dim2-1), where dim2 is the number of rows
;  ncolors   : the number of colors to use, default 256
;  ra        : reverse the x-axis labels to match right ascension
;  xtitle    : the title of the x-axis, default is "\Deltax (pixels)"
;  ytitle    : the title of the x-axis, default is "\Deltax (pixels)"
;  charsize  : the normal charsize keyword, default is 1.75
;  charthick : the normal charthick keyword, default is 3
;  xythick   : if set, then xthick and ythick are set to this value.
;  xthick    : the normal xthick keyword, default is 4. 
;  ythick    : the normal ythick keyword, default is 4.
;  bias      : the color table bias, same meaning as in ds9. Default = 0.5
;  contrast  : the color table contrast, same meaning as in ds9.  Default = 1.0
;  _extra    : other plot keywords are passed to plotimage unaltered
;
; OUTPUTS:
;  none
;
; EXAMPLE:
;  This plots the image im, with 0.00791 arcsec/pixel, with a width of 1.5 arcsec, and treats
;  the x-axis as right ascension
;
;  plot_imdata, im, 0.00791, 1.5, -5., 5., /ra, xtitle='\DeltaRA(arcsec)', xtitle='\DeltaDec(arcsec)'
; 
;
; HISTORY:
;  Written 2013 by Jared Males, jrmales@email.arizona.edu
;
;-
pro plot_imdata, im, platescale, plot_dim, minc, maxc,$
                  x0=x0, y0=y0,$
                  ncolors=ncolors, ra=ra, xtitle=xtitle, ytitle=ytitle, $
                  charsize=charsize, charthick=charthick, xythick=xythick, xthick=xthick, ythick=ythick, $
                  bias=bias, contrast=contrast,$
                  _extra=_extra
   
   
if(n_elements(xtitle) ne 1) then xtitle = textoidl('\Deltax (pixels)')
if(n_elements(ytitle) ne 1) then ytitle = textoidl('\Deltay (pixels)')

if(n_elements(charsize) ne 1) then charsize=1.75
if(n_elements(charthick) ne 1) then charthick=3

if(n_elements(ncolors) ne 1) then ncolors=256

if(n_elements(xythick) ne 1) then xythick=4

if(n_elements(xthick) ne 1) then xthick=xythick
if(n_elements(ythick) ne 1) then ythick=xythick


dim1 = double((size(im))[1])
dim2 = double((size(im))[2])


; ------ Stretch the color table --------
if(n_elements(bias) ne 1) then bias = 0.5
if(n_elements(contrast) ne 1) then contrast = 1.

;always call so that bottom color is black
rescale_ct, bias, contrast

; ------ Scale the Image -------------
RESULT = BYTSCL(im, MIN=minc, MAX=maxc, TOP=(ncolors-1))



xmin = -0.5d*(dim1)*platescale
xmax = xmin + ((dim1)*platescale)

ymin = -0.5d*(dim2)*platescale
ymax = ymin + ((dim2)*platescale)


if(n_elements(x0) ne 1) then x0 = 0.
if(n_elements(y0) ne 1) then y0 = 0.

plot_xmin = x0 - 0.5*plot_dim
plot_ymin = y0 - 0.5*plot_dim

if(keyword_set(ra)) then begin
   xmin = -1.*xmin
   xmax = -1.*xmax
   xrange=[plot_xmin+plot_dim, plot_xmin]
endif else begin
   xrange=[plot_xmin, plot_xmin + plot_dim]
endelse

print, xmin, xmax
print, ymin, ymax

plotimage, result, xrange=xrange,yrange=[plot_ymin,plot_ymin+plot_dim], $
imgxrange=[xmin, xmax],imgyrange=[ymin, ymax], /preserve_asp, $
ytitle=ytitle, xtitle=xtitle, xthick=xthick, ythick=ythick, charthick=charthick, charsize=charsize, $
xstyle=xstyle, ystyle=ystyle, _extra=_extra



end


