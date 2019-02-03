;+
; NAME: psplot_start
;
; PURPOSE:
;  Start postscript plotting with nice choices for font, etc.  To stop ps plotting call psplot_stop.
;
; DESCRIPTION:
;  Sets the device to 'ps'.  Also turns off vector fonts, and changes to bold courier font, which
;  typically looks good on a projector display.  Additional wrapping for colortable management is provided.
;
; INPUTS:
;  fname:   the name of the output file
;
; INPUT KEYWORDS:
;  colortable:  name of the colortable, to be passed to mxloadct.pro
;  ncolors:   number of colors to use
;  landscape:  set the landscape flag.  if not set the portrait flag is set
;  bits_per_pixel: bits per pixel to set on the device
;
; OUTPUTS:
;  none
;
; EXAMPLE:
;  psplot_start, 'goodplot.ps', colortable='idl:13', ncolors=256, /landscape
;
; HISTORY:
;  Written 2013 by Jared Males, jrmales@email.arizona.edu
;
;-
pro psplot_start, fname, colortable=colortable,  ncolors=ncolors,  landscape=landscape, bits_per_pixel=bits_per_pixel

;So that psplot_stop knows how to change the filename in the banner
common psplot, filename
filename = fname

if(n_elements(colortable) eq 0) then colortable='idl:0'

if(n_elements(ncolors) ne 1) then ncolors=256

if(n_elements(bits_per_pixel) ne 1) then bits_per_pixel=8

set_plot, 'ps'

device, filename=fname

;Turn off vector fonts, use postscript fonts
!p.font=0
device,/COURIER,/BOLD,ISOLATIN1=1 ;Jared's default font.

device, BITS_PER_PIXEL=bits_per_pixel, /COLOR


if(keyword_set(landscape)) then begin
   device, /landscape, /cmyk
endif else begin
   device, /portrait, /cmyk
endelse

;loadct,colortable,ncolors=ncolors
mxloadct, colortable, ncol=ncolors

end

