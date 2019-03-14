;+
; NAME:
;	TV_GOHSS
;
; PURPOSE:
;	Display any image with various automatic gray levels.
;
;
; DESCRIPTION:
; 	Display any kind of 2-D image with automatic ZSCALE gray levels,
;	LOGARITHM scaling, histogram gaussfit ranging, stretch ranging, or
;	user defined range and zoom scale.
;
; CATEGORY:
;	Image display
;
; CALLING SEQUENCE:
;	TV_GOHSS, image [,wind=wind, Nsigma=Nsigma, range=range, log=log,
;			scale=scale, title=title, help=help, fit=fit,
;			zscale=zscale, contrast=contrast, restore=restore]
;
; INPUTS:
;	Image:	The image to be displayed: a 2-D array.

;
;
; KEYWORDS:
;
;	WIND:	Window number for the visualization. If negative, the image will be
;		visualized into the current window, if one is open.
;
;	TITLE: 	Title for the window.
;
;	SCALE:	Zoom factor for the image: e.g. 2.0 to half the image size.
;
;	POS: 	2-element vector specifying the position of the image
;		within the window (i.e. the device coordinate of
;		the bottom-left corner of the image). Units of pixels.
;
;	RANGE:	2-element vector [min_level, max_level] with the low cut and high cut
;		of the visualization.
;
;	ZSCALE:	If set, uses the Zscale algorithm (from IRAF) to autoscale with the
;		best visualization range.
;
;	CONTRAST: Contrast to use with the Zscale algorithm (default = 1)
;
;	FIT:	Uses the same range of the image: [min(image), max(image)]
;
;	NSIGMA:	If not zero, the levels will range for +/- NSIGMA times the sigma of
;		a gaussian fit of the histogram around the median.
;		The default, if no keywords are specified, is to use NSIGMA=3.
;
;	LOG:	Logarithm display.
;
;	HELP:	Prints this documentation.
;
;
; CALLS:
;
;	ZSCALE_RANGE is called if ZSCALE is set.
;	STAT is called if NSIGMA is set or no keywords are specified.
;
; MODIFICATION HISTORY:
;
;	G. Li Causi - Rome Astronomical Observatory, Dec 2003.
;
;		Planned improvements:
;			- Restore functionality of INTERACTIVE and RESTORE keywords.
;
;-

PRO tv_gohss, immagine, wind=wind, Nsigma=Nsigma, range=range, $
			LOG=LOG, scale=scale, title=title, $
			help=help, $
			fit=fit, zscale = zscale, contrast=contrast, $
			pos=pos


;*************
;Help Message:
;*************
IF KEYWORD_SET(help) THEN BEGIN
	DOC_LIBRARY, 'tv_gohss'
	GOTO, fine
ENDIF

;*************
;Not defined image:
;*************
IF n_elements(immagine) EQ 0 THEN RETURN

;*************
;Not 2-Dim image:
;*************
s = size(immagine)
IF s[0] NE 2 THEN RETURN

;***************************
;Not 2-element vector range:
;***************************
IF n_elements(range) NE 0 THEN IF n_elements(range) NE 2 THEN RETURN


;*************
;Defaults:
;*************

;Position:
IF n_elements(pos) NE 2 THEN pos = [0, 0]

;Contrast:
IF NOT KEYWORD_SET(contrast) THEN contrast = 1.

;Display scale:
IF NOT keyword_set(scale) THEN scale = 1. ELSE IF scale LE 0 THEN scale = 1.

;Window Number:
free_win=0
IF n_elements(wind) EQ 0 THEN BEGIN
	wind = 0
	free_win = 1
ENDIF
IF wind LT 0 AND !d.window LT 0 THEN BEGIN	;se wind e' negativo o se non ci sono finestre attive
	wind = 0
	free_win = 1
ENDIF

;Window title:
IF n_elements(title) EQ 0 THEN title = ""

;;Interactive:
;IF KEYWORD_SET(interactive) AND NOT KEYWORD_SET(Nsigma) THEN Nsigma=3.

;Nsigma:
IF n_elements(range) EQ 0 AND NOT KEYWORD_SET(fit) AND NOT KEYWORD_SET(Nsigma) AND NOT KEYWORD_SET(zscale) THEN Nsigma=3.


;************
;Make Window:
;************

;*************
;Image size:
;*************
xsize = s[1]
ysize = s[2]

;*************
;Display size:
;*************
xsize2 = xsize / scale
ysize2 = ysize / scale

;*************
;Window size:
;*************
min_wind_xsize = 100
xwind = xsize2 > min_wind_xsize
ywind = ysize2

IF wind GT 31 THEN free_win=1
IF wind GE 0 THEN window, wind, free=free_win, xsize=xwind, ysize=ywind, title=title	;se wind = -1 disegna sulla finestra corrente



;**************
;DISPLAY IMAGE:
;**************

;*****************************
;Not Finite elements in image:
;*****************************
finite_ind = where(FINITE(immagine) NE 0, count)
IF count EQ 0 THEN Message, 'Image has no finite elements'
finite_img = immagine[finite_ind]


;*****************************************
;Constant image or zero range (mean gray):
;*****************************************
zero_range = 0
IF n_elements(range) NE 0 THEN IF range[0] EQ range [1] THEN zero_range = 1
IF min(finite_img) EQ max(finite_img) OR zero_range EQ 1 THEN BEGIN
	;tv, immagine*0.+ (min(immagine)-min_sig)/(max_sig-min_sig)*255>0<255
	tv, immagine*0. + 127, pos[0], pos[1]
	range = [0, 255]
	print, 'WARNING: Constant Image (Max and Min are equals).'
	RETURN
ENDIF

;****************
;Z-Scale Display:
;****************
IF KEYWORD_SET(zscale) THEN BEGIN
	range = Zscale_Range(finite_img, contrast)
	Nsigma = 0.
ENDIF


;**************
;Display range:
;**************
IF KEYWORD_SET(Nsigma) THEN BEGIN
	stat = Stat(finite_img, /silent)
	min_sig = (stat.med - Nsigma * stat.gauss_sig)
	max_sig = (stat.med + Nsigma * stat.gauss_sig)
ENDIF
IF KEYWORD_SET(fit) THEN BEGIN
	min_sig = min(finite_img)
	max_sig = max(finite_img)
ENDIF
IF n_elements(range) NE 0 THEN BEGIN
	min_sig = range[0]
	max_sig = range[1]
ENDIF
IF KEYWORD_SET(restore) THEN BEGIN
	restore, "TV_GOHSS.sav"
	min_sig = mmin
	max_sig = mmax
ENDIF


;*********
;TV image:
;*********

i2 = congrid(immagine, xsize2, ysize2, /center)

;Logarithm display:
IF KEYWORD_SET(LOG) THEN tvscl, alog((i2 < max_sig > min_sig) - min_sig + 1), pos[0], pos[1]

;Linear display:
IF NOT KEYWORD_SET(LOG) THEN tvscl, i2 < max_sig > min_sig, pos[0], pos[1]




;********************
;INTERACTIVE DISPLAY:
;********************

IF KEYWORD_SET(INTERACTIVE) THEN BEGIN

	init = 1
	zoom_fac = 1.
	xorig = 0.
	yorig = 0.
	mmax = max_sig
	mmin = min_sig
	x0 = 0.5
	y0 = 0.5

	tvcrs, x0,y0, /normal	;plot the cursor at the center

	REPEAT BEGIN

		cursor, x,y, 2, /normal					;read mouse position

		IF !MOUSE.BUTTON EQ 1 THEN BEGIN
			x = x0
			y = y0
			tvcrs, x0,y0, /normal	;plot the cursor at the center
			init = 1
		ENDIF

		WHILE !MOUSE.BUTTON EQ 1 DO BEGIN

			cursor, x,y, 2, /normal				;get coordinates (wait for change)

			contr = x0 / x /2.					;da 1 a 0 per Nsigma (orizzontale)
			bright = y0 / y /2.					;da 0 a 1 per punto medio (verticale)

			IF init EQ 1 THEN BEGIN
				contr = 0.5
				bright = 0.5
				init = 0
			ENDIF

			px = [0.,.5,1.]
			py = [.01,1,10.]						;volte per Nsigma a 0,.5 e 1
			coef = exp_3_pts(px,py)
			nNsigma = coef[0]*exp(coef[1]*contr)+coef[2]

			px = [0.,.5,1.]
			py = [-100.,0,100.]						;volte per Nsigma a 0,.5 e 1
			nmed = (bright - px[1])*py[2]/(px[2]-px[1])

			mmax = (stat.med + nmed*Nsigma + nNsigma*Nsigma* stat.gauss_sig) > (stat.min + 1)
			mmin = (stat.med + nmed*Nsigma - nNsigma*Nsigma * stat.gauss_sig) < (stat.max -1)

			IF KEYWORD_SET(LOG) THEN BEGIN
				tvscl, alog((i2 < mmax > mmin) - mmin + 1)
			ENDIF ELSE BEGIN
				tvscl, i2 < mmax > mmin
			ENDELSE

		ENDWHILE

		key = get_kbrd(0)
		IF key EQ "" THEN GOTO, n

		IF key EQ "s" THEN save, mmin, mmax, filename="TV_GOHSS.sav"

		IF key EQ "+" OR key EQ "-" THEN BEGIN
			xc = x * xsize2 * scale / zoom_fac + xorig  	;in pixel originali
			yc = y * ysize2 * scale / zoom_fac + yorig
			IF key EQ "+" THEN zoom_fac = zoom_fac * 2.
			IF key EQ "-" THEN zoom_fac = zoom_fac / 2. > 1
			lx = xsize2 / zoom_fac * scale
			ly = ysize2 / zoom_fac * scale
			x1 = (xc - lx/2) > 0
			y1 = (yc - ly/2) > 0
			x2 = (x1 + lx - 1) < (xsize - 1)
			y2 = (y1 + ly - 1) < (ysize - 1)
			x1 = (x2 - lx + 1) > 0
			y1 = (y2 - ly + 1) > 0
			x2 = (x1 + lx - 1) < (xsize - 1)
			y2 = (y1 + ly - 1) < (ysize - 1)
			i2 = congrid(immagine[x1:x2, y1:y2], xsize2, ysize2, /center)
			xorig = x1
			yorig = y1
			window, !d.window, xsize=xwind, ysize=ywind, title=title + string(fix(zoom_fac)) + "x"
			IF KEYWORD_SET(LOG) THEN BEGIN
				tvscl, alog((i2 < mmax > mmin) - mmin + 1)
			ENDIF ELSE BEGIN
				tvscl, i2 < mmax > mmin
			ENDELSE
			tvcrs, x0,y0, /normal
		ENDIF
n:
	ENDREP UNTIL key EQ "q"

	max_sig = mmax
	min_sig = mmin

ENDIF



;******************
;Output used range:
;******************
IF n_elements(range) EQ 0 THEN range = [min_sig, max_sig]	;per far uscire i valori correnti




;;	IF n_elements(range) ne 0 THEN BEGIN
;;		;mmrange = float(max_sig - min_sig)
;;		;tv, ((i2 < max_sig > min_sig) - min_sig) / mmrange * 255.
;;		tvscl, i2 < max_sig > min_sig
;;	ENDIF ELSE BEGIN



fine:
END
