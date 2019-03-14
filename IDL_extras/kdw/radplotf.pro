pro radplotf, image, x, y, out, fwhm, radpts, $
              inrad=inrad, outrad=outrad, drad=drad, $
              insky=insky, outsky=outsky, sky=sky, $
              ellratio=ellratio, ellpa=ellpa, $
              ang1=ang1, ang2=ang2, $
              badpix=badpix, fixpix = fixpix, $
              zpt = zpt, pixscl = pixscl, $
              plot=plot,  $
              outfile = outfile, $
              iterstat = iterstat, $
              silent=silent, verbose=verbose, $
			  nosky=nosky, median=median

;+
; PURPOSE:
; Program to calculate radial/elliptical profile of an image
; given aperture location, range of sizes, and inner and 
; outer radius for sky subtraction annulus. 
; NOTES:
; Calculates sky by
; median or user can specify a value.  Can extract profile
; from only a restricted range of angles (i.e. a sector) as well.
;
; User can enter a bad pixel mask; if so, bad pixels are
; excluded (not interpolated); total flux and differential
; flux are scaled from the fraction of pixels used, though
; the number of pixels is not.
;
; INPUTS
;  	image	
;	x,y	center location 
;
; OUTPUTS
;	out	output array for radial profile (see below)
;	fwhm	FWHM as determined by spline fit to radial profile
;	radpts
;
; KEYWORD INPUTS {defaults}
; 	inrad	inner radius for photometry {0.5}
;	outrad	outer radius for photometry {20}
;	drad	radial increment for photometry {1.0}
;	insky	inner radius for sky annulus {outrad+drad}
;	outsky	outer radius for sky annulues {insky+20}
;	sky	user specified constant sky value 
;	ellratio  for elliptical phot, ratio of major to minor axes 
;	ellpa	for elliptical phot, position angle CCW from up
;	        (with elliptical phot, all the radii become
;		the semi major axis)
;	ang1	starting angle for profile, measured CCW from up
;	ang2 	ending angle for profile (going CCW)
;       zpt     zeropoint (DN/sec from a 0 mag star) - for output file
;       pixscl  pixel scale (units/pixel) - for output file
;               (used only if 'zpt' is set)
;	badpix	bad pixel mask (0=bad, !0=good)
;	/silent	say nothing
;	/verbose print complete listing of results
;	/nosky	 don't do sky subtraction at all.
;
; OUTPUT FORMAT
; out(i,0) = radius (measured from middle of radial bin)
; within each radius:
; 	out(i,1)  = total flux 
;	out(i,2)  = net flux
; 	out(i,3)  = number of pixels (integrated)
;	out(i,4)  = differential flux (in radial bin bounded by out(0))
;	out(i,5)  = number of pixels (in radial bin)
;	out(i,6)  = average value in radial bin (sky subtracted)
;	out(i,7)  = stddev of pixel values in bin
;	out(i,8)  = median sky level
;	out(i,9)  = number of pixel in sky annulus
;	out(i,10) = stddev in sky annulus
;	out(i,11) = stddev of avg pixel values in radial bin
;		    (quadrature sum of out(7) and out(10))
;   if the user passes a zeropoint & pixel scale:
;       out(i,12)= radius in pixel scale units
;       out(i,13)= integrated magnitude
;       out(i,14)= avg surface brightness in radial bin (mag/units^2)
;
; radpts is a 2 x n array with the pixel values and distances for every
;   point inside the largest aperture (unsorted) where
;   radpts(0,*) are the distances and radpts(1,*) are the pixel values
;
; NOTES
; - to disable sky subtraction, set the /nosky keyword
; - implicitly assumes all noise is sky noise and gaussian
; - sector photometry is done only for the object, all flux in 
;   the sky annulus is assumed to be useable
;
; USES
;	inarc, splinefwhm
;
; HISTORY
; Written by M. C. Liu (UCB) 08/13/94
; 09/25/95 (MCL): added badpixel mask capability
; 10/25/95 (MCL): added elliptical photometry ability and
;	          ability to extract profile from a sector
; 10/06/98 (MCL): radii with no good pixels now return 0 flux in out array
; 10/22/98 (MCL): added /fixpix
; 02/01/99 (MCL): *** placed under RCS, version 1.0 ***
; 02/01/99 (MCL): can write results to 'outfile' w/explanatory header 
; 02/22/99 (MCL): adjusted output format of radial profile plot
;                 added 'zpt' and 'pixscl' keywords (used in output file)
; 03/15/99 (MCL): now uses ITERSTAT.PRO to determine median value in
;                   sky annulus, instead of ordinary MEDIAN
; 06/11/01 (MCL): added /iterstat for computing std dev in annuli
; 05/07/02 (MDP): added /nosky for use with simulations w/ no sky.
; 2002-12-05 (MDP): Added checks to ignore NaN pixels. NEEDS TO BE TESTED
; 2006-08-03 MDP: Added /MEDIAN option to compute median in each annulus.
; 					This is stored into out[*,5], in place of # of pixels.
; 					WARNING - potentially very very slow.

; 
;
; UNRESOLVED ISSUES
; - what are the proper values for inrad and drad?
; - would be slightly nicer if output was in form out a structure, instead
;   of having to know which numbers correspond to what
; - no accomodation for pixels which are partially in aperture; preferably
;   would like to use a simple weighting scheme (like in IRAF) to 
;   handle this (treats pixel coords at pixel corners, may be
;   preferable to work with pixel centers, i.e. (0,0) ->
;   (0.5,0.5), when accomodating fractional pixels)
; - output array should store # of bad pixels in aperture & radial
;   bins
; -> does not check for BADVAL pixels, but it should!
;-


; set defaults
if n_elements(inrad) eq 0 then inrad = 0.5*sqrt(2)
if keyword_set(outrad) eq 0 then outrad = 20.
if keyword_set(drad) eq 0 then drad=1.
if keyword_set(insky) eq 0 then insky = outrad+drad
if keyword_set(outsky) eq 0 then outsky = insky+drad+20.
if keyword_set(sky) then begin
	insky = outrad
	outsky = insky + 1.0
endif
if (keyword_set(ang1) xor keyword_set(ang2)) then begin
	print,'must set both ANG1 and ANG2 or set neither!'	
	retall
endif


; user explanation
if n_params() lt 3 then begin
	print,'pro radplotf,image,x,y,out,fwhm,radpts,'
	print,'             [inrad='+strc(inrad)+'],[outrad='+strc(outrad)+'],[drad='+strc(drad)+'],'
	print,'             [insky='+strc(insky)+'],[outsky='+strc(outsky)+'],[sky=],'
	print,'             [ellratio=],[ellpa=]'
	print,'             [ang1=],[ang2=]'
        print,'             [zpt=],[pixscl=],'
	print,'             [badpix=],[fixpix]'
	print,'             [plot],[silent],[verbose]'
	return
endif


; minimum number of sky pixels
MINSKY = 20


; sanity check of the photometry radii
if (outrad gt insky) or (outrad gt outsky) or (insky gt outsky) $
  or (inrad le 0.0) or (outsky le 0.0) or (insky le 0.0) or $
  (inrad ge outrad) then begin 
	message,'** your radii are screwed up! **'
endif
if not(keyword_set(silent)) then $
  print, '  inskyrad = ', strc(insky), ',  outskyrad = ', strc(outsky) 


; initialize arrays
inrad = float(inrad)
outrad = float(outrad)
drad = float(drad)
nrad = ceil((outrad-inrad)/drad) + 1
out = fltarr(nrad,12)


; extract relevant image subset (may be rectangular), translate coord origin,
;   bounded by edges of image
;   (there must be a cute IDL way to do this neater)
sz = size(image)
x0 = floor(x-outsky) 
x1 = ceil(x+outsky)   ; one pixel too many?
y0 = floor(y-outsky) 
y1 = ceil(y+outsky)
if (x0 lt 0.0) or (x1 gt (sz(1)-1)) or (y0 lt 0.0) or (y0 gt (sz(2)-1)) then $
  message,'sky apertures off image!', /info
x0 = x0 > 0.0
x1 = x1 < (sz(1)-1)
y0 = y0 > 0.0
y1 = y1 < (sz(2)-1)
nx = x1 - x0 + 1
ny = y1 - y0 + 1


; check phot apertures are on image
x0p = floor(x-outrad) 
x1p = ceil(x+outrad)   ; one pixel too many?
y0p = floor(y-outrad) 
y1p = ceil(y+outrad)
if (x0p lt 0.0) or (x1p gt (sz(1)-1)) or  $
  (y0p lt 0.0) or (y0p gt (sz(2)-1)) then $
  message,'phot apertures off image!!!', /info


; trim the image, translate coords
;print,'image trim: (',strc(fix(x0)),':',strc(fix(x1)),',',$
;		      strc(fix(y0)),':',strc(fix(y1)),')'
img = image(x0:x1,y0:y1)
xcen = x - x0
ycen = y - y0
if keyword_set(badpix) then bp = badpix(x0:x1,y0:y1)


; for debugging, can make some masks showing different regions
skyimg = fltarr(nx,ny)			; don't use /nozero!!
photimg = fltarr(nx,ny)			; don't use /nozero!!


; fix bad pixels if desired in image sub-region
; but preserve bad pixel mask info for output array info
if keyword_set(fixpix) then begin
    fixpix, img, bp ne 0, img, /silent
    ; according, set the bad pixel mask to all good
    ;bp = (bp-bp) + 1B
endif


; makes an array of (distance)^2 from center of aperture
;   where distance is the radial or the semi-major axis distance.
;   based on DIST_CIRCLE and DIST_ELLIPSE in Goddard IDL package,
;   but deals with rectangular image sections
distsq = fltarr(nx,ny,/nozero)
if keyword_set(ellratio) or keyword_set(ellpa) then begin

	if not(keyword_set(ellratio)) then ellratio=1.0
	if not(keyword_set(ellpa)) then ellpa = 0.0
	if not(keyword_set(silent)) then $
	  message,'elliptical phot, ratio='+strc(ellratio)+$
		', pa='+strc(ellpa),/info

	ang = ellpa /!RADEG                      
	cosang = cos(ang)
	sinang = sin(ang)
	xx = findgen(nx) - xcen
	yy = findgen(ny) - ycen

;	rotate pixels to match ellipse orientation
	xcosang = xx*cosang
	xsinang = xx*sinang

	for i = 0L,ny-1 do begin
	  xtemp =  xcosang + yy(i)*sinang
	  ytemp = -xsinang + yy(i)*cosang
	  distsq(0,i) = (xtemp*ellratio)^2 + ytemp^2 
	endfor

endif else begin

	xx = findgen(nx)
	yy = findgen(ny)
	x2 = (xx - xcen)^(2.0)
	y2 = (yy - ycen)^(2.0)
	for i = 0L,(ny-1) do $			; row loop
		distsq(*,i) = x2 + y2(i)

endelse


; mask out bad pixels by setting their distances too large for use
if keyword_set(badpix) and not(keyword_set(fixpix)) then begin
	w = where(bp eq 0,nbad)
	if (nbad gt 0) then distsq(w) = outsky^2.0+1
	if not(keyword_set(silent)) then $
          message, 'nbad= '+strc(nbad), /info
endif


;----------------------------------------------------------------------
; get sky level by masking and then medianing remaining pixels
; note use of "gt" to avoid picking same pixels as flux aperture
;  (checked graphically and there is no overlap between the 2 regions)
;  can disable sky subtraction then by setting insky=outsky
;----------------------------------------------------------------------
ns = 0
msky = 0.0
errsky = 0.0

if not(keyword_set(nosky)) then begin
if keyword_set(sky) then begin

	if not(keyword_set(silent)) then $
	  message,'using constant sky value',/info
	ns = 0
	msky = sky
	errsky = 0.0
	errsky2 = 0.0
	skyann = 0.0

endif else begin

	if not(keyword_set(silent)) then $
	  message,'finding sky level',/info
	in2 = insky^(2.0)
	out2 = outsky^(2.0)
	w = where((distsq gt in2) and (distsq le out2) and finite(img),ns)

	if ns ge 3 then begin
		; need at least 3 pixels to compute stddev of sky
		skyann = img(w)
                iterstat, skyann, istat, /silent
                msky = istat(2)
		;msky = median(skyann)
		errsky = stddev(skyann)
		skyimg(w) = -5.0
		photimg = skyimg
		if ns lt MINSKY then $
		  message,'* only '+strc(ns)+' sky pixels! *',/info
	endif else begin
		message,'* no pixels in sky annulus! *',/info
		ns = 0
		msky = 0.0
		errsky = 1.0e10
		skyann=0.0
	endelse


endelse
endif ; for if not(keyword_set(nosky))

errsky2 = errsky * errsky
out(*,8) = msky
out(*,9) = ns
out(*,10)= errsky
if keyword_set(verbose) then begin
	print,format='("sky stats: ",5(A," "))',strc(msky),strc(ns),strc(errsky),strc(min(skyann)),strc(max(skyann))
	print,'radii: ',inrad,outrad,drad,nrad
	print
endif


;----------------------------------------------------------------------
; now loop through photometry radii, finding the total flux, differential
;	flux, and differential average pixel value along with 1 sigma scatter
; 	relies on the fact the output array is full of zeroes
;----------------------------------------------------------------------
for i = 0,nrad-1 do begin

	dr = drad
	if i eq 0 then begin
		rin =  0.0
		rout = inrad
		rin2 = -0.01
	endif else begin
		rin = inrad + drad *(i-1)	
		rout = (rin + drad) < outrad
		rin2 = rin*rin
	endelse
	rout2 = rout*rout

; 	get flux and pixel stats in annulus, wary of counting pixels twice
;	checking if necessary if there are pixels in the sector
	w = where(distsq gt rin2 and distsq le rout2 and finite(img),np)

	pfrac = 1.0	; fraction of pixels in each annulus used

	if keyword_set(ang1) then begin
		whereis,img,w,xann,yann
		wang=where(inarc(xann,yann,ang1,ang2,xcen,ycen) ne 0,nang)
;		print,'angles removed ',strc(np-nang),' points'
		if (nang gt 0) then begin
		  w = w(wang) 
		  np = nang
		endif $
		  else np = 0
		pfrac = nang/np
	endif

	if np gt 0 then begin
		ann = img(w)
		dflux = total(ann,/NaN) * 1./pfrac
		if keyword_set(median) then dmed = median(ann)
		dnpix = np
		dnet = dflux - (dnpix * msky) * 1./pfrac
        davg = dnet / (dnpix * 1./pfrac)
                if np gt 1 then begin
                    if keyword_set(iterstat) then begin
                        iterstat, ann, istat, /silent
                        dsig = istat(3)
                    endif else $
                      dsig = stddev(ann) 
                endif else  $
                  dsig = 0.00

;		std dev in each annulus including sky sub error
		derr = sqrt(dsig*dsig + errsky2)
;		derr = sqrt(dsig*dsig + errsky2/dnpix)

		photimg(w) = rout2
	
		out(i,0) = (rout+rin)/2.0
		out(i,1) = out(i-1>0,1) + dflux
		out(i,2) = out(i-1>0,2) + dnet
		out(i,3) = out(i-1>0,3) + dnpix
		out(i,4) = dflux
		out(i,5) = dnpix
		if keyword_set(median) then out[i,5] = dmed
		out(i,6) = davg
		out(i,7) = dsig
		out(i,11) = derr
	endif else if (i ne 0) then begin
		out(i,0)= rout
		out(i,1:3) = out(i-1,1:3)
                out(i, 4:7) = 0.0
		out(i,11) = 0.0
	endif else begin
                out(i, 0) = rout
        endelse

;	print it out	
	if keyword_set(verbose) then begin
	  print,format='($,3(A," "))',strc(i),strc(rin),strc(rout)
	  print,format='($,3(A," "))',strc(out(i,1)),strc(out(i,2)),strc(out(i,3))
	  print,format='($,3(A," "))',strc(out(i,4)),strc(out(i,5)),strc(out(i,6))
	  print,format='($,2(A," "))',strc(out(i,7)),strc(out(i,11))
	  print
	endif

endfor


; fill radpts array after done with differential photometry
w = where(distsq ge 0.0 and distsq le outrad*outrad)
if keyword_set(ang1) then begin
	whereis,img,w,xann,yann
	wang=where(inarc(xann,yann,ang1,ang2,xcen,ycen) ne 0,nang)
	if (nang gt 0) then begin
	  w = w(wang) 
	  np = nang
	endif $
	  else np = 0
endif
radpts = dblarr(2,n_elements(w))
radpts(0,*) = sqrt(distsq(w))
radpts(1,*) = img(w)

if arg_present(fwhm) or keyword_set(plot) then begin 
	; compute FWHM via spline interpolation of radial profile
	fwhm = splinefwhm(out[*,0],out[*,6],/nopeakcheck,/fixpeak)
	if fwhm eq 999 then begin
		message,/info,"Unable to find FWHM!"
		;stop
	endif
endif

; plot
if keyword_set(plot) eq 1 then begin

	yy = !y.margin
	!y.margin = [4,4]

        title = '!Crad='+strc(inrad)+'-'+strc(outrad)+$
		', skyrad='+strc(round(insky))+'-'+strc(round(outsky))+$
		'!Cfwhm='+string(fwhm, '(F5.2)')+', sky='+strc(out(0,8))+$
		', center='+printcoo(x, y)
		;'!Cfwhm='+strc(fwhm)+', sky='+strc(out(0,8))+$
		;', center=('+strc(x)+','+strc(y)+')'

	if keyword_set(ellratio) then $
	  title = title + '!C ellratio='+strc(ellratio)+$
		  ', ellpa='+strc(ellpa)
	if keyword_set(ang1) then $
	  title = title + ', angles=['+strc(ang1)+','+strc(ang2)+']'
	
        if n_elements(radpts(1, *)) gt 100 then pp = 3 else pp = 1

        if alog10(max(radpts(1, *)) - min(radpts(1, *))) lt 4 then begin
            plot, radpts(0, *), radpts(1, *), psym = pp, xtitle = 'radius', $
              ytitle = 'raw counts', title = title
            oploterr, out(*, 0), out(*, 6)+out(*, 8), $
              out(*, 11)/sqrt(out(*, 5)), -4
            
;	    oplot,out(*,0),out(*,6)+out(*,8)
;	    oploterr,out(*,0),out(*,6),out(*,7)

        endif else begin
            plot_io, radpts(0, *), radpts(1, *), psym = pp, $
              xtitle = 'radius', ytitle = 'raw counts', title = title
            oplot, out(*, 0), out(*, 6)+out(*, 8), ps = -4
            oploterr, out(*, 0), out(*, 6)+out(*, 8), $
              out(*, 11)/sqrt(out(*, 5)), -4
        endelse

        ; overplot the phot radii
        plots, [outrad, outrad], !y.crange, line = 1, /clip
        plots, [insky, insky], !y.crange, line = 1, /clip
        plots, [outsky, outsky], !y.crange, line = 1, /clip
        
	!y.margin = yy
endif


; store results in output file, if desired
if keyword_set(outfile) then begin
    if filecheck(outfile) then begin
        openw, unit0, outfile, /get_lun
        printf, unit0, '# RADPLOTF.PRO: '+systime()
        printf, unit0, '# extracted radial profile'        
        printf, unit0, '# M. C. Liu (UCB)'        
        printf, unit0, '# '        
        printf, unit0, '# <INPUT PARMETERS>'        
        printf, unit0, '#   center = ', printcoo(x, y)      
        printf, unit0, '#   inner radius (pixels) = ', strc(inrad)
        printf, unit0, '#   outer radius (pixels) = ', strc(outrad)
        printf, unit0, '#   radius step size = ', strc(drad)
        if keyword_set(ellratio) then begin
            printf, unit0, '# '        
            printf, unit0, '#   elliptical aperture: axis ratio = ', ellratio
            printf, unit0, '#                        pos angle  = ', ellpa
        endif
        if keyword_set(ang1) then begin
            printf, unit0, '# '        
            printf, unit0, '#   using limited range of angle = ', printcoo(ang1, ang2)
        endif
        if keyword_set(fixpix) then $
          printf, unit0, '# ** FIXPIX used to fix bad pixels **'
        if keyword_set(sky) then  $
          printf, unit0, '# ** using user-entered constant sky level **'
        printf, unit0, '# '        
        printf, unit0, '# <OUTPUT FORMAT>'
        printf, unit0, '# out(i,0) = radius in pixels (each row)'
        printf, unit0, '#            (measured from middle of radial bin)'
        printf, unit0, '# and within each radius (row):'
        printf, unit0, '#   out(i,1) = total flux '
        printf, unit0, '#   out(i,2) = net flux'
        printf, unit0, '#   out(i,3) = number of pixels (integrated)'
        printf, unit0, '#   out(i,4) = differential flux in radial bin'
        printf, unit0, '#   out(i,5) = number of pixels in radial bin'
        printf, unit0, '#   out(i,6) = average value in radial bin (sky subtracted)'
        printf, unit0, '#   out(i,7) = stddev of pixel values in bin'
        printf, unit0, '#   out(i,8) = median sky level'
        printf, unit0, '#   out(i,9) = number of pixel in sky annulus'
        printf, unit0, '#   out(i,10)= stddev in sky annulus'
        printf, unit0, '#   out(i,11)= stddev of avg pixel values in radial bin'
        printf, unit0, '#              (quadrature sum of out(7) and out(10))'
        printf, unit0, '#'
        if keyword_set(zpt) then begin
            if not(keyword_set(pixscl)) then pixscl = 1.0
            printf, unit0, '# photometry inputs:'
            printf, unit0, '#   zeropoint (DN from 0 mag star) = ', strc(zpt), ' mag'
            printf, unit0, '#   pixel scale (units/pixel) = ', strc(pixscl)
            printf, unit0, '#'
            printf, unit0, '#   out(i,12)= radius in pixel scale units'
            printf, unit0, '#   out(i,13)= integrated magnitude'
            printf, unit0, '#   out(i,14)= avg surface brightness in radial bin (mag/units^2)'
            printf, unit0, '#'
            
            ; compute 3 extra columns with info in real units (mags, arcsec)
            outarr = fltarr(nrad, 15)-999.
            outarr(*, 0:11) = out
            outarr(*, 12) = outarr(*, 0)*pixscl
            wg = where(out(*, 2) gt 0.0, ng)
            if (ng gt 0) then  $
              outarr(wg, 13) = zpt - 2.5*alog10(out(wg, 2))
            wg = where(out(*, 4) gt 0.0, ng)
            if (ng gt 0) then  $
              outarr(wg, 14) = zpt - 2.5*alog10(out(wg, 6)/pixscl^2.0)
        endif else $
          outarr = out
        free_lun, unit0
        printarray, outarr, outfile = outfile, /silent
    endif
endif

end
