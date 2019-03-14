; subpix_centration_allangles.pro
; Katie Morzinski     2013 May 19     ktmorz@arizona.edu

; Centers images to subpixel accuracy.  (for saturated images)
; shifts by subpixel (tenth of pixel, from +/- 1 pixel around center
; For all angles:  Tests from 5 to 355 deg. rotational symmetry in 10 deg. increments
;                  and finds the average best centroid -- returns the image centered there.
; Just like subpix_centration.pro except does all angles.
; Set "satradius" to put NaN at the center radius you supply, for saturated pixels
;   It will place the circle at the center of the image so you have to keep iterating.
; If keyword set "perfect" it will iterate until the exact exact center is found at all angles

function subpix_centration_allangles, image, $
		satradius=satradius, $ ;radius of saturated data to ignore
		perfect=perfect, $     ;iterate until the exact center is found to fraction of pix specified perfect=0.1
		boxsize=boxsize, $     ;box size (diameter) within which to measure stddev (default 128 pix diam)
		debug=debug            ;print debug info messages

	targim = image
	target_image = image

	if keyword_set(perfect) then tol=perfect;tolerance -- size of offset allowed to be ~ zero

	n_angles = 36
	anglearr = (dindgen(n_angles)+0.5)*10.
	oxb_arr = dblarr(n_angles);best offset amount x
	oyb_arr = dblarr(n_angles);best offset amount y

	for k=0,n_angles-1 do begin
		angle = anglearr[k]
		ref_image = rot(targim,angle,cubic=-0.5)

		;subarrays
		if keyword_set(boxsize) then box=boxsize else box=128;diameter
		r1 = box/2.-box/4d0
		r2 = box/2.+box/4d0-1
		nx = (size(target_image))[1]
		ny = (size(target_image))[2]
		target = target_image[nx/2.-box/2.:nx/2.+box/2.-1,ny/2.-box/2.:ny/2.+box/2.-1]
		refim = ref_image[nx/2.-box/2.:nx/2.+box/2.-1,ny/2.-box/2.:ny/2.+box/2.-1]
		;; We carved out a subarray but we'll put it back into the full image at the end
		lcx = (nx-1)/2d0;large array center x
		lcy = (ny-1)/2d0;large array center y

		; Line up arrays to sub-pixel accuracy
		; subpix.pro
		numx = (size(refim))[1]
		numy = (size(refim))[2]
		nsh = 21 ;number of shifts

		tcx = (numx-1)/2d0;target center x
		tcy = (numy-1)/2d0;target center y
		rsh = dblarr(numx,numy);refim shifted
		this_rsh = dblarr(numx,numy,nsh,nsh);this refim shifted
		this_diff = dblarr(numx,numy,nsh,nsh)
		this_stddev = dblarr(nsh,nsh)
		for i=0,nsh-1 do begin;shift in x direction
			for j=0,nsh-1 do begin;shift in y direction
				this_refim = refim
				ox = (i-10)/10d0;offset x
				oy = (j-10)/10d0;offset y
				if keyword_set(satradius) then begin
					mask = circle(numx,numy, tcx-ox,tcy-oy, satradius, 1.0)
					mask = 1-mask ;invert
					mask[where(mask lt 1)] = float('NaN')
					;;Circle actually puts it at nx/2. instead of (nx-1)/2. like we want.
					;; and couldn't fix that -- so shift it here
					mask = rot(mask,0,1,tcx-ox-0.5,tcy-oy-0.5,cubic=-0.5,missing=1.0)
					this_refim = this_refim*mask
				endif
				this_rsh[*,*,i,j] = rot(this_refim,0d0,1d0,tcx-ox*1d0,tcy-oy*1d0,cubic=-0.5)
				this_diff[*,*,i,j] = this_rsh[*,*,i,j] - target
				this_stddev[i,j] = stddev(this_diff[r1:r2,r1:r2,i,j],/nan,/double)
			endfor;j
		endfor;i
		bs = array_indices([nsh,nsh],where(this_stddev eq min(this_stddev)),/dim);best shift
		oxb = ((bs[0]-10)/10.)/2d0;best offset amount x
		oyb = ((bs[1]-10)/10.)/2d0;best offset amount y
		oxb_arr[k] = oxb
		oyb_arr[k] = oyb
	endfor ;k

	; Plot it to see if centers changed drastically
	if keyword_set(debug) then begin
		;print,anglearr,lcx+oxb_arr,lcy+oyb_arr
		;imstat,lcx+oxb_arr,/v
		;imstat,lcy+oyb_arr
		plot,anglearr,lcx+oxb_arr,$
				xtitle='Angle of rotational symmetry tested (deg.)',ytitle='Centroid (pix.)',$
				yrange=[min([lcx+oxb_arr,lcy+oyb_arr])-.1,max([lcx+oxb_arr,lcy+oyb_arr])+.1],$
				xrange=[0,360],/xstyle,charsize=1.2,thick=3
		oplot,anglearr,lcy+oyb_arr,linestyle=2,thick=3
		al_legend,['x','y'],linestyle=[0,2],charsize=1,box=0,/bottom,/left
	endif

	;So now we have arrays of oxb and oyb, the best offsets in x and y,
	;for each angle of rotational symmetry tried, from 5 around to 355.
	; We simply take the mean x and y centroids as the best ones.
	; See illustration in subpix_centration_rotational_symmetry.tiff
	oxb = mean(oxb_arr,/double) ;Offset x Best
	oyb = mean(oyb_arr,/double) ;Offset y Best

	result = rot(target_image,0d0,1d0,lcx+oxb*1d0,lcy+oyb*1d0,cubic=-0.5)
	print,'Center of rotational symmetry: ',lcx+oxb,lcy+oyb



; 	if keyword_set(perfect) then begin
; 		target_image = result
; 		targim = target_image
; 		count=0
; 		while abs(oxb) gt tol or abs(oyb) gt tol do begin;if not perfect yet (zero offsets)

; 			for k=0,n_angles-1 do begin
; 				angle = anglearr[k]
; 				ref_image = rot(targim,angle*1d0,cubic=-0.5)

; 				;subarrays
; 				box = 128
; 				r1 = box/2.-box/4d0
; 				r2 = box/2.+box/4d0-1
; 				nx = (size(target_image))[1]
; 				ny = (size(target_image))[2]
; 				target = target_image[nx/2.-box/2.:nx/2.+box/2.-1,ny/2.-box/2.:ny/2.+box/2.-1]
; 				refim = ref_image[nx/2.-box/2.:nx/2.+box/2.-1,ny/2.-box/2.:ny/2.+box/2.-1]
; 				;; We carved out a subarray but we'll put it back into the full image at the end
; 				lcx = (nx-1)/2d0;large array center x
; 				lcy = (ny-1)/2d0;large array center y

; 				; Line up arrays to sub-pixel accuracy
; 				; subpix.pro
; 				numx = (size(refim))[1]
; 				numy = (size(refim))[2]
; 				nsh = 21 ;number of shifts

; 				tcx = (numx-1)/2d0;target center x
; 				tcy = (numy-1)/2d0;target center y
; 				rsh = fltarr(numx,numy);refim shifted
; 				this_rsh = fltarr(numx,numy,nsh,nsh);this refim shifted
; 				this_diff = fltarr(numx,numy,nsh,nsh)
; 				this_stddev = fltarr(nsh,nsh)
; 				for i=0,nsh-1 do begin;shift in x direction
; 					for j=0,nsh-1 do begin;shift in y direction
; 						this_refim = refim
; 						ox = (i-10)/10d0;offset x
; 						oy = (j-10)/10d0;offset y
; 						if keyword_set(satradius) then begin
; 							mask = circle(numx,numy, tcx-ox,tcy-oy, satradius, 1.0)
; 							mask = 1-mask ;invert
; 							mask[where(mask lt 1)] = float('NaN')
; 							;;Circle actually puts it at nx/2. instead of (nx-1)/2. like we want.
; 							;; and couldn't fix that -- so shift it here
; 							mask = rot(mask,0,1,tcx-ox-0.5,tcy-oy-0.5,cubic=-0.5,missing=1.0)
; 							this_refim = this_refim*mask
; 						endif
; 						this_rsh[*,*,i,j] = rot(this_refim,0,1,tcx-ox*1d0,tcy-oy*1d0,cubic=-0.5)
; 						this_diff[*,*,i,j] = this_rsh[*,*,i,j] - target
; 						this_stddev[i,j] = stddev(this_diff[r1:r2,r1:r2,i,j],/nan,/double)
; 					endfor;j
; 				endfor;i
; 				bs = array_indices([nsh,nsh],where(this_stddev eq min(this_stddev)),/dim);best shift
; 				oxb = ((bs[0]-10)/10.)/2d0;best offset amount x
; 				oyb = ((bs[1]-10)/10.)/2d0;best offset amount y
; 				oxb_arr[k] = oxb
; 				oyb_arr[k] = oyb
; 			endfor;k

; 			;So now we have arrays of oxb and oyb, the best offsets in x and y,
; 			;for each angle of rotational symmetry tried, from 5 around to 355.
; 			; We simply take the mean x and y centroids as the best ones.
; 			; See illustration in subpix_centration_rotational_symmetry.tiff
; 			oxb = mean(oxb_arr) ;Offset x Best
; 			oyb = mean(oyb_arr) ;Offset y Best

; 			result = rot(target_image,0d0,1d0,lcx+oxb*1d0,lcy+oyb*1d0,cubic=-0.5)
; 			;print,'Perfect while loop, Center of rotational symmetry: ',lcx+oxb,lcy+oyb
			
; 			count=count+1
; 			print,count,' = #Additional attempts to get perfect (within tolerance: ',tol,'pix)'

; 			target_image = result
; 			targim = target_image

; 			if abs(oxb) le tol and abs(oyb) le tol then $
; 				print,'Exiting Perfect while loop, Center of rotational symmetry: ',lcx+oxb,lcy+oyb
; 		endwhile ;while not perfect yet
; 	endif ;perfect



; 	return,result
; end
