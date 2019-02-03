; subpix_centration_allangles.pro
; Katie Morzinski     2013 May 19     ktmorz@arizona.edu

; Centers images to subpixel accuracy.  (for saturated images)
; shifts by subpixel (tenth of pixel, from +/- 1 pixel around center
; For all angles:  Tests from 5 to 355 deg. rotational symmetry in 10 deg. increments
;                  and finds the average best centroid -- returns the image centered there.
; Just like subpix_centration.pro except does all angles.

function subpix_centration_allangles, image, oxb, oyb, box=box, debug=debug, n_angles=n_angles

   ;r = rarr(128,128,/pix)
   ;std_idx = where(r gt 12 and r lt 64)
   
   if(n_elements(n_angles) ne 1) then n_angles = 36
	targim = image
	target_image = image

;	n_angles = 36
	anglearr = (findgen(n_angles)+0.5)*10.
	oxb_arr = fltarr(n_angles);best offset amount x
	oyb_arr = fltarr(n_angles);best offset amount y

	if(n_elements(box) ne 1) then box = 128.
	
	for k=0,n_angles-1 do begin
		angle = anglearr[k]
		ref_image = rot(targim,angle,cubic=-0.5)

		;subarrays
		;box = 128
		r1 = box/2.-box/4.
		r2 = box/2.+box/4.-1
		nx = (size(target_image))[1]
		ny = (size(target_image))[2]
		target = target_image[nx/2.-box/2.:nx/2.+box/2.-1,ny/2.-box/2.:ny/2.+box/2.-1]
		refim = ref_image[nx/2.-box/2.:nx/2.+box/2.-1,ny/2.-box/2.:ny/2.+box/2.-1]
		;; We carved out a subarray but we'll put it back into the full image at the end
		lcx = (nx-1)/2.;large array center x
		lcy = (ny-1)/2.;large array center y

		; Line up arrays to sub-pixel accuracy
		; subpix.pro
		numx = (size(refim))[1]
		numy = (size(refim))[2]
		nsh = 21 ;number of shifts

		tcx = (numx-1)/2.;target center x
		tcy = (numy-1)/2.;target center y
		rsh = fltarr(numx,numy);refim shifted
		this_rsh = fltarr(numx,numy,nsh,nsh);this refim shifted
		this_diff = fltarr(numx,numy,nsh,nsh)
		this_stddev = fltarr(nsh,nsh)
		for i=0,nsh-1 do begin;shift in x direction
			for j=0,nsh-1 do begin;shift in y direction
				ox = (i-10)/10.;offset x
				oy = (j-10)/10.;offset y
				this_rsh[*,*,i,j] = rot(refim,0,1,tcx-ox,tcy-oy,cubic=-0.5)
				this_diff[*,*,i,j] = this_rsh[*,*,i,j] - target
				imtest = this_diff[r1:r2,r1:r2,i,j]
				this_stddev[i,j] = stddev(imtest,/nan)
			endfor;j
		endfor;i
		bs = array_indices([nsh,nsh],where(this_stddev eq min(this_stddev)),/dim);best shift
		oxb = ((bs[0]-10)/10.)/2.;best offset amount x
		oyb = ((bs[1]-10)/10.)/2.;best offset amount y
		oxb_arr[k] = oxb
		oyb_arr[k] = oyb
	endfor;k

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
	oxb = mean(oxb_arr) ;Offset x Best
	oyb = mean(oyb_arr) ;Offset y Best

	result = rot(target_image,0,1,lcx+oxb,lcy+oyb,cubic=-0.5)
	print,'Center of rotational symmetry: ',lcx+oxb,lcy+oyb

	return,result
end
