;+
; NAME: visao_reg
;
; PURPOSE:
;  separate, register and align all of the raw science images and create an image cube
;
; INPUTS:
;  sci_imlist :
;  ref      :  the number of an image with a high-quality PSF to register against
;
; INPUT KEYWORDS:
;  unsharp  :  unsharp masks using gaussian of FWHM specified here and subtracts before aligning
;  clip     :  crops the array to specified width (square)
;  scl      :  rescale the continuum images by this factor
;
; OUTPUTS:
;
; OUTPUT KEYWORDS:
;    none
;
; EXAMPLE:
;
;
; HISTORY:
;  Written 2014-06-17 by Kate Follette, kfollette@as.arizona.edu
;  Modified April 2016 to include SDI images, write clip### into filename to avoid overwriting if want multiple FOV
;  Modified 4/27/16 to register against a single high quality science image and split out separating channels, dark and flat to separate procedure visao_separate_sdi
;  Modified April 2019 to pad arbitrarily, refine centers with circular symmetry, and interpolate over missing data (dust spots)
;-

pro visao_reg, ref, clip=clip, flat=flat, fwhm=fwhm, sdi=sdi, indiv=indiv, scl=scl, stp=stp, $
  fixpix=fixpix, refine_cen=refine_cen, mask=mask, rmax=rmax, pad=pad

  ;;naming specifications
  if keyword_set(flat) then namestr='_flat_' else namestr='_'
  if keyword_set(clip) then outstr = '_clip'+string(clip, format='(i03)')+string(namestr) else outstr=namestr

  ;;default values for rmax, pad, clip, and fwhm if not set in call
  if not keyword_set(rmax) then rmax=25
  if not keyword_set(pad) then pad=0
  if not keyword_set(FWHM) then fwhm=10.
  ;;if clip keyword not defined (not recommended), image size is full array (10254x512)
  if not keyword_set(clip) then clip=dim1

  ;; read in channel cubes from visao_separate_sdi, which are dark subtracted and flat fielded
  Line=readfits('Line'+string(namestr)+'preproc.fits', Linehead)
  Cont=readfits('Cont'+string(namestr)+'preproc.fits', Conthead)

  dim1=(size(Line))[1]
  dim2=(size(Line))[2]
  nims=(size(Line))[3]

  ;;pad the preprocessed images for later clipping
  ;;padding with NaNs. This makes the registration process a little slower but should avoid edge issues.
  Line_pad = dblarr(dim1+pad*2, dim2+pad*2, nims)*!Values.F_NAN
  Cont_pad = dblarr(dim1+pad*2, dim2+pad*2, nims)*!Values.F_NAN
  Line_pad[pad:dim1+pad-1,pad:dim2+pad-1,*]=Line
  Cont_pad[pad:dim1+pad-1,pad:dim2+pad-1,*]=Cont

  ;;release line and continuum cubes from memory
  delvar, Line
  delvar, Cont

  ;;redefine size and center geometry with pad
  dim1=dim1+pad*2
  dim2=dim2+pad*2
  xcen=(dim1-1)/2.
  ycen=(dim2-1)/2.

  ;;define geometry for clipping
  ;;shift into a single pixel for odd clip (+0.5,+0.5) from geometric center of 1024x512 array
  if clip mod 2 eq 0 then offset=0 else offset = -0.5
  xcen=xcen+offset
  ycen=ycen+offset
  trim = clip/2.-0.5

  ;;grab reference image to register against
  center_ref_line=Line_pad[*,*,ref-1]
  center_ref_cont=Cont_pad[*,*,ref-1]
  print, 'regsitering against image number', ref

  ;create a gaussian at reference pixel with FWHM~star FWHM (measure!) to register against
  ;;default FWHM is 10 unless set with keyword
  gauss_cen=psf_gaussian(npixel=[dim1, dim2], centroid=[xcen, ycen], FWHM=fwhm)

  ;smooth reference image to avoid registering on cosmic ray
  center_ref_line_smooth=gauss_smooth(center_ref_line, width=5.)*1.
  center_ref_cont_smooth=gauss_smooth(center_ref_cont, width=5.)*1.

  ;;move this reference image to center by cx correlating with gaussian
  method='F' ;for fourier cross-correlation
  subreg, gauss_cen, center_ref_line_smooth, sftline, method=string(method)
  subreg, gauss_cen, center_ref_cont_smooth, sftcont, method=string(method)

  ;;shift reference image to computed center
  center_ref_line=shift_interp(center_ref_line, sftline-1, spline=-0.5)
  center_ref_cont=shift_interp(center_ref_cont, sftcont-1, spline=-0.5)

  ;;create blank arrays for centering loop
  Line_smooth=dblarr(dim1, dim2, nims)
  Cont_smooth=dblarr(dim1, dim2, nims)
  Line_reg=dblarr(clip, clip, nims)
  Cont_reg=dblarr(clip,clip, nims)

  peak_ratio=dblarr(nims)
  expt=dblarr(nims)
  Line_shift_arr=dblarr(2,nims)
  Cont_shift_arr=dblarr(2,nims)
  crops=dblarr(4,nims)

  ;;image loop
  for i=0, nims-1 do begin
    ;smooth out cosmic rays for computing shifts
    Line_smooth[*,*,i]=gauss_smooth(Line_pad[*,*,i], width=5.)*1.
    Cont_smooth[*,*,i]=gauss_smooth(Cont_pad[*,*,i], width=5.)*1.

    ;; calculate offset from center of image by cx correlating with reference image
    subreg, center_ref_line, Line_pad[*,*,i], Line_shift, method=string(method)
    subreg, center_ref_cont, Cont_pad[*,*,i], Cont_shift, method=string(method)

    ;; record shifts for later use
    Line_shift_arr[*,i]=Line_shift
    Cont_shift_arr[*,i]=Cont_shift

    ;;shift images to center based on shifts computed with smoothed versions
    ;;overwrites original image in pad cube to save memory
    Line_pad[*,*,i]=shift_interp(Line_pad[*,*,i], Line_shift-1, spline=-0.5)
    Cont_pad[*,*,i]=shift_interp(Cont_pad[*,*,i], Cont_shift-1, spline=-0.5)

    ;;rescale continuum image by scl if keyword is set
    if keyword_set(scl) then begin
      Cont_pad[*,*,i]=rot(Cont_pad[*,*,i],0.,scl,cubic=-0.5)
    endif

    ;;trim centered images to specified clip
    Line_reg[*,*,i] = Line_pad[xcen-trim:xcen+trim,ycen-trim:ycen+trim,i]
    Cont_reg[*,*,i] = Cont_pad[xcen-trim:xcen+trim,ycen-trim:ycen+trim,i]
    print, xcen, ycen, trim

    print, 'processed image', i+1, '        of', nims

  endfor

  if keyword_set(fixpix) then begin
    print, "interpolating over NaN pixels using radial profile"
    for i=0, nims-1 do begin
      Lineim=Line_reg[*,*,i]
      Contim=Cont_reg[*,*,i]
      Line_rp = radprofim(Lineim)
      Cont_rp = radprofim(Contim)
      Line_idx = where(finite(Lineim) ne 1)
      Cont_idx = where(finite(Contim) ne 1)
      Lineim[Line_idx]=Line_rp[Line_idx]
      Contim[Cont_idx]=Cont_rp[Cont_idx]
      Line_reg[*,*,i]=Lineim
      Cont_reg[*,*,i]=Contim
    endfor
    ;;add interpolation method to headers
    sxaddpar, Linehead, 'FIXPIX', 'rpinterp'
    sxaddpar, Conthead, 'FIXPIX', 'rpinterp'
  endif

  ;;add to headers
  sxaddpar, Linehead, 'REFIM', ref
  sxaddpar, Conthead, 'REFIM', ref
  sxaddpar, Linehead, 'FWHM', fwhm
  sxaddpar, Conthead, 'FWHM', fwhm

  ;;write out registered cube
  writefits, 'Line'+string(outstr)+'reg.fits', Line_reg, Linehead
  writefits, 'Cont'+string(outstr)+'reg.fits', Cont_reg, Conthead

  ;;if keyword is specified, run center of circular symmetry algorithm to refine centers
  if keyword_set(refine_cen) then begin
    print, 'refining centers - this will take a few minutes'

    ;;reread in original images
    Line=readfits('Line'+string(namestr)+'preproc.fits', Linehead)
    Cont=readfits('Cont'+string(namestr)+'preproc.fits', Conthead)

    ;;original geometry
    dim1=(size(Line))[1]
    dim2=(size(Line))[2]
    nims=(size(Line))[3]

    ;;pad the original images again
    Line_pad = dblarr(dim1+pad*2, dim2+pad*2, nims)*!Values.F_NAN
    Cont_pad = dblarr(dim1+pad*2, dim2+pad*2, nims)*!Values.F_NAN
    Line_pad[pad:dim1+pad-1,pad:dim2+pad-1,*]=Line
    Cont_pad[pad:dim1+pad-1,pad:dim2+pad-1,*]=Cont

    ;;release line and continuum cubes from memory
    delvar, Line
    delvar, Cont

    ;;redefine geometry with pad
    dim1=dim1+pad*2
    dim2=dim2+pad*2
    xcen=(dim1-1)/2.
    ycen=(dim2-1)/2.

    ;;define geometry for clipping
    ;;shift into a single pixel for odd clip (-0.5,-0.5) from geometric center of 1024x512 array
    if clip mod 2 eq 0 then offset=0 else offset = -0.5
    xcen=xcen-offset
    ycen=ycen-offset
    trim = clip/2.-0.5

    ;;make medians of registered images from earlier loop
    Linemed=median(Line_reg, dim=3)
    Contmed=median(Cont_reg, dim=3)

    ;;make grid of possible centers to examine
    xr=indgen(101.)-100/2.
    yr=indgen(101.)-100/2.

    ;;run center of circular symmetry
    if keyword_set(mask) then begin ;for saturated images
      center_circlesym, Linemed, xr, yr, rmax, Line_xcen, Line_ycen, Line_grid, mask=mask
      center_circlesym, Contmed, xr, yr, rmax, Cont_xcen, Cont_ycen, Cont_grid, mask=mask
    endif else begin
      center_circlesym, Linemed, xr, yr, rmax, Line_xcen, Line_ycen, Line_grid
      center_circlesym, Contmed, xr, yr, rmax, Cont_xcen, Cont_ycen, Cont_grid
    endelse

    ;;compute shifts from returned centers
    Line_censhift=[(clip-1)/2.-Line_xcen,(clip-1)/2.-Line_ycen]
    Cont_censhift=[(clip-1)/2.-Cont_xcen,(clip-1)/2.-Cont_ycen]

    print, 'center of circular symmetry for median Line image is', Line_xcen, Line_ycen
    print, 'shifting Line images by', Line_censhift
    print, 'center of circular symmetry for median Continuum image is', Cont_xcen, Cont_ycen
    print, 'shifting Continuum images by', Cont_censhift

    for i=0, nims-1 do begin
      ;;overwrite full cube with shifted versions to avoid memory issues
      ;;shifts are originally computed shifts + center refinement
      Line_pad[*,*,i]=shift_interp(Line_pad[*,*,i], Line_shift_arr[*,i]-1+Line_censhift, spline=-0.5)
      Cont_pad[*,*,i]=shift_interp(Cont_pad[*,*,i], Cont_shift_arr[*,i]-1+Cont_censhift, spline=-0.5)
      Line_reg[*,*,i] = Line_pad[xcen-trim:xcen+trim,ycen-trim:ycen+trim,i]
      Cont_reg[*,*,i] = Cont_pad[xcen-trim:xcen+trim,ycen-trim:ycen+trim,i]
    endfor

    ;;add things to headers
    sxaddpar, Linehead, 'RMAX', rmax
    sxaddpar, Conthead, 'RMAX', rmax
    sxaddpar, Linehead, 'CENSHIFTX', Line_censhift[0]
    sxaddpar, Conthead, 'CENSHIFTX', Cont_censhift[0]
    sxaddpar, Linehead, 'CENSHIFTY', Line_censhift[1]
    sxaddpar, Conthead, 'CENSHIFTY', Cont_censhift[1]
    sxaddpar, Linehead, 'PAD', pad
    sxaddpar, Conthead, 'PAD', pad

    if keyword_set(fixpix) then begin
      print, "interpolating over NaN pixels using radial profile"
      for i=0, nims-1 do begin
        Lineim=Line_reg[*,*,i]
        Contim=Cont_reg[*,*,i]
        Line_rp = radprofim(Lineim)
        Cont_rp = radprofim(Contim)
        Line_idx = where(finite(Lineim) ne 1)
        Cont_idx = where(finite(Contim) ne 1)
        Lineim[Line_idx]=Line_rp[Line_idx]
        Contim[Cont_idx]=Cont_rp[Cont_idx]
        Line_reg[*,*,i]=Lineim
        Cont_reg[*,*,i]=Contim
      endfor
      ;;add interpolation method to headers
      sxaddpar, Linehead, 'FIXPIX', 'rpinterp'
      sxaddpar, Conthead, 'FIXPIX', 'rpinterp'
    endif

    ;;write files
    writefits, 'Line'+string(outstr)+'reg_refined.fits', Line_reg, Linehead
    writefits, 'Cont'+string(outstr)+'reg_refined.fits', Cont_reg, Conthead
  endif

  if keyword_set(stp) then  stop

  print, 'registration complete'

end
