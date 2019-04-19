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
;
;-

pro visao_reg, ref, clip=clip, flat=flat, fwhm=fwhm, sdi=sdi, indiv=indiv, scl=scl, stp=stp, $
  fixpix=fixpix, refine_cen=refine_cen, mask=mask, rmax=rmax, pad=pad

  if keyword_set(flat) then namestr='_flat_' else namestr='_'
  if keyword_set(refine_cen) and not keyword_set(rmax) then begin
    print, "Are you sure you don't want to set rmax? Type .c if you really want to proceed"
    stop
  endif
  if not keyword_set(pad) then pad=0

  ;; read in channel cubes from visao_separate_sdi, which are dark subtracted and flat fielded
  Line=readfits('Line'+string(namestr)+'preproc.fits', Linehead)
  Cont=readfits('Cont'+string(namestr)+'preproc.fits', Conthead)

  dim1=(size(Line))[1]
  dim2=(size(Line))[2]
  nims=(size(Line))[3]

  ;;pad the preprocessed images for later clipping
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
  
  print, dim1, dim2, xcen, ycen

  ;;grab reference image to register against
  center_ref_line=Line_pad[*,*,ref-1]
  center_ref_cont=Cont_pad[*,*,ref-1]
  print, 'regsitering against image number', ref

  ;create a gaussian at very center of dummy array  with FWHM~star FWHM (measure!) to register against
  ;;default FWHM is 10 unless set with keyword
  if not keyword_set(FWHM) then fwhm=10.
  gauss_cen=psf_gaussian(npixel=[dim1, dim2], centroid=[xcen+offset, ycen+offset], FWHM=fwhm)

  ;smooth to avoid registering on cosmic ray
  center_ref_line_smooth=gauss_smooth(center_ref_line, width=5.)*1.
  center_ref_cont_smooth=gauss_smooth(center_ref_cont, width=5.)*1.

  ;;move this reference image to be approximately centered
  method='F' ;for fourier cross-correlation
  subreg, gauss_cen, center_ref_line_smooth, sft, method=string(method)
  subreg, gauss_cen, center_ref_cont_smooth, sft, method=string(method)
  
  center_ref_line=shift_interp(center_ref_line, sft-1+offset, spline=-0.5)
  center_ref_cont=shift_interp(center_ref_cont, sft-1+offset, spline=-0.5)

  ;;if clip keyword not defined (not recommended), image size is full array (10254x512)
  if not keyword_set(clip) then clip=dim1

  ;;create blank arrays

  Line_smooth=dblarr(dim1, dim2, nims)
  Cont_smooth=dblarr(dim1, dim2, nims)
  Line_reg=dblarr(clip, clip, nims)
  Cont_reg=dblarr(clip,clip, nims)
  SDI_im=dblarr(clip, clip, nims)

  peak_ratio=dblarr(nims)
  expt=dblarr(nims)
  Line_shift_arr=dblarr(2,nims)
  Cont_shift_arr=dblarr(2,nims)
  crops=dblarr(4,nims)

  ;;image loop
  for i=0, nims-1 do begin
    ;smooth out cosmic rays
    Line_smooth[*,*,i]=gauss_smooth(Line_pad[*,*,i], width=5.)*1.
    Cont_smooth[*,*,i]=gauss_smooth(Cont_pad[*,*,i], width=5.)*1.

    ;; calculate approximate offset from center of image
    subreg, center_ref_line, Line_pad[*,*,i], Line_shift, method=string(method)
    subreg, center_ref_cont, Cont_pad[*,*,i], Cont_shift, method=string(method)

    ;; record shifts
    Line_shift_arr[*,i]=Line_shift
    Cont_shift_arr[*,i]=Cont_shift

    ;;make pad values NaN 
    ;Line_pad[where(Line_pad eq -999)]='NaN'
    ;Cont_pad[where(Cont_pad eq -999)]='NaN'

    ;;shift images
    Line_pad[*,*,i]=shift_interp(Line_pad[*,*,i], Line_shift-1-offset*2, spline=-0.5)
    Cont_pad[*,*,i]=shift_interp(Cont_pad[*,*,i], Cont_shift-1-offset*2, spline=-0.5)

    ;;rescale continuum image by scl if keyword is set
    if keyword_set(scl) then begin
      Cont_pad[*,*,i]=rot(Cont_pad[*,*,i],0.,scl,cubic=-0.5)
    endif
    Line_reg[*,*,i] = Line_pad[xcen-trim:xcen+trim,ycen-trim:ycen+trim,i]
    Cont_reg[*,*,i] = Cont_pad[xcen-trim:xcen+trim,ycen-trim:ycen+trim,i]

    print, 'processed image', i+1, '        of', nims

  endfor

  ;;if keyword is specified, run center of circular symmetry algorithm to refine centers
  if keyword_set(refine_cen) then begin
    print, 'refining centers - this will take a few minutes'

    ;;reread in original images to avoid interpolating multiple times
    Line=readfits('Line'+string(namestr)+'preproc.fits', Linehead)
    Cont=readfits('Cont'+string(namestr)+'preproc.fits', Conthead)
    dim1=(size(Line))[1]
    dim2=(size(Line))[2]
    nims=(size(Line))[3]
    
    ;;pad the preprocessed images for later clipping
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

    ;;make medians of registered images
    Linemed=median(Line_reg, dim=3)
    Contmed=median(Cont_reg, dim=3)
    
    ;;make grid of possible centers to examine
    xr=indgen(51.)-50/2.
    yr=indgen(51.)-50/2.

    if keyword_set(rmax) then lim=rmax else lim=dim1/2.

    ;;compute center
    if keyword_set(mask) then begin
      center_circlesym, Linemed, xr, yr, lim, Line_xcen, Line_ycen, Line_grid, mask=mask
      center_circlesym, Contmed, xr, yr, lim, Cont_xcen, Cont_ycen, Cont_grid, mask=mask
    endif else begin
      center_circlesym, Linemed, xr, yr, lim, Line_xcen, Line_ycen, Line_grid
      center_circlesym, Contmed, xr, yr, lim, Cont_xcen, Cont_ycen, Cont_grid
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
      Line_pad[*,*,i]=shift_interp(Line_pad[*,*,i], Line_shift_arr[*,i]-1-offset*2+Line_censhift, spline=-0.5)
      Cont_pad[*,*,i]=shift_interp(Cont_pad[*,*,i], Cont_shift_arr[*,i]-1-offset*2+Cont_censhift, spline=-0.5)
      Line_reg[*,*,i] = Line_pad[xcen-trim:xcen+trim,ycen-trim:ycen+trim,i]
      Cont_reg[*,*,i] = Cont_pad[xcen-trim:xcen+trim,ycen-trim:ycen+trim,i]
    endfor
        
  endif
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
  endif

  sxaddpar, Linehead, 'REFIM', ref
  sxaddpar, Conthead, 'REFIM', ref
  if keyword_set(fixpix) then begin
    sxaddpar, Linehead, 'FIXPIX', 'rpinterp'
    sxaddpar, Conthead, 'FIXPIX', 'rpinterp'
  endif

  if keyword_set(refine_cen) then begin
    if keyword_set(clip) then begin
      writefits, 'Line_clip'+string(clip, format='(i03)')+string(namestr)+'reg_refined.fits', Line_reg, Linehead
      writefits, 'Cont_clip'+string(clip, format='(i03)')+string(namestr)+'reg_refined.fits', Cont_reg, Conthead
    endif else begin
      writefits, 'Line'+string(namestr)+'reg_refined.fits', Line_reg, Linehead
      writefits, 'Cont'+string(namestr)+'reg_refined.fits', Cont_reg, Conthead
    endelse
  endif else begin
    if keyword_set(clip) then begin
      writefits, 'Line_clip'+string(clip, format='(i03)')+string(namestr)+'reg.fits', Line_reg, Linehead
      writefits, 'Cont_clip'+string(clip, format='(i03)')+string(namestr)+'reg.fits', Cont_reg, Conthead
    endif else begin
      writefits, 'Line'+string(namestr)+'reg.fits', Line_reg, Linehead
      writefits, 'Cont'+string(namestr)+'reg.fits', Cont_reg, Conthead
    endelse
  endelse

  if keyword_set(stp) then  stop

  print, 'registration complete'

end
