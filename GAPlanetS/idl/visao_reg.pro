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
  fixpix=fixpix, refine_cen=refine_cen, mask=mask, rmax=rmax

  if keyword_set(flat) then namestr='_flat_' else namestr='_'

  ;; read in channel cubes from visao_separate_sdi, which are dark subtracted and flat fielded
  if not keyword_set(indiv) then begin
    Line=readfits('Line'+string(namestr)+'preproc.fits', Linehead)
    Cont=readfits('Cont'+string(namestr)+'preproc.fits', Conthead)

    ;;grab reference image to register against
    center_ref_line=Line[*,*,ref-1]
    center_ref_cont=Cont[*,*,ref-1]
    print, 'regsitering against image number', ref
  endif else begin
    center_ref_line=readfits('./indiv/Line'+string(namestr)+string(ref, format='(i04)')+'.fits')
    center_ref_cont=readfits('./indiv/Cont'+string(namestr)+string(ref, format='(i04)')+'.fits')
    print, 'regsitering against', './indiv/Line/Cont'+string(namestr)+string(ref, format='(i04)')+'.fits'
  endelse

  dim1=(size(center_ref_line))[1]
  dim2=(size(center_ref_line))[2]
  if keyword_set(indiv) then begin
    spawn, 'ls indiv/Line* | wc -l', nims
    nims=long(nims[0])
  endif else begin
    nims=(size(Line))[3]
  endelse
  xcen=(dim1-1)/2.
  ycen=(dim2-1)/2.

  ;create a gaussian at very center of dummy array  with FWHM~star FWHM (measure!) to register against
  ;;default FWHM is 10 unless set with keyword
  if not keyword_set(FWHM) then fwhm=10.
  gauss_cen=psf_gaussian(npixel=[dim1, dim2], centroid=[xcen, ycen], FWHM=fwhm)

  ;smooth to avoid cosmic ray
  center_ref_line_smooth=gauss_smooth(center_ref_line, width=5.)*1.
  center_ref_cont_smooth=gauss_smooth(center_ref_cont, width=5.)*1.

  ;;move this reference image to be approximately centered
  method='F' ;for fourier cross-correlation
  subreg, gauss_cen, center_ref_line_smooth, sft, method=string(method)
  center_ref_line=shift_interp(center_ref_line, sft-1, spline=-0.5)
  subreg, gauss_cen, center_ref_cont_smooth, sft, method=string(method)
  center_ref_cont=shift_interp(center_ref_cont, sft-1, spline=-0.5)

  ;;if clip keyword not defined (not recommended), image size is full array (10254x512)
  if not keyword_set(clip) then clip=dim1

  ;;create blank arrays
  if not keyword_set(indiv) then begin
    cubedim=nims
  endif else begin
    cubedim=1
  endelse

  Line_smooth=dblarr(dim1, dim2, cubedim)
  Cont_smooth=dblarr(dim1, dim2, cubedim)
  Line_reg=dblarr(clip, clip, nims)
  Cont_reg=dblarr(clip,clip, nims)
  SDI_im=dblarr(clip, clip, nims)
  ;SDI_im2=dblarr(clip, clip, nims)
  peak_ratio=dblarr(nims)
  expt=dblarr(nims)
  Line_shift_arr=dblarr(2,nims)
  Cont_shift_arr=dblarr(2,nims)
  crops=dblarr(4,nims)

  ;;image loop
  for i=0, nims-1 do begin

    if keyword_set(indiv) then begin
      Line=readfits('./indiv/Line'+string(namestr)+string(i+1, format='(i04)')+'.fits', Linehead)
      Cont=readfits('./indiv/Cont'+string(namestr)+string(i+1, format='(i04)')+'.fits', Conthead)
      j=0
    endif else begin
      j=i
    endelse
    ;smooth out cosmic rays
    Line_smooth[*,*,j]=gauss_smooth(Line[*,*,j], width=5.)*1.
    Cont_smooth[*,*,j]=gauss_smooth(Cont[*,*,j], width=5.)*1.

    ;; calculate approximate offset from center of image
    subreg, center_ref_line, Line[*,*,j], Line_shift, method=string(method)
    subreg, center_ref_cont, Cont[*,*,j], Cont_shift, method=string(method)

    ;; record shifts
    Line_shift_arr[*,i]=Line_shift
    Cont_shift_arr[*,i]=Cont_shift

    ;if clip mod 2 eq 0 then offset=1 else offset=1.5
    Line[*,*,j]=shift_interp(Line[*,*,j], Line_shift-1, spline=-0.5)
    Cont[*,*,j]=shift_interp(Cont[*,*,j], Cont_shift-1, spline=-0.5)

    ;;rescale continuum image by scl if keyowrd is set
    if keyword_set(scl) then begin
      Cont[*,*,j]=rot(Cont[*,*,j],0.,scl,cubic=-0.5)
    endif

    ;;clip may overshoot edge of frame in some images. Adjust clip to be smaller in these cases
    clip=float(clip)

    ;;Line adjustments to clip
    if dim1/2-abs(Line_shift[0]) lt clip/2 then begin
      print, 'Line overshoots by', dim1/2-abs(Line_shift[0])-clip/2, 'pixels in x'
      Line_cropx=floor(dim1/2-abs(Line_shift[0]))
    endif else begin
      Line_cropx=clip/2
    endelse

    if dim2/2-abs(Line_shift[1]) lt clip/2 then begin
      print, 'Line overshoots by', dim2/2-abs(Line_shift[0])-clip/2, 'pixels in y'
      Line_cropy=floor(dim2/2-abs(Line_shift[1]))
    endif else begin
      Line_cropy=clip/2
    endelse

    ;;Cont adjustments to clip
    if dim1/2-abs(Cont_shift[0]) lt clip/2 then begin
      print, 'Cont overshoots by', dim1/2-abs(Line_shift[0])-clip/2, 'pixels in x'
      Cont_cropx=floor(dim1/2-abs(Line_shift[0]))
    endif else begin
      Cont_cropx=clip/2
    endelse
    if dim2/2-abs(Cont_shift[1]) lt clip/2 then begin
      print, 'Cont overshoots by', dim2/2-abs(Line_shift[0])-clip/2, 'pixels in y'
      Cont_cropy=floor(dim2/2-abs(Cont_shift[1]))
    endif else begin
      Cont_cropy=clip/2
    endelse
    ;stop
    ;record crops
    crops[*,i]=[Line_cropx, Line_cropy, Cont_cropx, Cont_cropy]

    Line_reg[clip/2-Line_cropx:clip/2+Line_cropx-1,clip/2-Line_cropy:clip/2+Line_cropy-1,i]= $
      Line[xcen-Line_cropx+1:xcen+Line_cropx,ycen-Line_cropy+1:ycen+Line_cropy,j]
    Cont_reg[clip/2-Cont_cropx:clip/2+Cont_cropx-1,clip/2-Cont_cropy:clip/2+Cont_cropy-1,i]= $
      Cont[xcen-Cont_cropx+1:xcen+Cont_cropx, ycen-Cont_cropy+1:ycen+Cont_cropy,j]

    print, 'processed image', i+1, '        of', nims

  endfor

  ;;if keyword is specified, run center of circular symmetry algorithm to refine centers
  if keyword_set(refine_cen) then begin
    print, 'refining centers - this will take a few minutes'
    
    ;;reread in original images to avoid interpolating multiple times
    Line=readfits('Line'+string(namestr)+'preproc.fits', Linehead)
    Cont=readfits('Cont'+string(namestr)+'preproc.fits', Conthead)

    ;;make medians of registered images
    Linemed=median(Line_reg, dim=3)
    Contmed=median(Cont_reg, dim=3)

    ;;make grid of possible centers to examine
    xr=indgen(51.)-50/2.
    yr=indgen(51.)-50/2.

    if keyword_set(rmax) then lim=rmax else lim=dim1/2.

    ;;compute center
    if keyword_set(mask) then begin
      center_circlesym, Linemed, xr, yr, lim, Line_xc, Line_yc, Line_grid, mask=mask
      center_circlesym, Contmed, xr, yr, lim, Cont_xc, Cont_yc, Cont_grid, mask=mask
    endif else begin
      center_circlesym, Linemed, xr, yr, lim, Line_xc, Line_yc, Line_grid
      center_circlesym, Contmed, xr, yr, lim, Cont_xc, Cont_yc, Cont_grid
    endelse

    ;;compute shifts from returned centers
    Line_censhift=[(clip-1)/2.-Line_xc,(clip-1)/2.-Line_yc]
    Cont_censhift=[(clip-1)/2.-Cont_xc,(clip-1)/2.-Cont_yc]

    print, 'center of circular symmetry for median Line image is', Line_xc, Line_yc
    print, 'shifting Line images by', Line_censhift
    print, 'center of circular symmetry for median Continuum image is', Cont_xc, Cont_yc
    print, 'shifting Continuum images by', Cont_censhift

    for i=0, nims-1 do begin
      ;;overwrite full cube with shifted versions to avoid memory issues
      ;;shifts are originally computed shifts + center refinement
      Line[*,*,i]=shift_interp(Line[*,*,i], Line_shift_arr[*,i]-1+Line_censhift, spline=-0.5)
      Cont[*,*,i]=shift_interp(Cont[*,*,i], Cont_shift_arr[*,i]-1+Cont_censhift, spline=-0.5)
      ;;read crops from before
      Line_cropx=crops[0,i]
      Line_cropy=crops[1,i]
      Cont_cropx=crops[2,i]
      Cont_cropy=crops[3,i]
      ;;crop images - will overwrite Line_reg from first pass at registration
      Line_reg[clip/2-Line_cropx:clip/2+Line_cropx-1,clip/2-Line_cropy:clip/2+Line_cropy-1,i]= $
        Line[xcen-Line_cropx+1:xcen+Line_cropx,ycen-Line_cropy+1:ycen+Line_cropy,i]
      Cont_reg[clip/2-Cont_cropx:clip/2+Cont_cropx-1,clip/2-Cont_cropy:clip/2+Cont_cropy-1,i]= $
        Cont[xcen-Cont_cropx+1:xcen+Cont_cropx, ycen-Cont_cropy+1:ycen+Cont_cropy,i]
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
