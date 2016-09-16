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

pro visao_reg, ref, clip=clip, flat=flat, fwhm=fwhm, sdi=sdi, indiv=indiv, scl=scl, stp=stp

  if keyword_set(flat) then namestr='_flat_' else namestr='_'

  ;; read in channel cubes from visao_separate_sdi, which are dark subtracted and flat fielded
  if not keyword_set(indiv) then begin
  Line=readfits('Line'+string(namestr)+'preproc.fits')
  Cont=readfits('Cont'+string(namestr)+'preproc.fits')
 
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

  ;;image loop
  for i=0, nims-1 do begin
  
   if keyword_set(indiv) then begin
    Line=readfits('./indiv/Line'+string(namestr)+string(i+1, format='(i04)')+'.fits')
    Cont=readfits('./indiv/Cont'+string(namestr)+string(i+1, format='(i04)')+'.fits')    
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

    ;; record shifts as a check
    Line_shift_arr[*,i]=Line_shift
    Cont_shift_arr[*,i]=Cont_shift

    Line[*,*,j]=shift_interp(Line[*,*,j], Line_shift-1, spline=-0.5)
    Cont[*,*,j]=shift_interp(Cont[*,*,j], Cont_shift-1, spline=-0.5)

    ;;rescale continuum image by scl if keyowrd is set
    if keyword_set(scl) then begin
          Cont[*,*,j]=rot(Cont[*,*,j],0.,scl,cubic=-0.5)
    endif

    ;;clip may overshoot edge of frame in some images. Adjust clip to be smaller in these cases

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
    Line_reg[clip/2-Line_cropx:clip/2+Line_cropx-1,clip/2-Line_cropy:clip/2+Line_cropy-1,i]= $
      Line[xcen-Line_cropx+1:xcen+Line_cropx,ycen-Line_cropy+1:ycen+Line_cropy,j]
    Cont_reg[clip/2-Cont_cropx:clip/2+Cont_cropx-1,clip/2-Cont_cropy:clip/2+Cont_cropy-1,i]= $
      Cont[xcen-Cont_cropx+1:xcen+Cont_cropx, ycen-Cont_cropy+1:ycen+Cont_cropy,j]

    ;;if clip is odd number, then need to shift by 0.5 to get to center of array. VisAO chip channels
    ;are even numbered, so it's expected that the center of the star will be between two pixels
    if clip/2. mod 1 gt 0 then begin
      Line_reg[*,*,i]=shift_interp(Line_reg[*,*,i], [0.5, 0.5], spline=-0.5)
      Cont_reg[*,*,i]=shift_interp(Cont_reg[*,*,i], [0.5, 0.5], spline=-0.5)
    endif

    print, 'processed image', i+1, '        of', nims

  endfor


    if keyword_set(clip) then begin
    writefits, 'Line_clip'+string(clip, format='(i03)')+string(namestr)+'reg.fits', Line_reg
    writefits, 'Cont_clip'+string(clip, format='(i03)')+string(namestr)+'reg.fits', Cont_reg
    endif else begin
      writefits, 'Line'+string(namestr)+'reg.fits', Line_reg
      writefits, 'Cont'+string(namestr)+'reg.fits', Cont_reg
    endelse
    
    if keyword_set(sdi) then begin
      if keyword_set(clip) then begin
      writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+string(namestr)+'_reg.fits', SDI_im
      ;writefits, 'SDI_sc'+string(1, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg.fits', SDI_im2
      endif else begin
        writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+string(namestr)+'_reg.fits', SDI_im
        ;writefits, 'SDI_sc'+string(1, format='(f05.2)')+'_reg.fits', SDI_im2 
      endelse
    endif

if keyword_set(stp) then  stop

print, 'registration complete'

end
