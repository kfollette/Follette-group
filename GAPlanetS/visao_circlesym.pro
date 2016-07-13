pro visao_circlesym, Line_cent, Cont_cent, msk=msk, clip=clip, sdi=sdi, rmax=rmax

  ;;finds center of circular symmetry of median combinations of registered images, and shifts this center to the center of the image cube
  ;;6-30-2016 KBF modification - do line and continuum separately, then SDI so don't hold too many large arrays in memory at once
  ;; also adding more print statements and removing /fits keyword so will write output to disk by default

  if keyword_set(clip) then Line=readfits('Line_clip'+string(clip,format='(i03)')+'_reg.fits') $
  else Line=readfits('Line_reg.fits')

  Linemed=median(Line, dim=3)

  if keyword_set(msk) then begin
    print, 'constructing mask with radius of', msk, ' pixels'
    mkmask, (size(Linemed))[1], (size(Linemed))[2], msk, mask
  endif

  dim1=(size(Line))[1]
  dim2=(size(Line))[2]

  xr=indgen(dim1/2.+1)-dim1/4.
  yr=indgen(dim2/2.+1)-dim2/4.

  if keyword_set(rmax) then lim=rmax else lim=dim1/2.

  print, 'calculating center of circular symmetry for median Line image'
  if keyword_set(msk) then begin
    center_circlesym, Linemed, xr, yr, lim, Line_xc, Line_yc, Line_grid, mask=mask
  endif else begin
    center_circlesym, Linemed, xr, yr, lim, Line_xc, Line_yc, Line_grid, mask=mask
  endelse

  Line_shift=[(dim1-1)/2.-Line_xc,(dim2-1)/2.-Line_yc]

  print, 'center of circular symmetry for median Line image is', Line_xc, Line_yc
  print, 'shifting all Line images by', Line_shift
  nims=(size(Line))[3]

  Line_cent=dblarr(dim1, dim2, nims)

  for i=0, nims-1 do begin
    print, 'shifting Line image', i+1, '   of', nims
    Line_cent[*,*,i]=shift_interp(Line[*,*,i], Line_shift, spline=-0.5)
  endfor

  print, 'writing centered line image cube'
  if keyword_set(clip) then begin
    writefits, 'Line_clip'+string(clip, format='(i03)')+'_reg_circsym.fits', Line_cent
  endif else begin
    writefits, 'Line_reg_circsym.fits', Line_cent
  endelse

  delvar, Line_cent, Linemed
  ;;;; release the Line cubes from memory and do the same thing all over again for the continuum

  if keyword_set(clip) then Cont=readfits('Cont_clip'+string(clip,format='(i03)')+'_reg.fits') $
  else Cont=readfits('Cont_reg.fits')

  Contmed=median(Cont, dim=3)

  if keyword_set(msk) then begin
    print, 'constructing mask with radius of', msk, ' pixels'
    mkmask, (size(Linemed))[1], (size(Linemed))[2], msk, mask
  endif

  print, 'calculating center of circular symmetry for median Continuum image'
  if keyword_set(msk) then begin
    center_circlesym, Contmed, xr, yr, lim, Cont_xc, Cont_yc, Cont_grid, mask=mask
  endif else begin
    center_circlesym, Contmed, xr, yr, lim, Cont_xc, Cont_yc, Cont_grid, mask=mask
  endelse

  Cont_shift=[(dim1-1)/2.-Cont_xc,(dim2-1)/2.-Cont_yc]

  print, 'center of circular symmetry for median Continuum image is', Cont_xc, Cont_yc
  print, 'shifting all Continuum images by', Cont_shift
  nims=(size(Line))[3]

  Cont_cent=dblarr(dim1, dim2, nims)

  for i=0, nims-1 do begin
    print, 'shifting Continuum image', i+1, '   of', nims
    Cont_cent[*,*,i]=shift_interp(Cont[*,*,i], Cont_shift, spline=-0.5)
  endfor

  print, 'writing centered continuum image cube'
  if keyword_set(clip) then begin
    writefits, 'Cont_clip'+string(clip, format='(i03)')+'_reg_circsym.fits', Cont_cent
  endif else begin
    writefits, 'Cont_reg_circsym.fits', Cont_cent
  endelse

  ;;now release continuum cubes from memory and shift SDI
  delvar, Cont_cent, Contmed

  if keyword_set(sdi) then begin
    if keyword_set(clip) then begin
      SDI1=readfits('SDI_sc'+string(sdi, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg.fits')
      ;SDI2=readfits('SDI_sc'+string(1, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg.fits')
    endif else begin
      SDI1=readfits('SDI_sc'+string(sdi, format='(f05.2)')+'_reg.fits')
      ; SDI2=readfits('SDI_sc'+string(1, format='(f05.2)')+'_reg.fits')
    endelse
    SDI1_cent=dblarr(dim1,dim2,nims)
    ;SDI2_cent=dblarr(dim1,dim2,nims)

    for i=0, nims-1 do begin
      print, 'shifting SDI image', i+1, '   of', nims
      SDI1_cent[*,*,i]=shift_interp(SDI1[*,*,i], avg_shift, spline=-0.5)
      ;SDI2_cent[*,*,i]=shift_interp(SDI2[*,*,i], avg_shift, spline=-0.5)
    endfor

    print, 'writing centered SDI image cube'
    if keyword_set(clip) then begin
      writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg_circsym.fits',SDI1_cent
      ;writefits, 'SDI_sc'+string(1, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg_circsym.fits',SDI2_cent
    endif else begin
      writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+'_reg_circsym.fits',SDI1_cent
      ;writefits, 'SDI_sc'+string(1, format='(f05.2)')+'_reg_circsym.fits',SDI2_cent
    endelse
  endif
  stop
end