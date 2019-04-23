pro visao_inventory, sci_imlist, dark_imlist, flat_imlist, rotoff_sciims, filt, sciims, wfe=wfe, mag1=mag1, stp=stp, totrot=totrot, dir=dir

  if not keyword_set(dir) then dir='raw'

  visao_getimtypes, fnames, imtypes, exptime=exptime, vfw3posn=vfw3posn, vfw2posn=vfw2posn, aoloopst=aoloopst, $
    rotoff=rotoff, avgwfe=avgwfe, region=region, am=am, gain=gain, mag1fwhm=mag1fehm, dimmfwhm=dimmfwhm, subdir=dir

  sdx = where( imtypes eq 0 and aoloopst eq 1, count1)
  if count1 gt 0 then print, 'Found ', n_elements(sdx), ' closed loop science frames.' else print, 'no science frames'

  ddx = where( imtypes eq 2 , count2 )
  if count2 gt 0 then print, 'Found ', n_elements(ddx), ' darks.' else print, 'no dark frames'

  fdx = where( imtypes eq 4 , count3 )
  if count3 gt 0 then print, 'Found ', n_elements(fdx), ' flats.' else print, 'no flat frames'

  list=(intarr(4))

  exp_sort=exptime[sort(exptime)]
  exp_sort=float(number_formatter(exp_sort[*],decimals=2))
  exp_use=exp_sort[uniq(exp_sort)]
  if n_elements(exp_use) gt 1 then print, 'Warning - more than one exposure time in this image set - separate before proceeding'

  if keyword_set(wfe) then wfemax=wfe else wfemax=1000.
  if keyword_set(mag1) then mag1max=mag1 else mag1max=1000.

  for l=0, n_elements(exp_use)-1 do begin ;exposure time loop
    print, 'l', l

    ;;select only closed loop images with same filter and exposure time, and with wfe cutoff if /WFE keyword is set

    if((where(strmatch(vfw3posn,'*alpha*')))[0] ne -1) then begin
      Haims = where(strmatch(vfw3posn, '*alpha*') and ( imtypes eq 0 ) and $
        ( aoloopst eq 1 ) and $
        (number_formatter(exptime,decimals=2) eq exp_use[l]) and (avgwfe lt wfemax) )
      if Haims[0] ne -1 then print, 'Found ', n_elements(Haims), 'closed loop H alpha with exposure time', exp_use[l], $
        'and wfe<', wfemax
      list[0]=1
    endif

    if((where(strmatch(vfw3posn,'*S II*')))[0] ne -1) then begin
      SIIims = where(strmatch(vfw3posn, '*S II*') and ( aoloopst eq 1 ) and ( imtypes eq 0 ) and $
        (exptime eq exp_use[l]) and (avgwfe lt wfemax) )
      if SIIims[0] ne -1 then print, 'Found ', n_elements(SIIims), 'closed loop [SII] with exposure time', exp_use[l]
      list[1]=1
    endif

    if((where(strmatch(vfw3posn,'*O I*')))[0] ne -1) then begin
      OIims = where( strmatch(vfw3posn, '*O I*') and ( aoloopst eq 1 ) and ( imtypes eq 0 ) and $
        (exptime eq exp_use[l]) and (avgwfe lt wfemax) )
      if OIims[0] ne -1 then print, 'Found ', n_elements(OIims), 'closed loop [OI] with exposure time', exp_use[l]
      list[2]=1
    endif

    if((where(strmatch(vfw2posn,'*z*')))[0] ne -1) then begin
      zpims = where(strmatch(vfw2posn, '*z*') and ( imtypes eq 0 ) and ( aoloopst eq 1 ) and $
        (exptime eq exp_use[l]) and (avgwfe lt wfemax) )
      if zpims[0] ne -1 then print, 'Found ', n_elements(Haims), 'closed loop zp with exposure time', exp_use[l], $
        'and wfe<', wfemax
      list[3]=1
    endif


  endfor
  ;stop
  nfilt=where(list ne 0)
  if n_elements(nfilt) gt 1 then print, 'more than one SDI filter in this image set - separate before proceeding'

  if list[0] eq 1 then sciims=Haims
  if list[0] eq 1 then filt='Ha'
  if list[1] eq 1 then sciims=SIIims
  if list[1] eq 1 then filt='SII'
  if list[2] eq 1 then sciims=OIims
  if list[2] eq 1 then filt='OI'
  if list[3] eq 1 then sciims=zpims
  if list[3] eq 1 then filt='zp'
  if n_elements(nfilt) gt 1 then sciims='NaN'
  ;stop
  if count1 ne 0 then begin
    rotoff_sciims=rotoff[sciims]
    rotoff_cont=requad_angles(rotoff_sciims)
    totrot=max(rotoff_cont)-min(rotoff_cont)
    print, 'total rotation of this dataset is   ', totrot, '   degrees'
    gain_sciims = gain[sciims]
  endif

  if count1 ne 0 then sci_imlist=fnames[sciims]
  if count2 ne 0 then dark_imlist=fnames[ddx]
  if count2 ne 0 then gain_darks = gain[ddx]
  if count3 ne 0 then flat_imlist=fnames[fdx]

  gains_sorted=gain[sort(gain)]
  gains_uniq = gains_sorted[uniq(gains_sorted)]
  print, gains_uniq
  if n_elements(gains_uniq) gt 1 then print, "WARNING: Multiple gains"

  if keyword_set(stp) then stop

  return

end
