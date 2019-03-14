pro visao_dark, dark_imlist, master_dark, writecube=writecube

  visao_inventory, sci_imlist, dark_imlist, flat_imlist, rotoff_sciims, filt, wfe=wfe, stp=stp

  dummy_im=readfits(dark_imlist[0])
  dim1=(size(dummy_im[*,*]))[1]
  dim2=(size(dummy_im[*,*]))[2]

  print, 'Creating master dark from', n_elements(dark_imlist), ' dark frames'
  darks=dblarr(dim1,dim2,n_elements(dark_imlist))
  expt=dblarr(n_elements(dark_imlist))
  level_ch1=dblarr(n_elements(dark_imlist))
  level_ch2=dblarr(n_elements(dark_imlist))
  gain=strarr(n_elements(dark_imlist))
  
  for i=0, n_elements(dark_imlist)-1 do begin
    darks[*,*,i]=readfits(dark_imlist[i], head, /silent)
    expt[i]=sxpar(head, 'EXPTIME')
    gain[i]=sxpar(head, 'V47GAIN')
    level_ch1[i]=median(darks[0:1023,0:511,i])
    level_ch2[i]=median(darks[0:1023,512:1023,i])
    ;print, i, level_ch1[i], level_ch2[i]
  endfor

  gains_sorted=gain[sort(gain)]
  gains_uniq = gains_sorted[uniq(gains_sorted)]
  print, gains_uniq
  
  ;;loop over gains
  for i=0, n_elements(gains_uniq) - 1 do begin
    gaindx = where(gain eq gains_uniq[i])
    darks_singlegain = darks[*,*,gaindx]
    print, 'writing files for ', n_elements(gaindx), ' darks with gain ', gains_uniq[i]
    ;;compute and print some summary statistics
    if n_elements(uniq(number_formatter(expt,decimals=2))) eq 1 then begin
      if (size(darks_singlegain))[3] eq 1 then master_dark=darks_singlegain else master_dark=median(darks_singlegain, dim=3)
      writefits, 'dark_'+strcompress(gains_uniq[i], /REMOVE_ALL)+'_ch1medians.fits', level_ch1[gaindx]
      writefits, 'dark_'+strcompress(gains_uniq[i], /REMOVE_ALL)+'_ch2medians.fits', level_ch2[gaindx]
      print, 'median channel 1 dark level for gain ', gains_uniq[i]+' is ', median(level_ch1[gaindx]), 'ADU and std deviation is', stddev(level_ch1[gaindx]), 'ADU'
      print, 'median channel 2 dark level for gain ', gains_uniq[i]+' is ', median(level_ch2[gaindx]), 'ADU and std deviation is', stddev(level_ch2[gaindx]), 'ADU'
      mkhdr, hdr, master_dark
      sxaddpar, hdr, 'EXPTIME', expt[0]
      sxaddpar, hdr, 'NDARKS', n_elements(dark_imlist)
      sxaddpar, hdr, 'CH1_MED', median(level_ch1[gaindx])
      sxaddpar, hdr, 'CH1_STD', stddev(level_ch1[gaindx])
      sxaddpar, hdr, 'CH2_MED', median(level_ch2[gaindx])
      sxaddpar, hdr, 'CH2_STD', stddev(level_ch2[gaindx])
      writefits, 'master_dark_'+strcompress(gains_uniq[i], /REMOVE_ALL)+'.fits', master_dark, hdr
    endif else print, 'more than one exposure time in dark list - no dark created'

    if keyword_set(writecube) then writefits, 'darks.fits', darks, hdr


  endfor
  if keyword_set(stp) then  stop

end
