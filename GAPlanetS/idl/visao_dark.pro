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

  for i=0, n_elements(dark_imlist)-1 do begin
    darks[*,*,i]=readfits(dark_imlist[i], head, /silent)
    expt[i]=sxpar(head, 'EXPTIME')
    level_ch1[i]=median(darks[0:1023,0:511,i])
    level_ch2[i]=median(darks[0:1023,512:1023,i])
    ;print, i, level_ch1[i], level_ch2[i]
  endfor

;;compute and print some summary statistics
  if n_elements(uniq(number_formatter(expt,decimals=2))) eq 1 then begin
    if n_elements(dark_imlist) eq 1 then master_dark=darks else master_dark=median(darks, dim=3)
    writefits, 'dark_ch1medians.fits', level_ch1
    writefits, 'dark_ch2medians.fits', level_ch2
    print, 'median channel 1 dark level is', median(level_ch1), 'ADU and std deviation is', stddev(level_ch1), 'ADU'
    print, 'median channel 2 dark level is', median(level_ch2), 'ADU and std deviation is', stddev(level_ch2), 'ADU'
    mkhdr, hdr, master_dark
    sxaddpar, hdr, 'EXPTIME', expt[0]
    sxaddpar, hdr, 'NDARKS', n_elements(dark_imlist)
    sxaddpar, hdr, 'CH1_MED', median(level_ch1)
    sxaddpar, hdr, 'CH1_STD', stddev(level_ch1)
    sxaddpar, hdr, 'CH2_MED', median(level_ch2)
    sxaddpar, hdr, 'CH2_STD', stddev(level_ch2)
    writefits, 'master_dark.fits', master_dark, hdr
  endif else print, 'more than one exposure time in dark list - no dark created'

if keyword_set(writecube) then writefits, 'darks.fits', darks, hdr
  
if keyword_set(stp) then  stop

end
