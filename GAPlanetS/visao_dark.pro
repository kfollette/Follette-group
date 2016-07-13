pro visao_dark, dark_imlist, master_dark

visao_inventory, sci_imlist, dark_imlist, flat_imlist, rotoff_sciims, filt, wfe=wfe

  dummy_im=readfits(dark_imlist[0])
  dim1=(size(dummy_im[*,*]))[1]
  dim2=(size(dummy_im[*,*]))[2]
  
  print, 'Creating master dark from', n_elements(dark_imlist), ' dark frames'
  darks=dblarr(dim1,dim2,n_elements(dark_imlist))
  expt=dblarr(n_elements(dark_imlist))

  for i=0, n_elements(dark_imlist)-1 do begin
    darks[*,*,i]=readfits(dark_imlist[i], head, /silent)
    expt[i]=sxpar(head, 'EXPTIME')
  endfor


  if n_elements(uniq(expt)) eq 1 then begin
    if n_elements(dark_imlist) eq 1 then master_dark=darks else master_dark=median(darks, dim=3)
      mkhdr, hdr, master_dark
      sxaddpar, hdr, 'EXPTIME', expt[0]
      sxaddpar, hdr, 'NDARKS', n_elements(dark_imlist)
    writefits, 'master_dark.fits', master_dark, hdr
  endif else print, 'more than one exposure time in dark list - no dark created'
  
  stop
end