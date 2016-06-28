pro visao_cube, nims, unsharp=unsharp, name=name, sdi=sdi, subdir=subdir

  ;;keyword sumims to do coadds rather than median combine when binning
  ;;note: coadds and binning here don't really make much sense. Better to coadd with visao_coadd, then make any wfe cuts, then cube

  ;; 4/11/16 - binning is defunct (replaced by visao_coadd), so removed - old version saved in visao_cube_old.pro. Also added sdi cubes. Unsharp now takes a value
  ;; as in unsharp=6 is a 6-pixel unsharp mask

  if not keyword_set(subdir) then subdir='aligned'

  x=readfits(string(subdir)+'/Ha_0001.fits')
  xdim=(size(x[*,*]))[1]
  ydim=(size(x[*,*]))[2]

  Ha_ims=dblarr(xdim,ydim,nims)
  Cont_ims=dblarr(xdim,ydim,nims)
  Ha_ims_smoot=dblarr(xdim,ydim,nims)
  Cont_ims_smoot=dblarr(xdim,ydim,nims)
  if keyword_set(sdi) then begin
    SDI_ims=dblarr(xdim,ydim,nims)
    SDI_ims2=dblarr(xdim,ydim,nims)
    SDI_ims_smoot=dblarr(xdim,ydim,nims)
    SDI_ims2_smoot=dblarr(xdim,ydim,nims)
  endif
  
  rotoff=dblarr(nims)

  for x=0, nims-1 do begin
    ;;insert statement skipping frames
    Ha_ims[*,*,x]=readfits(string(subdir)+'/Ha_'+string(x+1, format='(i04)')+'.fits', head, /silent)
    Cont_ims[*,*,x]=readfits(string(subdir)+'/Cont_'+string(x+1, format='(i04)')+'.fits', /silent)
    if keyword_set(sdi) then begin
      SDI_ims[*,*,x]=readfits(string(subdir)+'/SDI_sc'+string(sdi, format='(f05.2)')+'_'+string(x+1, format='(i04)')+'.fits', /silent)
      SDI_ims2[*,*,x]=readfits(string(subdir)+'/SDI_sc01.00_'+string(x+1, format='(i04)')+'.fits', /silent)
    endif
    if keyword_set(unsharp) then begin
      Ha_ims_smoot[*,*,x]=Ha_ims[*,*,x]-gauss_smooth(Ha_ims[*,*,x],unsharp)
      Cont_ims_smoot[*,*,x]=Cont_ims[*,*,x]-gauss_smooth(Cont_ims[*,*,x],unsharp)
      if keyword_set(sdi) then begin
        SDI_ims_smoot[*,*,x]=SDI_ims[*,*,x]-gauss_smooth(SDI_ims[*,*,x],unsharp)
        SDI_ims2_smoot[*,*,x]=SDI_ims2[*,*,x]-gauss_smooth(SDI_ims2[*,*,x],unsharp)
      endif
    endif
    rotoff[x]=sxpar(head, 'ROTOFF')
  endfor

  writefits, 'Line_'+string(name)+'.fits', Ha_ims
  writefits, 'Cont_'+string(name)+'.fits', Cont_ims
  if keyword_set(sdi) then begin
    writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+'_'+string(name)+'.fits', SDI_ims
    writefits, 'SDI_sc01.00_'+string(name)+'.fits', SDI_ims2
  endif
  if keyword_set(unsharp) then begin
    writefits, 'Line_'+string(name)+'_unsharp.fits', Ha_ims_smoot
    writefits, 'Cont_'+string(name)+'_unsharp.fits', Cont_ims_smoot
    if keyword_set(sdi) then begin
      writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+'_'+string(name)+'_unsharp.fits', SDI_ims_smoot
      writefits, 'SDI_sc01.00_'+string(name)+'_unsharp.fits', SDI_ims2_smoot
    endif
  endif
  writefits, 'rotoff_'+string(name)+'.fits', rotoff

end