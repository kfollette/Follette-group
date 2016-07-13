;+
; NAME: visao_coadd
;
; PURPOSE:
;  coadd images aligned with visao_reg and make a wavefront error cut, if desired.
;
; INPUTS:
;  ncoadds : number of images to coadd
;
; INPUT KEYWORDS:
;  wfecut  : maximum wavefront error to allow into final dataset
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
;  Written 02-18-16 by Kate Follette, kbf@stanford.edu
;
;-

pro visao_coadd_new, ncoadds, clip=clip, wfecut=wfecut, sdi=sdi

  visao_inventory, sci_imlist, dark_imlist, flat_imlist, rotoff_sciims, filt, wfe=wfe, mag1=mag1, totrot=totrot

  if keyword_set(clip) then begin
    imcube=readfits('Line_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics.fits')
    contcube=readfits('Cont_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics.fits')
    if keyword_set(sdi) then begin
      sdi1cube=readfits('SDI_sc'+string(sdi, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg_circsym_nocosmics.fits')
      sdi2cube=readfits('SDI_sc'+string(1, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg_circsym_nocosmics.fits')
    endif
  endif else begin
    imcube=readfits('Line_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics.fits')
    contcube=readfits('Cont_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics.fits')
    if keyword_set(sdi) then begin
      sdi1cube=readfits('SDI_sc'+string(sdi, format='(f05.2)')+'_reg_circsym_nocosmics.fits')
      sdi2cube=readfits('SDI_sc'+string(1, format='(f05.2)')+'_reg_circsym_nocosmics.fits')
    endif
  endelse

  dim1=(size(imcube))[1]
  dim2=(size(imcube))[2]
  head_dim=n_elements(header)

  ;;coadd aligned images and write into directory
  nims=(size(imcube))[3]
  nims_coadded=floor(nims/ncoadds)
  ;im_coadded=dblarr(dim1, dim2)
  toadd_Line=dblarr(dim1, dim2, ncoadds)
  toadd_Cont=dblarr(dim1, dim2, ncoadds)
  rotoff_list=readfits('rotoff_nocosmics.fits')
  wfe_list=readfits('avgwfe_nocosmics.fits')
  rotoff=dblarr(ncoadds)
  wfe=dblarr(ncoadds)
  rotoffs=dblarr(nims_coadded)
  wfes=dblarr(nims_coadded)

  exptime_list=readfits('exptime_nocosmics.fits')
  keep_matrix=dblarr(nims_coadded)
  Line_im_coadded=dblarr(dim1, dim2, nims_coadded)
  Cont_im_coadded=dblarr(dim1, dim2, nims_coadded)
  rotoffs=dblarr(nims_coadded)
  if keyword_set(sdi) then begin
    toadd_SDI1=dblarr(dim1, dim2, ncoadds)
    toadd_SDI2=dblarr(dim1, dim2, ncoadds)
    SDI_im1_coadded=dblarr(dim1, dim2, nims_coadded)
    SDI_im2_coadded=dblarr(dim1, dim2, nims_coadded)
  endif
  header= strarr(102,nims)

  for i=0, nims_coadded-1 do begin
    for j=0, ncoadds-1 do begin
      imno=(i*ncoadds+j)
      stat='image '+string(imno+1)+' of '+string(nims)
      statusline, stat
      toadd_Line[*,*,j]=imcube[*,*,imno]
      toadd_Cont[*,*,j]=contcube[*,*,imno]
      if keyword_set(sdi) then begin
        toadd_SDI1[*,*,j]=sdi1cube[*,*,imno]
        toadd_SDI2[*,*,j]=sdi2cube[*,*,imno]
      endif
      rotoff[j]=rotoff_list[imno]
      wfe[j]=wfe_list[imno]
      expt=exptime_list[j]
    endfor
    Line_im_coadded[*,*,i]=total(toadd_Line, 3)
    Cont_im_coadded[*,*,i]=total(toadd_Cont,3)
    if keyword_set(sdi) then begin
      SDI_im1_coadded[*,*,i]=total(toadd_SDI1, 3)
      SDI_im2_coadded[*,*,i]=total(toadd_SDI2,3)
    endif

    rotoffs[i]=mean(rotoff)
    wfes[i]=mean(wfe)

    if keyword_set(wfecut) then begin
      if mean(wfe) lt wfecut then keep_matrix[i]=1
    endif
  endfor
  writefits, 'rotoff_coadd'+string(ncoadds, format='(i02)')+'.fits', rotoffs
  writefits, 'Line_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics_coadd'+string(ncoadds, format='(i02)')+'.fits', Line_im_coadded
  writefits, 'Cont_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics_coadd'+string(ncoadds, format='(i02)')+'.fits', Cont_im_coadded

  if keyword_set(sdi) then begin
    writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+'_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics_coadd'+string(ncoadds, format='(i02)')+'.fits', SDI_im1_coadded
    writefits, 'SDI_sc'+string(1, format='(f05.2)')+'_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics_coadd'+string(ncoadds, format='(i02)')+'.fits' , SDI_im2_coadded
  endif
  ;;now cut on wavefront error
  if keyword_set(wfecut) then begin
    kept=where(keep_matrix gt 0)
    ;print, kept
    nkept=n_elements(kept)
    rotoff_kept=dblarr(nkept)
    wfe_kept=dblarr(nkept)
    Line_wfecut=dblarr(dim1, dim2, nkept)
    Cont_wfecut=dblarr(dim1, dim2, nkept)
    if keyword_set(sdi) then begin
      SDI_im1_wfecut=dblarr(dim1, dim2, nkept)
      SDI_im2_wfecut=dblarr(dim1, dim2, nkept)
    endif

    for k=0, nkept-1 do begin
      index=kept[k]
      ;print, index
      rotoff_kept[k]=rotoffs[index]
      wfe_kept[k]=wfes[index]
      Line_wfecut[*,*,k]=Line_im_coadded[*,*,index]
      Cont_wfecut[*,*,k]=Cont_im_coadded[*,*,index]
      if keyword_set(sdi) then begin
        SDI_im1_wfecut[*,*,k]=SDI_im1_coadded[*,*,index]
        SDI_im2_wfecut[*,*,k]=SDI_im2_coadded[*,*,index]
      endif
    endfor

    writefits, 'rotoff_coadd'+string(ncoadds, format='(i02)')+'.fits', rotoff_kept

    writefits, 'Line_clip'+string(clip,format='(i03)')+'_reg_circsym_coadd'+string(ncoadds, format='(i02)')+'_wfe'+string(wfecut,format='(f05.1)')+'.fits', Line_wfecut
    writefits, 'Cont_clip'+string(clip,format='(i03)')+'_reg_circsym_coadd'+string(ncoadds, format='(i02)')+'_wfe'+string(wfecut,format='(f05.1)')+'.fits', Cont_wfecut

    if keyword_set(sdi) then begin
      writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+'_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics_coadd'+string(ncoadds, format='(i02)')+'_wfe'+string(wfecut,format='(f05.1)')+'.fits', SDI_im1_wfecut
      writefits, 'SDI_sc'+string(1, format='(f05.2)')+'_clip'+string(clip,format='(i03)')+'_reg_circsym_nocosmics_coadd'+string(ncoadds, format='(i02)')+'_wfe'+string(wfecut,format='(f05.1)')+'.fits', SDI_im2_wfecut
    endif

    requad_rotoffs=requad_angles(rotoff_kept)
    requad_rotoffs=requad_angles(rotoff_kept)

    print, 'total coadded images with wfe <', wfecut, '=', total(keep_matrix), ' out of', nims_coadded, 'total'
    print, 'total rotation =', max(requad_rotoffs) - min(requad_rotoffs), 'out of', totrot
    print, 'total integration = ', expt*ncoadds*nkept/60, ' min out of', expt*ncoadds*nims_coadded/60, ' min total'
  endif
  ;stop
end