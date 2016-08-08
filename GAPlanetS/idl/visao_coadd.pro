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
;  Modified 08-01-16 by KBF to reflect GAPlanetS infrastructure developments
;
;-

pro visao_coadd, fname, rotfname, wfefname, ncoadds, wfecut=wfecut, dir=dir

  im_cube=readfits(string(fname)+'.fits')
  dim1=(size(im_cube))[1]
  dim2=(size(im_cube))[2]
  nims=(size(im_cube))[3]
  nims_coadded=floor(nims/ncoadds)

  rotoff_cube=readfits(string(rotfname)+'.fits')
  wfe_cube=readfits(string(wfefname)+'.fits')

  ;;temporary arrays to hold individual images and params for each coadd
  toadd_im=dblarr(dim1, dim2, ncoadds)
  rotoff=dblarr(ncoadds)
  wfe=dblarr(ncoadds)

  ;;arrays to hold coadded images
  keep_matrix=dblarr(nims_coadded)
  im_coadded=dblarr(dim1, dim2, nims_coadded)
  rotoffs=dblarr(nims_coadded)

  for i=0, nims_coadded-1 do begin
    for j=0, ncoadds-1 do begin
      stat='image '+string(i*ncoadds+j+1)+' of '+string(nims)
      statusline, stat
      toadd_im[*,*,j]=im_cube[*,*,i*ncoadds+j]
      rotoff[j]=rotoff_cube[i*ncoadds+j]
      wfe[j]=wfe_cube[i*ncoadds+j]
    endfor
    im_coadded[*,*,i]=total(toadd_im, 3)
    rotoffs[i]=mean(rotoff)
    if keyword_set(wfecut) then begin
      if mean(wfe) lt wfecut then keep_matrix[i]=1
    endif
    ;stop
  endfor
  ;stop

  ;;write into directory if keyword is set
  if keyword_set(dir) then begin
    spawn, 'mkdir '+string(dir)
    fname='./'+string(dir)+'/'+string(fname)
  endif

  ;;write arrays
  writefits, string(fname)+'_coadd'+string(ncoadds, format='(i02)')+'.fits', im_coadded

  ;;make wfe cut if desired

  if keyword_set(wfecut) then begin
    kept=where(keep_matrix gt 0)
    ;print, kept
    nkept=n_elements(kept)
    rotoff_kept=dblarr(nkept)
    im_wfecut=dblarr(dim1,dim2,nkept)

    for k=0, nkept-1 do begin
      index=kept[k]
      ;print, index
      rotoff_kept[k]=rotoffs[index]
      im_wfecut[*,*,k]=im_coadded[*,*,index]
    endfor

    writefits, string(fname)+'_coadd'+string(ncoadds, format='(i02)')+'_wfe'+string(wfecut, format='(f05.1)')+'.fits', im_wfecut
    writefits, string(rotfname)+'_coadd'+string(ncoadds, format='(i02)')+'_wfe'+string(wfecut, format='(f05.1)')+'.fits', rotoff_kept

    requad_rotoffs=requad_angles(rotoff_kept)
    requad_rotoffs=requad_angles(rotoff_kept)

    print, 'total coadded images with wfe <', wfecut, '=', total(keep_matrix), ' out of', nims_coadded, 'total'
    print, 'total rotation =', max(requad_rotoffs) - min(requad_rotoffs), 'out of', totrot
   ; print, 'total integration = ', expt*ncoadds*nkept/60, ' min out of', expt*ncoadds*nims_coadded/60, ' min total'

  endif
  ;stop
end
