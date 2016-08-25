;+
; NAME: visao_removecosmics
;
; PURPOSE:
;  remove cosmic rays from a designated image cube, and apply to other cubes
;  in the same directory if desired
;
; INPUTS:
; fname: string file name of image cube that will be used to identify cosmic rays
;
; INPUT KEYWORDS:
; all: cull the same image slices from all other cubes in the same directory
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
;  Written 2015 by Kate Follette, kbf@stanford.edu
;  07-29-2016 KBF. Revied to read a more generic set of image cubes.
;     Added all keyword and genericized to find and appy to any circsym cube
;-

pro visao_removecosmics, fname, all=all, nantest=nantest, stp=stp

  ;;read in image (should be SDI)
  im=readfits(string(fname)+'.fits')
  zdim=(size(im))[3]
  print, zdim

  if keyword_set(nantest) then begin
    lowpix=where(im lt -2000)
    print, n_elements(lowpix), 'low pixels'
    highpix=where(im gt 2000)
    print, n_elements(highpix), 'high pixels'
    im[lowpix]='NaN'
    im[highpix]='NaN'
  endif

  ds9_imselect, im, index=idx

  writefits, string(fname)+'_nocosmics.fits', im
  writefits, 'cosmic_arr.fits', idx
  ;idx = readfits('cosmic_arr_SDIcube.fits')

  ;;apply same image culling to all 
  if keyword_set(all) then begin
    rotoffs=readfits('rotoff_preproc.fits')
    wfe=readfits('avgwfe_preproc.fits')
    exptime=readfits('exptime_preproc.fits')

    rotoffs_noc=dblarr(n_elements(idx))
    wfe_noc=dblarr(n_elements(idx))
    exptime_noc=dblarr(n_elements(idx))

    for i=0, n_elements(idx)-1 do begin
      rotoffs_noc[i]=rotoffs[idx[i]]
      wfe_noc[i]=wfe[idx[i]]
      exptime_noc[i]=exptime[idx[i]]
    endfor

    writefits, 'rotoff_nocosmics.fits', rotoffs_noc
    writefits, 'avgwfe_nocosmics.fits', wfe_noc
    writefits, 'exptime_nocosmics.fits', exptime_noc

    ;;find all circsym processed cubes in directory
    srchstr = strcompress('*circsym.fits', /remove)
    fnames = file_search(srchstr)

    for i=0, n_elements(fnames)-1 do begin
      ;;check hasn't already been culled
      img=readfits(fnames[i])
      xdim=(size(img))[1]
      ydim=(size(img))[2]
      zdim=(size(img))[3]
      img_culled=dblarr(xdim,ydim,n_elements(idx))
      if n_elements(idx) lt zdim then begin
        img_culled[*,*,i]=img[*,*,idx[i]]
        fstring=strsplit(fnames[i],'.',/extract)
        writefits, string(fstring[0])+'_nocosmics.fits', img
      endif else begin
        print, 'image', string(fnames[i]), 'has only', zdim, 'images in its cube.'
      endelse
    endfor
    
  endif

if keyword_set(stp) then  stop

end
