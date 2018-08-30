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

pro visao_removecosmics, fname, namestr, stp=stp

  ;;read in image (should be SDI)
  im=readfits(string(fname)+'.fits')
  zdim=(size(im))[3]
  print, zdim

  ds9_imselect, im, index=idx

  writefits, string(fname)+'_nocosmics.fits', im
  writefits, namestr+'cosmics.fits', idx

  ;; cull rotoff cube as well
  rotoffs=readfits('rotoff_preproc.fits')
  rotoffs_noc=dblarr(n_elements(idx))

  for i=0, n_elements(idx)-1 do begin
    rotoffs_noc[i]=rotoffs[idx[i]]
  endfor

  writefits, 'rotoff_no'+namestr+'cosmics.fits', rotoffs_noc

  if keyword_set(stp) then  stop

end
