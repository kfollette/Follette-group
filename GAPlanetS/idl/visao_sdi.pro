;+
; NAME: visao_sdi
;
; PURPOSE:
;  apply scalings (either individually or in bulk) and subtract continuum from line to create
;  an SDI image
;
; INPUTS:
;
; INPUT KEYWORDS:
; indiv: apply calculated scalings one by one
; scale: a single scale to apply to all of the images
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
;  Written 07-26-16 by Kate Follette, kbf@stanford.edu
;
;-

pro visao_sdi, indiv=indiv, scale=scale, clip=clip, flat=flat, stp=stp

  if keyword_set(flat) then namestr='_flat_' else namestr='_'

  if keyword_set(clip) then begin
    Line=readfits('Line_clip'+string(clip,format='(i03)')+string(namestr)+'reg_circsym.fits')
    Cont=readfits('Cont_clip'+string(clip,format='(i03)')+string(namestr)+'reg_circsym.fits')
  endif else begin
    Line=readfits('Line'+string(namestr)+'reg_circsym.fits')
    Cont=readfits('Cont'+string(namestr)+'reg_circsym.fits')
  endelse

  if keyword_set(indiv) then begin
    scl = readfits('scale_factors.fits')
    sdi = 'indiv'
  endif

  if keyword_set(scale) then begin
    scl=dblarr((size(Line))[3])+scale
    sdi=string(scale, format='(f05.2)')
  endif

  if not keyword_set(scale) and not keyword_set(indiv) then begin
    print, 'either scale or indiv keyword must be set'
    stop
  endif

  if n_elements(scl) gt 1 then begin
    print, 'scaling continuum images by individually computed ratios and subtracting'
  endif else begin
    print, 'scaling continuum images by', scale, 'and subtracting'
  endelse

  SDIim=dblarr(clip, clip, n_elements(scl))

  for i =0, n_elements(scl) -1 do begin
    status='image number'+string(i)+'  of'+string(n_elements(scl))
    statusline, status, 0
    SDIim[*,*,i] = Line[*,*,i] - scl[i]*Cont[*,*,i]
  endfor

  ;; release Line and continuum cubes from memory
  delvar, Line, Cont

  writefits, 'SDI_'+string(sdi)+'_clip'+string(clip,format='(i03)')+string(namestr)+'reg_circsym.fits', SDIim

if keyword_set(stp) then stop

print, 'done creating SDI image cube'

end
