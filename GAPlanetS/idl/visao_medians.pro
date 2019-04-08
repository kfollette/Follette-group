  ;+
  ; NAME: visao_medians
  ;
  ; PURPOSE:
  ;  make median combination of registered visao SDI images
  ;
  ; INPUTS:
  ;
  ; INPUT KEYWORDS:
  ;  reg : set if making combo of coarsely registered images
  ;  circsymreg : set if making combo of finely registered images
  ;  clip : set to clip size so finds the right files
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
  ;  Written 2019-02-12 by Kate Follette kfollette@amherst.edu

pro visao_medians, reg=reg, circsymreg=circsymreg, medreg=medreg, ref=ref, clip=clip, stp=stp

  if keyword_set(reg) then namestr='_reg'
  if keyword_set(circsymreg) then namestr='_circsymreg'
  if keyword_set(medreg) then namestr='_reg_medcircsym'
  if keyword_set(ref) then namestr='_reg_refined'

  lineim = readfits('Line_clip'+string(clip, format='(i03)')+'_flat'+namestr+'.fits', linehead)
  linemed = median(lineim, dim=3)
  writefits, 'Line_clip'+string(clip, format='(i03)')+'_flat'+namestr+'_med.fits', linemed, linehead
  contim = readfits('Cont_clip'+string(clip, format='(i03)')+'_flat'+namestr+'.fits', conthead)
  contmed = median(contim, dim=3)
  writefits, 'Cont_clip'+string(clip, format='(i03)')+'_flat'+namestr+'_med.fits', contmed, conthead

  if keyword_set(stp) then stop
end