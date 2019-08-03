;+
; NAME: visao_combine_registered_ims
;
; PURPOSE:
;  combine registered images with different exposure times into one image cube
;
; INPUTS:
;  dir_list
;
; INPUT KEYWORDS:
;
; OUTPUTS:
;
; OUTPUT KEYWORDS:
;    none
;
; EXAMPLE:
;  visao_combine_registered_ims, ['merged_3sec/','merged_5sec/'], clip=451, /flat
;
; HISTORY:
;  Written 2019-08-03 by Kate Follette
;-

pro visao_combine_registered_ims, dir_list, clip=clip, flat=flat, stp=stp

  ;;naming specifications  
  if keyword_set(clip) then namestr = '_clip'+string(clip, format='(i03)') else namestr=''
  if keyword_set(flat) then namestr=string(namestr)+'_flat_' else namestr='_'
  
  nexp = n_elements(dir_list)
  nims=dblarr(nexp)
  expts=dblarr(nexp)
  total_ims=0
  
  for i=0, nexp-1 do begin
    dummyim = readfits(string(dir_list[i])+'Line'+namestr+'reg.fits')
    exp = readfits(string(dir_list[i])+'exptime_preproc.fits')
    nims[i]=(size(dummyim))[3]
    expts[i]=exp[0]
    print, 'found', nims[i], ' images in directory ', dir_list[i], ' with exposure time', expts[i]   
    if i gt 0 then total_ims=nims[i-1]+nims[i]
    print, 'total images', total_ims
    delvar, dummyim
  endfor
  
  ;make an array to store everything
  new_Line_arr = dblarr(clip, clip, total_ims)
  new_Cont_arr = dblarr(clip, clip, total_ims)
  rotoffs = dblarr(total_ims)
  last=0
  
  for i=0, nexp-1 do begin
    Lineim = readfits(string(dir_list[i])+'Line'+namestr+'reg.fits', Linehead)
    Contim = readfits(string(dir_list[i])+'Cont'+namestr+'reg.fits', Conthead)
    rots = readfits(string(dir_list[i])+'rotoff_preproc.fits')
    new_Line_arr[*,*,last:last+nims[i]-1]=Lineim/expts[i]
    new_Cont_arr[*,*,last:last+nims[i]-1]=Contim/expts[i]
    rotoffs[last:last+nims[i]-1]=rots
    last=last+nims[i]
  endfor
  
sxaddpar, Linehead, 'EXPT', 'normalized to 1'
sxaddpar, Conthead, 'EXPT', 'normalized to 1'

writefits, 'Line'+namestr+'reg_combo.fits', new_Line_arr, Linehead
writefits, 'Cont'+namestr+'reg_combo.fits', new_Cont_arr, Conthead
writefits, 'rotoffs_preproc.fits', rotoffs

if keyword_set(stp) then stop

end

