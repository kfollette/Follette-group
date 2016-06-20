pro visao_removecosmics, sdi, clip=clip

;;read in sdi image 
SDIim=readfits('SDI_sc'+string(sdi, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg_circsym.fits')

ds9_imselect, SDIim, index=idx

writefits, 'SDI_sc'+string(sdi, format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg_circsym_nocosmics.fits', SDIim
writefits, 'cosmic_arr.fits', idx

rotoffs=readfits('rotoff_preproc.fits')
wfe=readfits('avgwfe_preproc.fits')
exptime=readfits('exptime_preproc.fits')

Line=readfits('Line_clip'+string(clip, format='(i03)')+'_reg_circsym.fits')
Cont=readfits('Cont_clip'+string(clip, format='(i03)')+'_reg_circsym.fits')
SDI2=readfits('SDI_sc'+string(1., format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg_circsym.fits')

Line_noc=dblarr(clip,clip,n_elements(idx))
Cont_noc=dblarr(clip,clip,n_elements(idx))
SDI2_noc=dblarr(clip,clip,n_elements(idx))
rotoffs_noc=dblarr(n_elements(idx))
wfe_noc=dblarr(n_elements(idx))
exptime_noc=dblarr(n_elements(idx))


for i=0, n_elements(idx)-1 do begin
  print, i
  Line_noc[*,*,i]=Line[*,*,idx[i]]
  Cont_noc[*,*,i]=Cont[*,*,idx[i]]
  SDI2_noc[*,*,i]=SDI2[*,*,idx[i]]
  rotoffs_noc[i]=rotoffs[idx[i]]
  wfe_noc[i]=wfe[idx[i]]
  exptime_noc[i]=exptime[idx[i]]
endfor

writefits, 'Line_clip'+string(clip, format='(i03)')+'_reg_circsym_nocosmics.fits', Line_noc
writefits, 'Cont_clip'+string(clip, format='(i03)')+'_reg_circsym_nocosmics.fits', Cont_noc
writefits, 'SDI_sc'+string(1., format='(f05.2)')+'_clip'+string(clip, format='(i03)')+'_reg_circsym_nocosmics.fits', SDI2_noc
writefits, 'rotoff_nocosmics.fits', rotoffs_noc
writefits, 'avgwfe_nocosmics.fits', wfe_noc
writefits, 'exptime_nocosmics.fits', exptime_noc

end
