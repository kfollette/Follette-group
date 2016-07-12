  ;+
  ; NAME:
  ;   posplanets
  ;
  ; DESCRIPTION:
  ;   Injects a positive planet at given contrast, separation and PA or grid of these values. 
  ;
  ; INPUTS:
  ;   ims        : a cube of real images
  ;   rotoffs    : an array of rotoff angles
  ;   contrasts  : an array of contrast values
  ;   PAs        : an array of position angles, in degrees
  ;   seps       : an array of separation values, in pixels
  ;
  ; INPUT KEYWORDS:
  ;   pim        : a separate image array to inject in cases where image is saturated or radial profile has been subtracted from image cube
  ;   cont       : continuum image to subtract from final data product to create asdi image
  ;
  ; KEYWORDS:
  ;
  ; OUTPUTS:
  ;
  ; OUTPUT_KEYWORDS:
  ;    none
  ;
  ; MODIFICATION HISTORY:
  ;  Written 2015/04/07 by Katherine B. Follette (kfollette@stanford.edu)
  ;
  ;-
  
  pro posplanets, ims, rotoffs, contrasts, PAs, seps, typ, mr=mr, msk=msk, rz=rz, azz=azz, subdir=subdir, pim=pim, immask=immask

  if keyword_set(mr) then mr=mr else mr=7
  if keyword_set(msk) then msk=msk else msk=1
  if keyword_set(rz) then rz=rz else rz=ceil((size(ims))[1]/2.)
  if keyword_set(azz) then azz=azz else azz=360
  if keyword_set(subdir) then subdir=subdir else subdir='.'
  if keyword_set(pim) then pim=pim else pim=ims
  if not keyword_set(immask) then immask=ims*0.+1
  x=file_test('./'+string(subdir), /DIRECTORY)
  if x eq 0 then begin
    spawn, 'mkdir '+string(subdir)
  endif

  ;;create Halpha map with no planet emission for ASDI
;  haims=readfits('imsnop_final.fits')
;  harotoffs=readfits('Ha_rotoff_cube_wfe170pt5_nocosmics.fits')      
;  pca_regions, hafinim, haims, harotoffs+90-0.59, msk, rz, azz, [2,4,6,8,10,20,50,100], minrad=mr, sigma=3

cfinim=readfits('./KLIP_channels/KLIP_Cont_radprofsub_pt25mask_rad50_az360_minrad3_3sigclip_10modes.fits')

  for j=0, n_elements(PAs)-1 do begin
    for k=0, n_elements(seps)-1 do begin
      PA=PAs[j]
      sep=seps[k]
      adi_psfmodel, model, (size(pim))[1], (size(pim))[2], rotoffs+90, sep, PA, psf0=pim
      for i=0, n_elements(contrasts)-1 do begin
        contrast=contrasts[i]
        fkpl=ims+(model*contrast)
        pca_regions, finim, fkpl*immask, rotoffs+90-0.59, msk, rz, azz, [2,4,6,8,10,20,50,100], minrad=mr, sigma=3, mask=immask
        writefits, string(subdir)+'/'+string(typ)+'_POSPL_PCA_rad'+string(rz,format='(i03)')+'_'+string(azz,format='(i03)')+'_mask'+$
          string(msk,format='(i1)')+'_minrad'+string(mr,format='(i02)')+'_sep'+string(sep, format='(i03)')+'pix_PA'+$
          string(PA,format='(i03)')+'_contrast'+string(contrast,format='(E7.1)')+'.fits', finim
          ASDI=finim[*,*,4]-cfinim
        smpl=filter_image(ASDI, fwhm=2, /all)
        snrmap, smpl, fkplmap, mask=[136,123,5]
        writefits, string(subdir)+'/'+string(typ)+'_POSPL_PCA_rad'+string(rz,format='(i03)')+'_'+string(azz,format='(i03)')+'_mask'+$
          string(msk,format='(i1)')+'_minrad'+string(mr,format='(i02)')+'_sep'+string(sep, format='(i03)')+'pix_PA'+$
          string(PA,format='(i03)')+'_contrast'+string(contrast,format='(E7.1)')+'_10modes_02pixsm_snrmap.fits', fkplmap     
      endfor
    endfor
  endfor

end
