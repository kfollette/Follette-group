pro fkpl_explore, dataset

;  haim_2=readfits('Ha_im_'+dataset+'.fits')
;  subradprofim, haim_2
;  haim_1=readfits('Ha_im_'+dataset+'.fits')
;  cim_2=readfits('Cont_im_'+dataset+'.fits')
;  subradprofim, haim_2
;  cim_1=readfits('Cont_im_'+dataset+'.fits')
;  rotoffs=readfits('rotoff_'+dataset+'.fits')
  haim=readfits('Ha_10coadd_wfe98.5_nocosmics_subradprof_zeroes.fits')
  harotoffs=readfits('Ha_10coadd_wfe98.5_nocosmics_subradprof_zeroes.fits')
  cim=readfits('Cont_10coadd_wfe98.5_nocosmics_subradprof_zeroes.fits')
  crotoffs=readfits('Cont_10coadd_wfe98.5_nocosmics_subradprof_zeroes.fits')

  ;;inject fake planets
  contrast=7.D-3
  separation=0.25 ;arcsec
  platescale=0.00798
  sep=separation/platescale
  
  adi_psfmodel, ha_fakes, (size(haim))[1], (size(haim))[2], harotoffs+90-0.59, [sep, sep, sep, sep], [0, 90, 180, 270], fwhm=5.2
  adi_psfmodel, cont_fakes, (size(cim))[1], (size(cim))[2], crotoffs+90-0.59, [sep, sep, sep, sep], [0, 90, 180, 270], fwhm=5.2
  ;adi_psfmodel, ha_fakes, (size(haim_1))[1], (size(haim_1))[2], rotoffs+90-0.59, [sep, sep, sep, sep], [0, 90, 180, 270], psf0=haim_1
 ; adi_psfmodel, cont_fakes, (size(cim_1))[1], (size(cim_1))[2], rotoffs+90-0.59, [sep, sep, sep, sep], [0, 90, 180, 270], psf0=cim_1
  
  haims_in1=dblarr((size(haim_1))[1], (size(haim_1))[2],(size(haim_1))[3]) 
  haims_in1=haim_1+contrast*ha_fakes
  haims_in2=haim_2+contrast*ha_fakes
  cims_in1=cim_1*1.4  ;;1.4 is cloogy estimate of brightness diff
  cims_in2=cim_2*1.4
  
  ;;grid of radial zone sizes
  rzones=[50]
  azzones=[360]
  rotmask=[1,2,3]
  modes=[5,10,20,50,75,100]
  
    for i=0, n_elements(rzones)-1 do begin ; rad zone loop
      for j=0, n_elements(azzones)-1 do begin
        for k=0, n_elements(rotmask)-1 do begin
          instring='r'+string(rzones[i], format='(i03)')+'_az'+string(azzones[j], format='(i03)')+'_m'+string(rotmask[k], format='(f4.2)')+'_haref'
          instring2='r'+string(rzones[i], format='(i03)')+'_az'+string(azzones[j], format='(i03)')+'_m'+string(rotmask[k], format='(f4.2)')+'_allref'
          pca_regions, finim_haref, haims_in1, rotoffs+50-0.59, rotmask[k], rzones[i], azzones[j], [5,10,20,50,75,100], minrad=5, fitsfile=string(instring)
          pca_regions, finim_allref, haims_in1, rotoffs+50-0.59, rotmask[k], rzones[i], azzones[j], [5,10,20,50,75,100], minrad=5, ref_ims=cims_in1, fitsfile=string(instring2)
        endfor
        ;;no mask for contref only
        instring3='r'+string(rzones[i], format='(i03)')+'_az'+string(azzones[j], format='(i03)')+'_contref'
        pca_regions, finim_allref, haims_in1, rotoffs+50-0.59, 0, rzones[i], azzones[j], [5,10,20,50,75,100], ref_ims=cims_in1, /refonly, fitsfile=string(instring3)
      endfor
    endfor

;;same for radial profile subtracted version
    for i=0, n_elements(rzones)-1 do begin ; rad zone loop
      for j=0, n_elements(azzones)-1 do begin
        for k=0, n_elements(rotmask)-1 do begin
          instring='rps_r'+string(rzones[i], format='(i03)')+'_az'+string(azzones[j], format='(i03)')+'_m'+string(rotmask[k], format='(f4.2)')+'_haref'
          instring2='rps_r'+string(rzones[i], format='(i03)')+'_az'+string(azzones[j], format='(i03)')+'_m'+string(rotmask[k], format='(f4.2)')+'_allref'
;          pca_regions, finim_haref, haims_in2, rotoffs+50-0.59, rotmask[k], rzones[i], azzones[j], [5,10,20,50,75,100], minrad=5, fitsfile='Ha_'+string(instring)
          pca_regions, finim_allref, haims_in2, rotoffs+50-0.59, rotmask[k], rzones[i], azzones[j], modes, minrad=5, ref_ims=cims_in2, fitsfile='Ha_'+string(instring2)
          pca_regions, cfinim_allref, cims_in2, rotoffs+50-0.59, rotmask[k], rzones[i], azzones[j], modes, minrad=5, ref_ims=cims_in2, fitsfile='Cont_'+string(instring2)
          ASDI=finim_allref-cfinim_allref
          writefits, 'ASDI_'+string(instring2)+'.fits', ASDI
        endfor
        ;no mask for contref only
        instring3='rps_r'+string(rzones[i], format='(i03)')+'_az'+string(azzones[j], format='(i03)')+'_contref'
        pca_regions, finim_allref, haims_in2, rotoffs+50-0.59, 0, rzones[i], azzones[j], [5,10,20,50,75,100], ref_ims=cims_in2, /refonly, fitsfile=string(instring3)
      endfor
    endfor
end