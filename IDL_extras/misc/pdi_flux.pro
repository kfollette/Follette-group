pro pdi_flux, objname, expt, imtype, boxsz, stellar_flux

  ;;flux conversion

  
  ;;;Flux Calibration and scaling
  Unsat_Imed=readfits(string(objname)+'_'+strcompress(string(expt),/remove_all)+ $
    '_Imed.fits', junk, /NO_UNSIGNED)
    ;; sky subtract
  Unsat_sub=Unsat_Imed-median(Unsat_Imed)
  ndim1=(size(Unsat_Imed))[1]
  ndim2=(size(Unsat_Imed))[2]
  xcen=(ndim1-1)/2.
  ycen=(ndim2-1)/2.
  total_flux=total(Unsat_sub[xcen-boxsz:xcen+boxsz,ycen-boxsz:ycen+boxsz])
  unsat_conv=stellar_flux/(total_flux)
  print, 'unsat total', total_flux, 'unsat scale', unsat_conv
  
  Image=readfits(string(objname)+'_'+strcompress(string(expt),/remove_all)+'_'+string(imtype)+'.fits', junk, /NO_UNSIGNED)
  Image_bgd=Image-median(Image)
  
  ADU_to_Jy=(unsat_conv)/(.0095)^2*1000 ;;mJy per sq arcsec
  print, 'conversion', ADU_to_Jy
  
  writefits, string(objname)+'_'+strcompress(string(expt),/remove_all)+'_'+string(imtype)+'_flux.fits', Image_bgd*ADU_to_Jy
 
  stop
  
end
