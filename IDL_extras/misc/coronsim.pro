pro coronsim, npix, smp, teldiam, wl, aberr=aberr, pupil=pupil, apod=apod, fpm=fpm, lyot=lyot

  ;+
  ; PURPOSE:
  ;  Coronagraphic image simulator
  ;
  ; INPUTS:
  ;  npix     = number of pixels in simulated image
  ;  smp      = pixel sampling in multiples of Nyquist
  ;  teldiam  = telescope diameter (in meters)
  ;  wl       = output science wavelength (in meters)
  ;
  ; INPUT KEYWORDS:
  ; aberr     = phase map as a .fits file OR the stddev of a randomly generated aberrated phase map in units of radians
  ; pupil     = pupil transmission function as a .fits file
  ; fpm       = size of focal plane mask in arcseconds
  ; apod      = apodizer transmission function as a .fits file
  ; lyot      = lyot mask as a .fits file
  ;
  ; OUTPUTS:
  ;
  ;
  ; HISTORY:
  ;  Written 2015-05-29 by Kate Follette
  ;  Revised 2015-08-06 to include Fourier resampling
  ;-

  ;calculate pixel scale for array to be transformed
  pixscale = teldiam/(npix/2.)*smp
  print, 'pixel scale is   ', pixscale*100, '  centimeters per subaperture'
  teldiampix = round(teldiam/pixscale)
  print, teldiampix, '  pixels across telescope pupil'
  md=(teldiampix-1)/2.
  st=(npix-1)/2.

  ;; insert telescope pupil of appropriate size into array to be transformed
  ;if fits image provided, read in and resize
  if keyword_set(pupil) then begin
    pup=readfits(pupil)
    cen=((size(pup))[1]-1.)/2.
    if keyword_set(apod) then begin
      apodiz=readfits(apod)
      pup=pup*apodiz
    endif

    ;v1 resample pupil via interpolation - leaves lots of artifacts
    ;        pupmap=ft(pup)
    ;        puprsz=ft(pupmap[cen-st:cen+st,cen-st:cen+st], /inverse)
    ;        puprsz=congrid(pup, teldiampix, teldiampix, cubic=-0.5)

    ;v2
    ;;fourier downsample pupil image
    pupft=ft(pup)
    pupft_clip=pupft[cen-md:cen+md,cen-md:cen+md]
    puprsz=abs(ft(pupft_clip, /inverse))
    pupil=dblarr(npix, npix)
    pupil[st-md:st+md,st-md:st+md]=puprsz  ;;THIS SEEMS TO BE THE STEP WHERE WE GET OFF CENTER BY 1 PIXEL
    scl2=stdev(pup)/stdev(pupil)
    pupil=pupil*scl2
  endif else begin    ;otherwise assume filled circular aperture
    pupil=dblarr(npix, npix)
    pupil(where(radmap(pupil) le teldiampix/2.0))=1.0
  endelse

  ;; calculate pixel scale in image plane
  imscale=wl/(npix*pixscale)*206000
  print, 'pixel scale is   ', imscale, '  arcseconds per pixel in image plane'

  ;;create focal plane mask transmission function (0 inside FPM, 1 elsewhere)
  foc=dblarr(npix,npix)+1.
  if keyword_set(fpm) then foc(where(radmap(foc) le fpm/imscale))=0.0

  ;resize lyot mask if present
  if keyword_set(lyot) then begin
    lytraw=readfits(lyot)
    lyt=dblarr(npix, npix)

    ;;Fourier downsample lyot mask to match pupil
    lytft=ft(lytraw)
    cen=((size(lytraw))[1]-1.)/2.
    lytft_clip=lytft[cen-md:cen+md,cen-md:cen+md]
    lytrsz=abs(ft(lytft_clip, /inverse))
    lyt[st-md:st+md,st-md:st+md]=lytrsz
  endif

  ;create unaberrated PSF
  e_pupil_unaberr = pupil
  e_focal_unaberr = ft(e_pupil_unaberr)
  inten_focal_unaberr = abs(e_focal_unaberr)^2
  e_fpm_unaberr=e_focal_unaberr*foc
  inten_fpm_unaberr=abs(e_fpm_unaberr)^2
  if keyword_set(lyot) then begin
    e_lyot_unaberr=ft(e_fpm_unaberr, /inverse)
    e_final_unaberr=ft(e_lyot_unaberr*lyt)
    inten_final_unaberr=abs(e_final_unaberr)^2
  endif else begin
    inten_final_unaberr=inten_fpm_unaberr
  endelse

  ;;set up phase map

  phase=dblarr(npix,npix)

  ;;if input phase map is a fits image, read in and set up binning if 3D
  if keyword_set(aberr) and type(aberr) eq 7 then begin
    phaseraw=readfits(aberr)*5 ;;Remove this *5 later!!

    ;;if 3D
    if (size(phaseraw))[0] eq 3 then begin
      ndim3=(size(phaseraw))[3]
      nbin=ceil(ndim3/100)
      lastbinsz=ndim3-(100*(nbin-1))
      print, 'input phase map is 3D and will be processed in', nbin, ' image bins'
      ;;if 2D
    endif else begin
      nbin=1
      binsz=1
      ndim3=1
    endelse

  endif else begin
    ;;no phase aberration option
    if not keyword_set(aberr) then aberr=0.
    nbin=1
    binsz=1
  
  endelse
  
  phaseft_pad=dblarr(teldiampix, teldiampix)
  if nbin eq 1 then binsz=1 else binsz=100
  inten_final=dblarr(npix,npix,binsz)
  inten_final_meds=dblarr(npix,npix,nbin)
  phsz=47/2.

  for i=0, nbin-1 do begin
    for j=0, binsz-1 do begin
      if i eq nbin-1 and i ne 0 then binsz=lastbinsz else binsz=100
      if nbin eq 1 then binsz=1
      if keyword_set(aberr) and type(aberr) eq 7 then begin
        ;phase maps are always 48x48, but need to map to pupil
        ;;forier upsample pupil to appropriate size
        print, 'processing bin', i+1, ' image number', j+1, ' which is', i*100+j+1, ' of', ndim3, ' total'
        phaseft=ft(phaseraw[*,*,i*100+j])
        phaseft_pad[md-phsz:md+phsz,md-phsz:md+phsz]=phaseft
        phasersz=abs(ft(phaseft_pad, /inverse))
        phase[st-md:st+md,st-md:st+md]=phasersz/(wl*1.d6)*2*!dpi
        scl=stdev(phaseraw)/stdev(phase)
        phase=phase*scl
        ;stop
      endif else begin
        ;;random phase option
        phaseraw=randomn(42, 48,48)*aberr
        phaseft=ft(phaseraw)
        phaseft_pad[md-phsz:md+phsz,md-phsz:md+phsz]=phaseft
        phasersz=abs(ft(phaseft_pad, /inverse))
        phase[st-md:st+md,st-md:st+md]=phasersz
        scl=stdev(phaseraw)/stdev(phase)
        phase=phase*scl
      endelse

      ;Now comes the actual optics propagation

      ;electric field in pupil plane
      im=complex(0,1)
      e_pupil= pupil* exp(im*phase)
      
      ;stop
      ;electric field in image plane
      e_focal = ft(e_pupil)

      ;electric field intensity in image plane (~PSF)
      inten_focal=abs(e_focal)^2

      ;apply focal plane mask in image plane
      e_fpm=e_focal*foc
      inten_fpm=abs(e_fpm)^2

      if keyword_set(lyot) then begin
        ;apply Lyot stop in pupil plane
        e_lyot=ft(e_fpm, /inverse)
        e_final=ft(e_lyot*lyt)
        inten_final[*,*,j]=abs(e_final)^2
      endif else begin
        inten_final[*,*,j]=inten_fpm
      endelse
    endfor
    if nbin gt 1 then inten_final_meds[*,*,i]=total(inten_final, 3)
  endfor

  ;;create output plots
  if nbin eq 1 then final_image=inten_final else final_image=total(inten_final_meds, 3)

  cen=(npix-1.)/2

  if keyword_set(lyot) then !P.MULTI=[0,2,3] else !P.MULTI=[0,2,2]
  loadct, 1
  rng=teldiampix/2
  cgimage, pupil[cen-rng:cen+rng,cen-rng:cen+rng]*255
  cgtext, 0.5, 0.5, 'Pupil', color='white'
  cgimage, alog10(phase[cen-rng:cen+rng,cen-rng:cen+rng]);*255
  cgtext, 0.5, 0.5, 'Phase', color='white'
  FOV=2.7
  FOVas=FOV/imscale/2.
  cgimage, alog10(inten_final_unaberr[cen-FOVas:cen+FOVas,cen-FOVas:cen+FOVas]);*abs((alog10(max(inten_final_unaberr)-min(inten_final_unaberr)))/255)+alog10(min(inten_final_unaberr))
  cgtext, 0.5, 0.5, 'Unaberrated PSF', color='white'
  cgimage, alog10(inten_fpm[cen-FOVas:cen+FOVas,cen-FOVas:cen+FOVas]);*abs((alog10(max(inten_focal)-min(inten_focal)))/255)+alog10(min(inten_focal))
  cgtext, 0.5, 0.5, 'Focal Plane with FPM', color='white'
  ;  if keyword_set(lyot) then cgimage, lyt[cen-rng:cen+rng,cen-rng:cen+rng]*255
  if keyword_set(lyot) then cgimage, alog10(abs(e_lyot[cen-rng:cen+rng,cen-rng:cen+rng]))
  if keyword_set(lyot) then cgtext, 0.5, 0.5, 'Lyot Plane Before Stop', color='white'
  if keyword_set(lyot) then cgimage, alog10(final_image[cen-FOVas:cen+FOVas,cen-FOVas:cen+FOVas]);*abs((alog10(max(inten_focal)-min(inten_focal)))/255)+alog10(min(inten_focal))
  if keyword_set(lyot) then cgtext, 0.5, 0.5, 'Final Focal Plane Image with Lyot', color='white'
  stop
end