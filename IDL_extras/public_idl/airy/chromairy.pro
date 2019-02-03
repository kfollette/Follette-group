;+
; NAME: chromairy
;
; PURPOSE: Generate an ideal PSF given a pupil, bandpass, etc.
;
; DESCRIPTION:
;  Calculates a perfect PSF in a bandpass, can include effects of atmosphere and stellar blackbody.
;  If the fitsname keyword is set, writes a fits file, with a detailed header.
;
; INPUTS:
;  telD       ;  telescope diameter (meters)
;  eps        ;  central obscuration fraction  (0<=eps<1)
;  imsz       ;  linear size of the output image, e.g. 64 for a 64x64 images
;  ps         ;  platescale of the output image (arcsec/pixel)
;  passlam    ;  vector containing the filter curve wavelength scale (microns)
;  passtrans  ;  vector containing the filture curve transmission (0<= passtrans < 1)
;
; KEYWORDS:
;  teff       ;  the temperature of the stellar spectrum, in Kelvin
;  atmo       ;  apply atmospheric transmission
;  bpname     ;  name of the bandpass, only used in the fits header
;  fitsname   ;  string name of the fits file output
;  inairy     :  2D image of an airy pattern to use instead of the analytic one
;  inps       :  The platescale of inairy
;  inlam      :  The input wavelength of inairy (microns)
;
; OUTPUTS:
;  ac         ;  the PSF (ac = airy-chromatic)
;  meanlam    ;  the transmission weighted mean wavelength of the bandpass
;  pix        ;  array containing the radial position of each pixel, in arcsec
;  bstar      ;  the normalized planck function for teff
;  atmot      ;  the atmospheric transmission
;  tottrans   ;  the resultant total transmission
;
; EXAMPLE:
;  chromairy, 6.5, 0.29, 512, 0.008, iplam, iptrans, ac
;
; HISTORY:
;  Written 2012-10-30 by Jared Males, jrmales@email.arizona.edu
;  Updated 2012-11-06 Added fits file output. (Jared Males)
;  Updated 2013-01-28 Documentation clarified (Katie Morzinski, ktmorz@arizona.edu)
;
;-
pro chromairy, telD, eps, imsz, ps, passlam, passtrans, ac, meanlam, pix, bstar, atmot, tottrans, teff=teff, atmo=atmo, bpname=bpname, fitsname=fitsname, $
inairy=inairy, inps=inps, inlam=inlam

;Set flag on whether to use analytic airy, or interpolate the input airy pattern
ina = 0 
if(keyword_set(inairy) and keyword_set(inps) and keyword_set(inlam)) then begin
   ina = 1
endif

;Generate stellar spectrum
if(keyword_set(teff)) then begin
   ;I think these have to be right for this to work.
   h = 6.626068d-34 ;planck's constant
   c = 299792458.d; speed of light
   k = 1.3806503d-23

   bstar = 2.*h*c^2/(passlam*1d-6)^5/(exp(h*c/((passlam*1d-6)*k*teff) - 1.0))
   
   bstar = bstar/max(bstar) ;just normalize it
   
endif else begin

   bstar = dblarr(n_elements(passlam)) + 1.d
   
endelse

;Generate atmospheric transmission.  Code provided by Katie Morzinski.
if(keyword_set(atmo)) then begin
   ;atmosphere transmission (throughput) from Allen's Astrophysical Quantities p.264
   atmos_l = [10000., 7500, 5000, 4000, 3000, 2000, 1000, 900, 800, 700, 650, 600, 550, 500, 450, 400, 380, 360, 340, 320, 300, 280, 260, 230, 220, 200]
   atmos_t = [.759, .025, .294, .837, .376, .441, .846, .645, .767, .709, .699, .659, .637, .591, .527, .438, .394, .345, .287, .177, .0062, 0, 0, 0, 0, 0]

   linterp, atmos_l/1.d3, atmos_t, passlam, atmotf
   
   atmot = double(atmotf)
   
endif else begin

   atmot = dblarr(n_elements(passlam)) + 1.d
   
endelse
   
pix = dblarr(imsz,imsz)

cen = .5d*double(imsz)

for i=0.d,imsz-1.d do begin
   for j=0.d,imsz-1.d do begin
      pix[i,j] = sqrt((i-cen)^2. + (j-cen)^2.)*ps
   endfor
endfor

ac = dblarr(imsz,imsz)

tottrans = passtrans*bstar*atmot

for i=0, n_elements(passlam)-1 do begin

   lamd = pix/(0.2063d*double(passlam[i])/double(telD))

   if(ina eq 0) then begin
      ;Use the analytic airy pattern
      a = airy(lamd, eps)
   endif else begin
      ;Now re-sample the input airy pattern
      a = rot(inairy, 0., (inps/ps)*(passlam[i]/inlam), cubic=-.5)
   endelse
   
   ac = ac + a*tottrans[i]

endfor

meanlam = total(passlam*tottrans)/total(tottrans)

;Generate fits header and write fits file
if(keyword_set(fitsname)) then begin

   mkhdr, header, ac, /EXTEND

   headend = header[n_elements(header)-1-3:*]

   header[n_elements(header)-1-3] = 'COMMENT Chromatic Airy Pattern made with chromairy.pro'

   if(keyword_set(bpname)) then begin
      header[n_elements(header)-1-2] = 'BANDPASS=    ' + bpname
   endif else begin
     header[n_elements(header)-1-2] = 'BANDPASS=       unknown'
   endelse
   
   header[n_elements(header)-1-1] = 'TELDIAM =           ' + string(float(telD), format='(F10.2)') + ' / Telescope diameter (m)' 

   header[n_elements(header)-1-0] = 'CENOBS  =           ' + string(float(eps), format='(F10.3)') + ' / Central obscuration'
   
   header =  [header, 'PLTSCL  =           ' + string(float(ps), format='(F10.5)') + ' / Plate scale (arcsec/pix)']
   header = [header, 'MEANLAM =           ' + string(float(meanlam), format='(F10.5)') + ' / Mean wavelength (um)']
   if(keyword_set(teff)) then begin
      header = [header, 'TEFF    =           ' + string(float(teff), format='(F10.1)') + ' / Stellar effective temp (K)']
   endif else begin
      header = [header, 'TEFF    =                 none / Stellar effective temp (K)']
   endelse
   if(keyword_set(atmo)) then begin
      header = [header, 'ATMO    =                 yes  / atmosphere included']
   endif else begin
      header = [header, 'ATMO    =                  no  / atmosphere not included']
   endelse

   header = [header,  'ORIGFILE=    ' + fitsname + ' / original file name.']
   
   header = [header, headend]

   writefits, fitsname, ac, header

endif

end



