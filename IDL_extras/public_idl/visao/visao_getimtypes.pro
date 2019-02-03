;+
; NAME: visao_getimtypes
;
; PURPOSE:
;  Extracts image type and other headers from all the VisAO images in a directory
;
; DESCRIPTION:
;   Generates a list of files in a directory, which can be specified with the subdir keyword, which match the 
;   regex "V47_*.fits".  The prefix can be changes with the prefix keyword (e.g. to 'dsub_').  At minimum returns
;   the filenames and the image types of the files (science, dark, etc.).  If the ims keyword is given an argument, 
;   the images are returned as a cube, which by default is fltarr. To avoid re-generating the list of files, and instead use
;   the file names passed in as fnames, set the keyword usefnames.  Various header keyword values can also be extracted.
;
; INPUTS:
;  none
;
; INPUT KEYWORDS:
;   prefix      :  the filename prefix, default is 'V47_'
;   subdir      :  if set, operates on images in the specified directory, otherwise uses the current directory
;   region      :  a 4 element vector specifying an area of the images to return. [x0,x1, y0,y1]
;   goodims     :  a vector if indices for the images you want to retrieve, if not set all images are retrieved
;   usefnames   : load the images specified in fnames, rather than getting all images in subdir
;   double      : load images as doubles.  default is float.
;
; OUTPUTS:
;   fnames      :   The file names of the images
;   imtypes     :   Image types, 0=sci, 1=acq, 2=dark, 3=sky, 4=flat
;
; OUTPUT KEYWORDS:
;   ims         : a cube of images with format [size, size, no_images]
;   [KEYWORD]   : a vector containing the value of the keyword for each image, e.g. AOLOOPST gets the loop status.
;                 refer to a visao fits header for all available keywords.
;
; EXAMPLE:
;  visao_getimtypes, fnames, imtypes, exptime=exptime, ims=ims
;
; HISTORY:
;  Written 2012-11-14 by Jared Males, jrmales@email.arizona.edu
;          2013-09-01 documentation update by JRM.
;
; BUGS/WISH LIST:
;  Supported keywords are incomplete, need to expand.
;
;-
pro visao_getimtypes, fnames, imtypes, EXPTIME=EXPTIME, VFOCPOS=VFOCPOS, AOLOOPST=AOLOOPST, ROTOFF=ROTOFF, $
                      VGIMXPOS=VGIMXPOS, VGIMYPOS=VGIMYPOS, AVGWFE=AVGWFE, STDWFE=STDWFE, VFW1POSN=VFW1POSN, $
                      VFW2POSN=VFW2POSN, VFW3POSN =VFW3POSN, DATEOBS=DATEOBS, AOREC=AOREC, UT=UT, AM=AM, $
                      HA=HA,FRAMENO=FRAMENO, FGDMATS=FGDMATS, ORIGXCEN=ORIGXCEN, ORIGYCEN=ORIGYCEN, GOODCENT=GOODCENT, $
                      ROTANG=ROTANG, GAIN=GAIN, OBJECT=OBJECT, MAG1FWHM=MAG1FWHM, DIMMFWHM=DIMMFWHM, pixrt=pixrt, $
                      ims=ims, prefix=prefix, subdir=subdir, region=region, goodims=goodims, usefnames=usefnames, $
                      double=double

;** First thing we do is generate a list of fits files **

if(n_elements(subdir) eq 0) then begin
   subdir = "."
endif

if(~keyword_set(usefnames) or n_elements(fnames) eq 0) then begin 
   usefnames = 0
   if(n_elements(prefix) eq 0) then prefix='V47_'

   srchstr = strcompress(subdir + '/' + prefix+'*.fits', /remove)
   
   fnames = file_search(srchstr)
   
endif else begin
   usefnames = 1
endelse


if(n_elements(goodims) gt 0) then begin
   fnames = fnames[goodims]
endif 
  
;*** allocate outputs ***

nims = n_elements(fnames)
   

imtypes = intarr(nims)

;*** allocate optional outputs ***
get_aoloopst = 0
if( arg_present(AOLOOPST) ) then begin
   get_aoloopst = 1
   aoloopst = intarr(nims)
endif

get_exptime = 0
if( arg_present(EXPTIME) ) then begin
   get_exptime = 1
   exptime = dblarr(nims)
endif

get_vfocpos = 0
if( arg_present(VFOCPOS) ) then begin
   get_vfocpos = 1
   vfocpos = dblarr(nims)
endif

get_rotoff = 0
if( arg_present(ROTOFF) ) then begin
   get_rotoff = 1
   rotoff = dblarr(nims)
endif

get_rotang = 0
if( arg_present(ROTANG) ) then begin
   get_rotang = 1
   rotang = dblarr(nims)
endif

get_vgimxpos = 0
if( arg_present(vgimxpos) ) then begin
   get_vgimxpos = 1
   vgimxpos = dblarr(nims)
endif

get_vgimypos = 0
if( arg_present(vgimypos) ) then begin
   get_vgimypos = 1
   vgimypos = dblarr(nims)
endif

get_avgwfe = 0
if( arg_present(avgwfe) ) then begin
   get_avgwfe = 1
   avgwfe = dblarr(nims)
endif

get_stdwfe = 0
if( arg_present(stdwfe) ) then begin
   get_stdwfe = 1
   stdwfe = dblarr(nims)
endif

get_vfw1posn = 0 
if( arg_present(vfw1posn) ) then begin
   get_vfw1posn = 1
   vfw1posn = strarr(nims)
endif

get_vfw2posn = 0
if( arg_present(vfw2posn) ) then begin
   get_vfw2posn = 1
   vfw2posn = strarr(nims)
endif

get_vfw3posn = 0
if( arg_present(vfw3posn) ) then begin
   get_vfw3posn = 1
   vfw3posn = strarr(nims)
endif

get_dateobs = 0
if( arg_present(dateobs) ) then begin
   get_dateobs = 1
   dateobs = dblarr(nims)
endif

get_fgdmats = 0
if( arg_present(fgdmats) ) then begin
   get_fgdmats = 1
   fgdmats = dblarr(nims)
endif

get_ut = 0
if( arg_present(ut) ) then begin
   get_ut = 1
   ut = dblarr(nims)
endif

get_am = 0
if( arg_present(am) ) then begin
   get_am = 1
   am = dblarr(nims)
endif

get_ha = 0
if( arg_present(ha) ) then begin
   get_ha = 1
   ha = dblarr(nims)
endif

get_aorec=0
if( arg_present(aorec) ) then begin
   get_aorec = 1
   aorec = strarr(nims)
endif

get_gain = 0
if( arg_present(gain) ) then begin
   get_gain = 1
   gain = intarr(nims)
endif

get_frameno=0
if( arg_present(frameno) ) then begin
   get_frameno = 1
   frameno = lonarr(nims)
endif

get_origxcen=0
if( arg_present(origxcen) ) then begin
   get_origxcen = 1
   origxcen = fltarr(nims)
endif

get_origycen=0
if( arg_present(origycen) ) then begin
   get_origycen = 1
   origycen = fltarr(nims)
endif

get_goodcent=0
if( arg_present(goodcent) ) then begin
   get_goodcent = 1
   goodcent = intarr(nims)
endif

get_object = 0
if( arg_present(object) ) then begin
   get_object = 1
   object = strarr(nims)
endif

get_dimmfwhm = 0
if( arg_present(dimmfwhm) ) then begin
   get_dimmfwhm = 1
   dimmfwhm = fltarr(nims)
endif

get_mag1fwhm = 0
if( arg_present(mag1fwhm) ) then begin
   get_mag1fwhm = 1
   mag1fwhm = fltarr(nims)
endif

get_pixrt = 0
if( arg_present(pixrt) ) then begin
   get_pixrt = 1
   pixrt = intarr(nims)
endif


;*** Allocate output image if desired ***
retims = 0
if( arg_present(ims) ) then begin
   retims = 1
   
   imraw= mrdfits(fnames[0], 0, /silent)
   ndim1 = (size(imraw))[1]
   ndim2 = (size(imraw))[2]
   
   if(~keyword_set(region)) then begin
      region = [0, ndim1-1, 0, ndim2-1]
   endif
   
   rgx = region
   
   if(keyword_Set(double)) then begin
      ims = dblarr(rgx[1]-rgx[0]+1, rgx[3]-rgx[2]+1, nims)
   endif else begin
      ims = fltarr(rgx[1]-rgx[0]+1, rgx[3]-rgx[2]+1, nims)
   endelse
endif



;*** Now read in each image and header ***
for h=0l, nims-1 do begin ;n_elements(corrected_file)-1 do begin
   status = strcompress(string(h+1) + '/' + string(nims), /rem)
   statusline, status, 0
   
   if(retims eq 1) then begin
      imraw= mrdfits(fnames[h], 0, header, /silent)
      ims[*,*,h] = imraw[rgx[0]:rgx[1],rgx[2]:rgx[3]]
   endif else begin
      fits_read, fnames[h], 0, header, /header_only
   endelse
   
   hidx = strmatch(header,'VIMTYPE*')
   if((where(hidx gt 0))[0] ne -1) then begin
      imtstr = (strmid(header[where(hidx gt 0)], 10, 21))
      if( (strmatch(imtstr, '*SCIENCE*'))[0] gt  0) then imtypes[h] = 0
      if( (strmatch(imtstr, '*ACQUISITION*'))[0] gt 0) then imtypes[h] = 1
      if( (strmatch(imtstr, '*DARK*'))[0] gt 0) then imtypes[h] = 2
      if( (strmatch(imtstr, '*SKY*'))[0] gt 0) then imtypes[h] = 3
      if( (strmatch(imtstr, '*FLAT*'))[0] gt 0) then imtypes[h] = 4
   endif

   if( get_aoloopst) then begin
      if(strcompress(sxpar(header,'AOLOOPST'),/remove_all) eq 'CLOSED') then begin
         aoloopst[h] = 1
      endif else begin
         aoloopst[h] = 0
      endelse
   endif
   
   if( get_exptime ) then begin
      exptime[h] = sxpar(header, 'EXPTIME')
   endif
   
   if( get_vfocpos ) then begin
      vfocpos[h] = sxpar(header, 'VFOCPOS')
   endif
  
   if( get_rotoff ) then begin
      rotoff[h] = sxpar(header, 'ROTOFF')
   endif
   
   if( get_rotang ) then begin
      rotang[h] = sxpar(header, 'ROTANG')
   endif
   
   if( get_vgimxpos ) then begin
      vgimxpos[h] = sxpar(header, 'VGIMXPOS')
   endif
   
   if( get_vgimypos ) then begin
      vgimypos[h] = sxpar(header, 'VGIMYPOS')
   endif
   
   if( get_avgwfe ) then begin
      avgwfe[h] = sxpar(header, 'AVGWFE')
   endif
   
   if( get_stdwfe ) then begin
      stdwfe[h] = sxpar(header, 'STDWFE')
   endif
   
   if( get_vfw1posn ) then begin
      vfw1posn[h] = strcompress(sxpar(header, 'VFW1POSN'), /rem)
   endif
   
   if( get_vfw2posn ) then begin
      vfw2posn[h] = strcompress(sxpar(header, 'VFW2POSN'), /rem)
   endif
   
   if( get_vfw3posn ) then begin
      vfw3posn[h] = strcompress(sxpar(header, 'VFW3POSN'), /rem)
   endif
   
   if( get_dateobs ) then begin
      dobs = sxpar(header, 'DATE-OBS')
      dateobs[h] = fitsdate_conv(dobs, 'JULIAN')
   endif
   
   if( get_fgdmats ) then begin
      ts = sxpar(header, 'FGDMATS')
      fgdmats[h] = fitsdate_conv(ts)
   endif
   
   if( get_ut ) then begin
      ut[h] = sxpar(header, 'UT')
   endif
   
   if( get_am ) then begin
      am[h] = sxpar(header, 'AM')
   endif
   
   if( get_ha ) then begin
      ha[h] = sxpar(header, 'HA')
   endif
   
   if( get_aorec ) then begin
      aorec[h] = sxpar(header, 'AOREC')
   endif
   

   if( get_gain ) then begin
      g = strcompress(fxpar(header, 'V47GAIN'),/rem)
      if (g eq 'HIGH') then gain[h] = 0
      if (g eq 'MEDHIGH') then gain[h] = 1
      if (g eq 'MEDLOW') then gain[h] = 2
      if (g eq 'LOW') then gain[h] = 3
   endif
      
   if( get_frameno ) then begin
      frameno[h] = sxpar(header, 'FRAMENO')
   endif
   
   if( get_origxcen ) then begin
      origxcen[h] = sxpar(header, 'ORIGXCEN')
   endif
   
   if( get_origycen ) then begin
      origycen[h] = sxpar(header, 'ORIGYCEN')
   endif
   
   if( get_goodcent ) then begin
      goodcent[h] = sxpar(header, 'GOODCENT')
   endif
   
   if( get_object ) then begin
      object[h] = strcompress(sxpar(header, 'OBJECT'), /rem)
   endif
   
   if( get_dimmfwhm ) then begin
      dimmfwhm[h] = sxpar(header, 'DIMMFWHM')
   endif
   
   if( get_mag1fwhm ) then begin
      mag1fwhm[h] = sxpar(header, 'MAG1FWHM')
   endif
   
   if( get_pixrt ) then begin
      pixrt[h] = sxpar(header, 'V47PIXRT')
   endif
   
endfor

statusline, /clear

end
