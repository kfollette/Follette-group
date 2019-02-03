;+
; NAME: clio_getims
;
; PURPOSE:
;  Extracts image type and other headers from all the Clio images in a directory
;
; DESCRIPTION:
;   Generates a list of files in a directory, which can be specified with the subdir keyword, which match the 
;   regex "[prefix]_*.fits".  The default prefix is '', but it can be changed with the prefix keyword (e.g. to 'Clio_').  
;   At minimum returns the filenames of the files.  If the ims keyword is given an argument, 
;   the images are returned as a cube, which by default is fltarr. To avoid re-generating the list of files, and instead use
;   the file names passed in as fnames, set the keyword usefnames.  Various header keyword values can also be extracted.
;
; INPUTS:
;  none
;
; INPUT KEYWORDS:
;   prefix      :  the filename prefix, default is empty ''
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
;          2016-07-10
;
; BUGS/WISH LIST:
;  Supported keywords are incomplete, need to expand.
;
;-
pro clio_getims, prefix, fnames, EXPTIME=EXPTIME, AOLOOPST=AOLOOPST, ROTOFF=ROTOFF, $
                      AVGWFE=AVGWFE, STDWFE=STDWFE, $
                      DATEOBS=DATEOBS, DATE1=DATE1, AOREC=AOREC, UT=UT, AM=AM, $
                      HA=HA, DIMMFWHM=DIMMFWHM, $
                      OBJECT=OBJECT, $
                      ims=ims, prefix=prefix, subdir=subdir, region=region, goodims=goodims, usefnames=usefnames, $
                      double=double

;** First thing we do is generate a list of fits files **

if(n_elements(subdir) eq 0) then begin
   subdir = "."
endif

if(~keyword_set(usefnames) or n_elements(fnames) eq 0) then begin 
   usefnames = 0
   if(n_elements(prefix) eq 0) then prefix=''

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



get_dateobs = 0
if( arg_present(dateobs) ) then begin
   get_dateobs = 1
   dateobs = dblarr(nims)
endif

get_date = 0
if( arg_present(date1) ) then begin
   get_date = 1
   date1 = dblarr(nims)
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

get_dimmfwhm=0
if( arg_present(dimmfwhm) ) then begin
   get_dimmfwhm = 1
   dimmfwhm = dblarr(nims)
endif


get_object = 0
if( arg_present(object) ) then begin
   get_object = 1
   object = strarr(nims)
endif


;*** Allocate output image if desired ***
retims = 0
if( arg_present(ims) ) then begin
   retims = 1
   
   imraw= mrdfits(fnames[0], 0, /silent, /unsigned)
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
      imraw= mrdfits(fnames[h], 0, header, /silent, /unsigned)
      ims[*,*,h] = imraw[rgx[0]:rgx[1],rgx[2]:rgx[3]]
   endif else begin
      fits_read, fnames[h], 0, header, /header_only
   endelse
   

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
   
  
   if( get_rotoff ) then begin
      rotoff[h] = sxpar(header, 'ROTOFF')
   endif
   
   if( get_rotang ) then begin
      rotang[h] = sxpar(header, 'ROTANG')
   endif
   
   
   if( get_avgwfe ) then begin
      avgwfe[h] = sxpar(header, 'AVGWFE')
   endif
   
   if( get_stdwfe ) then begin
      stdwfe[h] = sxpar(header, 'STDWFE')
   endif
   
   
   if( get_dateobs ) then begin
      dobs = sxpar(header, 'DATE-OBS')
      dateobs[h] = fitsdate_conv(dobs, 'JULIAN')
   endif
   
   if( get_date ) then begin
      dobs = sxpar(header, 'DATE')
      date1[h] = fitsdate_conv(dobs, 'JULIAN')
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
   
   if( get_dimmfwhm ) then begin
      dimmfwhm[h] = sxpar(header, 'DIMMFWHM')
   endif
   
   if( get_aorec ) then begin
      aorec[h] = sxpar(header, 'AOREC')
   endif
   
      
   
   if( get_object ) then begin
      object[h] = strcompress(sxpar(header, 'OBJECT'), /rem)
   endif
   
endfor

statusline, /clear

end
