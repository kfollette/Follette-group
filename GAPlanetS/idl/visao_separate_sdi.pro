;+
; NAME: visao_separate_sdi
;
; PURPOSE:
;  subtract dark, divide by flat, and separate two channels of sdi data
;
; INPUTS:
;
;
; INPUT KEYWORDS:
; indiv: when datasets are too large to make one huge 1024 x 1024 x nim cube, align and write them to an aligned directory rather than storing in memory
;
; OUTPUTS:
;
;
; OUTPUT KEYWORDS:
;    none
;
; EXAMPLE:
;
;
; HISTORY:
; 2017-08-16 KBF modified to populate headers
; 2019-03-14 KBF modified to handle multiple gains
;-

pro visao_separate_sdi, Line, Cont, avgwfe, rotoff, flat=flat, indiv=indiv, stp=stp, wfe=wfe, suffix=suffix, toelectrons=toelectrons

  visao_inventory, sci_imlist, dark_imlist, flat_imlist, rotoff_sciims, filt, wfe=wfe, mag1=mag1
  ;;create aligned directory if doesn't already exist

  dummy_im=readfits(sci_imlist[0])
  dim1=(size(dummy_im[*,*]))[1]
  dim2=(size(dummy_im[*,*]))[2]
  nims= n_elements(sci_imlist)

  ;;test whether there is already a dark in this directory
  ;y=file_test('master_dark.fits')
  ;if y eq 0 then print, 'You need to make a master dark first - please run visao_dark'
  ;if y eq 1 then master_dark=readfits('master_dark.fits', darkhead, /silent)

  ;;create blank arrays
  if not keyword_set(indiv) then begin
    raw=dblarr(dim1,dim2, nims)
    Line=dblarr(dim1, dim2/2, nims)
    Cont=dblarr(dim1, dim2/2, nims)
  endif else begin
    raw=dblarr(dim1,dim2,1)
    Line=dblarr(dim1, dim2/2,1)
    Cont=dblarr(dim1, dim2/2,1)
    x=file_test('./indiv', /DIRECTORY)
    if x eq 0 then spawn, 'mkdir indiv'
    if x eq 1 then begin 
      ok_str = 'ok'
      ;;count files that begin with the string 'Line'
      spawn, 'ls -l indiv/* | wc -l', nfile
      read, ok_str, prompt='there are currently'+string(nfile)+$
        'files in your indiv directory, and I am going to delete them. OK? [y/n]'
      if(ok_str eq 'y' or ok_str eq 'Y') then spawn, 'rm indiv/*'
      if(ok_str eq 'n' or ok_str eq 'N') then stop
    endif
  endelse
  
  expt=dblarr(nims)
  avgwfe=dblarr(nims)
  rotoff=dblarr(nims)
  object=strarr(nims)
  vfw3posn=strarr(nims)
  dateobs=strarr(nims)
  ra = strarr(nims)
  dec = strarr(nims)
  inst = strarr(nims)
  gain = strarr(nims)
  master_dark=dblarr(1024,1024)
  gainnamelist = ["LOW","MEDLOW","MEDHIGH","HIGH"]
  gainlist = [12.3,3.34,1.77,0.47]

  if keyword_set(flat) then flatim=readfits(flat)

  ;;image loop
  for i=0,nims-1 do begin
    if keyword_set(indiv) then typ='individually' else typ='in bulk'
    status='processing image number'+string(i+1)+'  of'$
           +string(nims)+'  '+string(typ)
    statusline, status, 0
    if not keyword_set(indiv) then begin
      j=i
      raw[*,*,i]=readfits(sci_imlist[i], head, /silent)*1.
    endif else begin
      j=0
      raw[*,*,0]=readfits(sci_imlist[i], head, /silent)*1
    endelse
    dateobs[i]=sxpar(head, 'DATE-OBS')
    expt[i]=sxpar(head, 'EXPTIME')
    ;;round exposure time to nearest hundredth of second to avoid conflicts between gain modes
    expt[i]=round(expt[i]*100)/100.0d
    avgwfe[i]=sxpar(head, 'AVGWFE')
    rotoff[i]=sxpar(head, 'ROTOFF')
    object[i]=sxpar(head, 'OBJECT')
    vfw3posn[i]=sxpar(head, 'VFW3POSN')
    ra[i]=sxpar(head, 'RA')
    dec[i]=sxpar(head, 'DEC')
    inst[i]=sxpar(head, 'INSTRUME')
    gain[i]=sxpar(head, 'V47GAIN')
    
    ;;gains in e-/ADU, used for normalization
    modecheck = sxpar(head, 'V47PIXRT')
    if modecheck ne 250 then print, 'CCD speed not 250kHz - gain table invalid' and stop
    
    ;;if gain has changed, read in new master dark
    if (i eq 0) or (gain[i] ne gain[i-1]) then begin
      print, 'reading new dark for gain ', gain[i]
      master_dark=readfits('master_dark_'+strcompress(gain[i], /REMOVE_ALL)+'.fits', darkhead)
      if keyword_set(toelectrons) then begin
        scalefactor = gainlist[where(gainnamelist eq strcompress(gain[i], /REMOVE_ALL))]
        print, 'scale factor for this gain is ', scalefactor, ' e- per ADU'
      endif else begin
        scalefactor = 1.
        print, 'not converting to electrons per ADU'
      endelse
    endif
    
    ;;some catches in case dark not read properly
    if (size(master_dark))[1] ne 1024 or median(master_dark) lt 300 then print, 'something is wrong here' and stop
    
    ;dark subtract and flat field (if specified)
    if keyword_set(flat) then begin
      ;raw[*,*]=((raw[*,*]-master_dark)/expt[i])/flatim
      raw[*,*,j]=(raw[*,*,j]-master_dark)*scalefactor[0]/flatim
    endif else begin
      ;raw[*,*,i]=((raw[*,*,i]-master_dark)/expt[i])
      raw[*,*,j]=((raw[*,*,j]-master_dark))*scalefactor[0]
    endelse

    ;separate channels
    ;;for Ha top channel is line
    if filt eq 'Ha' then begin
      Line[*,*,j]=raw[0:dim1-1,0:dim2/2-1,j]
      Cont[*,*,j]=raw[0:dim1-1,dim2/2:dim2-1,j]
      ;;for all other filters, bottom is line
    endif else begin
      Cont[*,*,j]=raw[0:dim1-1,0:dim2/2-1,j]
      Line[*,*,j]=raw[0:dim1-1,dim2/2:dim2-1,j]
    endelse

    if keyword_set(indiv) then begin
      if keyword_set(flat) then begin
        writefits, './indiv/Line_flat_'+string(i+1, format='(i04)')+'.fits', Line[*,*,j], head
        writefits, './indiv/Cont_flat_'+string(i+1, format='(i04)')+'.fits', Cont[*,*,j], head
      endif else begin
        writefits, './indiv/Line_'+string(i+1, format='(i04)')+'.fits', Line[*,*,j], head
        writefits, './indiv/Cont_'+string(i+1, format='(i04)')+'.fits', Cont[*,*,j], head
      endelse
    endif

  endfor

  ;;BASIC CHECKS

  ;;check that all images have same exposure time
  if n_elements(uniq(expt)) ne 1 then print, 'WARNING - more than one exposure time in this cube'

  ;;check that darks and images have same exposure time
  dark_expt=sxpar(darkhead, 'EXPTIME')
  
  ;;round exposure time to nearest hundredth of second to be consistent with images
  dark_expt=round(dark_expt)

  if dark_expt ne expt[0] then print, 'WARNING - dark does not have same exposure time as images'

  ;;check that the object is the same in all images
  if n_elements(uniq(object)) ne 1 then print, 'WARNING - more than one object in this cube'

  ;;check that the filter wheel was in the same place in all images
  if n_elements(uniq(vfw3)) ne 1 then print, 'WARNING - more than one SDI filter in this cube'

  ;; add parameters to header
  mkhdr, head_new, Line
  sxaddpar, head_new, 'INSTRUME', inst[uniq(inst)]
  sxaddpar, head_new, 'OBJECT', object[uniq(object)]
  sxaddpar, head_new, 'RA', median(ra)
  sxaddpar, head_new, 'DEC', median(dec)
  if keyword_set(flat) then begin
    sxaddpar, head_new, 'FLAT', flat
  endif
  sxaddpar, head_new, 'EXPTIME', expt[uniq(expt)]
  if not keyword_set(wfe) then wfe='No cut'
  sxaddpar, head_new, 'WFE_CUT', wfe
  sxaddpar, head_new, 'MED_WFE', median(avgwfe)
  sxaddpar, head_new, 'STD_WFE', stdev(avgwfe)
  sxaddpar, head_new, 'VFW3POSN', vfw3posn[uniq(vfw3posn)] 
  ;wavelength keyword required for PyKLIP
  sxaddpar, head_new, 'WLENGTH', 0.656
  sxaddpar, head_new, 'DATE-START', strmid(dateobs[0],0,10)
  sxaddpar, head_new, 'TIME-START', strmid(dateobs[0],11,8)
  sxaddpar, head_new, 'DATE-END', strmid(dateobs[nims-1],0,10)  
  sxaddpar, head_new, 'TIME-END', strmid(dateobs[nims-1],11,8)
  if keyword_set(toelectrons) then unittyp = 'e-' else unittyp = 'ADU'
  sxaddpar, head_new, 'UNITS', unittyp

    if not keyword_set(indiv) then begin
      if keyword_set(flat) then begin
        writefits, 'Line_flat_preproc.fits', Line, head_new
        ;modify WLENGTH before writing continuum header
        sxaddpar, head_new, 'WLENGTH', 0.642
        writefits, 'Cont_flat_preproc.fits', Cont, head_new
      endif else begin
        writefits, 'Line_preproc.fits', Line, head_new
        sxaddpar, head_new, 'WLENGTH', 0.642
        writefits, 'Cont_preproc.fits', Cont, head_new
      endelse
    endif
    writefits, 'avgwfe_preproc.fits', avgwfe
    writefits, 'rotoff_preproc.fits', rotoff
    writefits, 'exptime_preproc.fits', expt

if keyword_set(stp) then stop

print, '        separate sdi complete            '

end
