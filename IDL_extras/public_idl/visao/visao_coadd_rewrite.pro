;+
; NAME: 
;   visao_coadd_rewrite
;
; DESCRIPTION:
;   coadds a set of images, and writes the results to disk with filenames coadd_XXXXX.fits.  coadds sequential images
;   until a maximum elapsed time or maximum change in angle is reached.  Unless meancomb is set, images are median combined.
;
; INPUTS:
;   outpath   : the output path
;   fnames    : string array containing the names of files to coadd
;   dateobs   : vector of date-obs, in units of days
;   rotoff    : vector of rotator offsets, in degrees
;   exptime   : the exposure time of these images, scalar
;   maxdt     : the maximum delta-t to coadd across (seconds)
;   maxdq     : the maximum delta-angle to coadd across (degrees)
; 
; INPUT KEYWORDS:
;   num      : if set, then just coadds in chunks of num images [NOT IMPLEMENTED YET]
;   meancomb : if set, then images are mean combined instead of median
;   region   : only coadd and write a specific region.  if not set, whole image is coadded
;   xregion  : if set, images are registered to first image in each coadd. if xregion is a 4 element vector, it specifies
;              the region to use for registration.
;   xstar    : the x shifts applied in registration
;   ystar    : the y shifts applied in registration
;
; OUTPUTS:
;   none
;
;
; MODIFICATION HISTORY:
;  Written 2013/01/20 by Jared Males (jrmales@email.arizona.edu)
;
; BUGS/WISH LIST:
;
;-
pro visao_coadd_rewrite, outpath, fnames, dateobs, rotoff, exptime, maxdt, maxdq, num=num, meancomb=meancomb, dark=dark,$
                          region=region, xregion=xregion, xstar=xstar, ystar=ystar, $
                          fakeinj=fakeinj, fakesep=fakesep, fakepa=fakepa, fakefw=fakefw, fakeadu=fakeadu

domean = 0
if(keyword_set(meancomb)) then domean = 1

nims = n_elements(fnames)

curr_im = 0l
im_num = 0l

;normalize the angles
rqrotoff = requad_angles(rotoff)

testim = double(mrdfits(fnames[0], /silent))

dim1 = (size(testim))[1]
dim2 = (size(testim))[2]
   
if(n_elements(region) ne 4) then begin
   region = [0, dim1, 0, dim2]
endif else begin
   dim1 = region[1]-region[0]
   dim2 = region[3]-region[2]
endelse

if(n_elements(dark) eq 0) then begin
   dark = fltarr(dim1, dim2)
endif

;if xregion is set and is only single nonzero value, turn it into a 4 vector
if(n_elements(xregion) eq 1) then begin
   if(keyword_set(xregion)) then begin
      xregion=[0., dim1, 0., dim2]
   endif
endif
   
while (curr_im lt nims) do begin

   included_files = fnames[curr_im]
   ims = float(mrdfits(fnames[curr_im], 0, header, /silent))
   
   ims = ims[region[0]:region[1]-1, region[2]:region[3]-1]
      
   avg_dateobs = dateobs[curr_im]
   avg_rotoff = rqrotoff[curr_im]
   totexp = exptime
   
   start_im = curr_im

   curr_im = curr_im + 1

   if(curr_im lt nims) then begin
   
   while((dateobs[curr_im] - dateobs[start_im] lt maxdt/86400.) and (rotoff[curr_im] - rotoff[start_im] lt maxdq)) do begin
   
      status = strcompress(string(curr_im+1) + '/' + string(nims), /rem) + ' (' + string(im_num) + ' written)'
      statusline, status, 0
      
      included_files = [included_files, fnames[curr_im]]

      ;read in the next image, and extract the region   
      im = float(mrdfits(fnames[curr_im], /silent))
      im = im[region[0]:region[1]-1, region[2]:region[3]-1]
                
      ims = [[[ims]], [[im]]]
 
      avg_dateobs = avg_dateobs + dateobs[curr_im]
      avg_rotoff = avg_rotoff + rqrotoff[curr_im]
      totexp = totexp + exptime
   
      curr_im = curr_im + 1
      
      
      if(long(curr_im) ge nims) then break
      
   endwhile
   endif
   
   avg_dateobs = avg_dateobs / double(n_elements(included_files))
   avg_rotoff  = avg_rotoff  / double(n_elements(included_files))

   if(n_elements(xregion) eq 4 and (size(ims))[0] eq 3 ) then begin
      visao_reg_cube, ims, xst, yst,rbox=xregion, /doshift
      if(n_elements(xstar) eq 0) then begin
         xstar = xst
         ystar = yst
      endif else begin
         xstar = [xstar, xst]
         ystar = [ystar, yst]
      endelse
   endif
   
   if(domean eq 0) then begin
      if((size(ims))[0] eq 3) then begin
         outim = median(ims, dim=3)
      endif else begin
         outim = ims
      endelse
   endif else begin
      if((size(ims))[0] eq 3) then begin
         outim = total(ims,3) / double(n_elements(included_files))
      endif else begin
         outim = ims;
      endelse
   endelse

   outfile = 'coadd_' + string(im_num, format='(I05)')
   outname = strcompress(outpath + '/' + outfile + '.fits')

   
   sxaddpar, header, 'ROTOFF', avg_rotoff
   sxaddpar, header, 'DATE-OBS', fitsdate_conv(avg_dateobs, 'F')
   sxaddpar, header, 'NCOADDS', n_elements(included_files)
   sxaddpar, header, 'TOTEXP', totexp
   
   if(domean eq 0) then begin
      sxaddpar, header, 'HISTORY', 'files median combined by visao_coadd_rewrite.pro:'
   endif else begin
      sxaddpar, header, 'HISTORY', 'files mean combined by visao_coadd_rewrite.pro:'
   endelse
   for h=0,n_elements(included_files)-1 do sxaddpar, header, 'HISTORY', included_files[h]
   
   
   if(keyword_set(fakeinj)) then begin
      
      if(n_elements(fakeadu gt 1)) then begin

         for z=0,n_elements(fakeadu)-1 do begin
         
            adi_psfmodel, model, dim1, dim2, visao_getderot(avg_rotoff), fakesep[z], fakepa[z], fwhm=fakefw
      
            outim = outim+model*fakeadu[z]
            
         endfor
      endif else begin
      
         adi_psfmodel, model, dim1, dim2, visao_getderot(avg_rotoff), fakesep, fakepa, fwhm=fakefw
      
         outim = outim+model*fakeadu
         
      endelse
   endif
   
   mwrfits, outim-dark, outname, header, /silent

   
   im_num = im_num + 1
   
 
endwhile

statusline, /clear

print,  'visao_coadd_rewrite: ' + strcompress(string(im_num),/rem) + ' images written.'

end     
      









