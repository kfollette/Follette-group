;+
; NAME: 
;   ds9_imexam
;
; PURPOSE:
;   image examination using the ds9 image viewer
;
; DESCRIPTION:
;   Using keystrokes inside the ds9 image viewer, calculates various statistics, photometry, and/or plots image data.
;   The keystrokes are:
;         m  -  calculates and displays number of pixels, min, max, mean, median, stand. dev.
;         a  -  not implemented yet
;         c  -  centroids using gcntrd, and presents location and sep/pa from image center
;         j
;         k
;         r  -  radial profile plot.  first centroids.
;
; INPUTS:
;   none
; 
; INPUT KEYWORDS:
;   image : an image to display and examine
;
; OUTPUTS:
;   none
;
;
; MODIFICATION HISTORY:
;  Written 2013/04/12 by Jared Males (jrmales@email.arizona.edu)
;  2013/05/26 added centroiding (Jared Males)
;
; BUGS/WISH LIST:
;   Need to add fitting  (a, j, k, r)
;   Need to add photometry (a)
;-
pro ds9_imexam, image


common ds9env, $    ;this is the ds9 environment common block
       ds9obj    ;the ds9 object

;Check if initialized
ds9_init

;Display the image if supplied
if(n_elements(image) gt 1) then ds9, image

stats = dblarr(6)

key = ''

;controls whether the header for stats is printed or not.
headline = 0


ds9obj->getimage, ds9im

B = size(ds9im, /dim)

while key ne 'q' do begin

   ds9obj->basic_imexam, tkey, x, y, z
   key = tkey

   if key eq 'm' then begin
   
      ds9_subim, x, y, z, ds9obj->statsw(), subim
      
      stats[0] = n_elements(subim)
      stats[1] = min(subim)
      stats[2] = max(subim)
      stats[3] = mean(subim)
      stats[4] = median(subim)
      stats[5] = stddev(subim)
      
      if (headline ne 1) then begin
         print, 'Statistics:'
         print, format='(A8,A8,A8,A8,A8,A8)', '#pix','min', 'max', 'mean', 'median', 'stddev'
         headline = 1
      endif
      
      print, format='(F8.0,F8.1,F8.1,F8.1,F8.1,F8.1)', stats[0],stats[1],stats[2],stats[3],stats[4],stats[5]
   endif
   
   if key eq 'j' then begin
   
      ds9_subim, x, y, z, ds9obj->plotw(), subim
      
      plot, floor(y-0.5*ds9obj->plotw()) + findgen(ds9obj->plotw()), subim[0.5*ds9obj->plotw(), *], xrange=[y-.5*ds9obj->plotw(),y+.5*ds9obj->plotw()], /xst, xtitle='Row', ytitle='Value'
      empty
   endif
   
   if key eq 'k' then begin
   
      ds9_subim, x, y, z, ds9obj->plotw(), subim
            
      plot, floor(x-0.5*ds9obj->plotw())+findgen(ds9obj->plotw()), subim[*,.5*ds9obj->plotw()], xrange=[x-.5*ds9obj->plotw(),x+.5*ds9obj->plotw()], /xst, xtitle='Column', ytitle='Value'
      empty
   endif
   
   if key eq 'r' then begin
   
      corefit, ds9im, xcen, ycen, peak, fwhm, xguess=x, yguess=y, fwguess=ds9obj->fwhm(), fitim=fitim, fitrad=fitrad, subim=subim, skyrad = 50., skywidth=5., /circular
       
      print, format='(A8,A8,A8,A8)', 'xcen', 'ycen', 'peak', 'fwhm'
      print, xcen, ycen, peak, fwhm
      
      plot, fitrad, subim,  psym=1, xtitle='Radius (Pixels)', ytitle='Value'
      idx = sort(fitrad)
      oplot, fitrad[idx], fitim[idx]
      
      empty
   endif
   
   if key eq 'c' then begin
      ds9_centroid, xout, yout, sep, pa, x, y, z
      
      if (headline ne 2) then begin
         print, 'Centroid:'
         print, format='(A8,A8,A8,A8)', 'x','y', 'sep', 'PA'
         headline = 2
      endif
      
      print, xout, yout, sep, pa, format='(F8.2,F8.2,F8.2,F8.2)'
   endif
   
endwhile

end
