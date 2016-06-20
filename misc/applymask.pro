pro applymask, im, mask, maskedim, read=read, write=write, nan=nan

  ;; replaces pixels under mask with zeroes
  ;; mask should contain 1s in pixels to be masked
  ;; read keyword = set if you are inputting .fits files for im and mask
  ;; if nan set, replaces with nans, not zeroes

if keyword_set(read) then begin
  im=readfits(string(im)+'.fits')
  maskim=readfits(string(mask)+'.fits')
endif

  xdim=(size(im))[1]
  ydim=(size(im))[2]

  for x=0, xdim-1 do begin
    for y=0, ydim-1 do begin
      if keyword_set(nan) then begin
        if mask[x,y] eq 1 then im[x,y]='NaN'
      endif else begin
      if mask[x,y] eq 1 then im[x,y]=0
      endelse
    endfor
  endfor

maskedim=im

if keyword_set(write) then begin
writefits, string(file)+'_mask.fits', im
endif

end