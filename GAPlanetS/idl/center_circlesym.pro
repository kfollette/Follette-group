pro center_circlesym, im, xr, yr, rmax, xc, yc, grid, mask=mask
;+
; Finds the center of a star image assuming circular symmetry
; im   :  the input image
; xr   : vector of x-offsets from the center of the image
; yr   : vector of y-offsets from the center of the image
; rmax : maximum radius to consider
; xc   : the center x position
; yc   : the center y position
; grid : grid of results (xr vs yr vs stddev)
; mask : optional 1/0 mask, 0 specifies pixels to ignore.
;-

   get_cubedims, im, dim1, dim2
   ;print, dim1, dim2
   
   grid = fltarr(n_elements(xr), n_elements(yr))
   for i=0,n_elements(xr)-1 do begin
      for j=0,n_elements(yr)-1 do begin
         ;print, i, j, format='((x, I0),$)'
         if dim1 ne dim2 then r = rarr(dim1, dim2, xc=xr[i], yc=yr[j]/2, /pix) else r = rarr(dim1, dim2, xc=xr[i], yc=yr[j], /pix)
         ;stop
         for k=0,rmax do begin
          status = string(i) + ' /' + string(n_elements(xr)) + $
            string(j) + ' /' + string(n_elements(yr)) + $
            string(k) + ' /' + string(rmax)
          statusline, status, 0
            if(n_elements(mask) eq n_elements(im)) then begin
               idx = where( r ge k and r lt k+1 and mask gt 0 )
            endif else begin
               idx = where( r ge k and r lt k+1 )
            endelse
            
            if(idx[0] eq -1) then continue
            sd = stddev(im[idx], /nan)
            if(~finite(sd)) then sd = 0
            ;; modified to be sum of squares (squared second term) to match Jared's most recent C version
            grid[i,j] = grid[i,j] + (sd/abs(median(im[idx])))^2
         endfor
      endfor
   endfor

   writefits, 'circsym_grid_idl.fits', grid
   ming = min(grid, idx)
   pos = array_indices(grid, idx)
   
   
   gcntrd, -1*grid, pos[0], pos[1], xcc, ycc, 0.5*n_elements(xr)
   
   ;print, xcc
   ;print, ycc
   xc = xr[0] + xcc*(xr[1]-xr[0]) + 0.5*(dim1-1)
   yc = yr[0] + ycc*(yr[1]-yr[0]) + 0.5*(dim2-1)
   
;   xc = xr[xcc+0.5] + 0.5*(dim1-1)
;   yc = yr[ycc+0.5] + 0.5*(dim2-1)


end


