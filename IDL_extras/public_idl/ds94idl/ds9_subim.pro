pro ds9_subim, x, y, z, width, subim
;extracts a sub-image from ds9image


common ds9env, $    ;this is the ds9 environment common block
       ds9obj    ;the ds9 object


ds9obj->getimage, ds9im
       
if ( (size(ds9im))[0] ne 3) then begin
   subim =  ds9im[floor(x-0.5*width):floor(x-0.5*width)+width-1, $          
         floor(y-0.5*width):floor(y-0.5*width)+width-1]  
endif else begin
   subim =  ds9im[floor(x-0.5*width):floor(x-0.5*width)+width-1, $          
         floor(y-0.5*width):floor(y-0.5*width)+width-1, z]
endelse

end

