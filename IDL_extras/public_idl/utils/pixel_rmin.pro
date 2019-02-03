function pixel_rmin, xp, yp, xc, yc

rmin = xp*0. - 1e25 ;allocation


idx = where(abs(xp-xc) le 0.5)
if(idx[0] ne -1) then rmin[idx] = (abs(yp[idx]-yc) - 0.5)

idx = where(abs(yp-yc) le 0.5)
if(idx[0] ne -1) then rmin[idx] = (abs(xp[idx]-xc) - 0.5)

idx = where(rmin eq - 0.5)
if(idx[0] ne -1) then rmin[idx] = 0. ;centered

idx = where(rmin eq - 1e25)
if(idx[0] ne -1) then rmin[idx] = $
sqrt( (xp[idx] - 0.5*sign(xp[idx]-xc) - xc )^2 +  (yp[idx] - 0.5*sign(yp[idx]-yc) - yc )^2) ;The nearest corner.
;sqrt( (xp[idx] - 0.5*sign(xp[idx]-xc) - xc )^2 +  (yp[idx] - 0.5*sign(yp[idx]-yc) - yc )^2) ;The nearest corner.


return, rmin

end
