pro matlabct, ctname, ncol=ncol, red, green, blue, get=get


if(n_elements(ncol) ne 1) then ncol = 256


if(ctname eq 'bone') then matlabct_bone, red_x, red_c, green_x, green_c, blue_x, blue_c




x = findgen(ncol)

linterp, red_x, red_c*(ncol-1.), x/(ncol-1), red
linterp, green_x, green_c*(ncol-1.), x/(ncol-1), green
linterp, blue_x, blue_c*(ncol-1.), x/(ncol-1), blue

;red = smooth(red, 55)
;green = smooth(green, 55)
;blue = smooth(blue, 55)
;red = spline(red_x, red_c*(ncol-1.), x/(ncol-1), 10.)
;green = spline(green_x, green_c*(ncol-1.), x/(ncol-1), 10.0)
;blue = spline(blue_x, blue_c*(ncol-1.), x/(ncol-1), 10.0)

if ~keyword_set(get) then begin
   tvlct, red, green, blue
endif

end

