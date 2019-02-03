pro saoct, ctname, ncol=ncol, red, green, blue, get=get


if(n_elements(ncol) ne 1) then ncol = 256


if(ctname eq 'grey') then saoct_grey, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'red') then saoct_red, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'green') then saoct_green, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'blue') then saoct_blue, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'a') then saoct_a, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'b') then saoct_b, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'bb') then saoct_bb, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'he') then saoct_he, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'i8') then saoct_i8, red_x, red_c, green_x, green_c, blue_x, blue_c

if(ctname eq 'heat') then saoct_heat, red_x, red_c, green_x, green_c, blue_x, blue_c

if(ctname eq 'cool') then saoct_cool, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'coolgreen') then saoct_coolgreen, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'rainbow') then saoct_rainbow, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'standard') then saoct_standard, red_x, red_c, green_x, green_c, blue_x, blue_c
if(ctname eq 'sls') then saoct_sls, red_x, red_c, green_x, green_c, blue_x, blue_c




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

