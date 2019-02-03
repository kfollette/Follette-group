function pixel_circle_intcpt, x, y, r

w = 0

;intercept(s) of top
rp = sqrt(r^2-(0.-(y+0.5))^2)
if(rp ne 0.) then begin

   xcirc = rp - (x-0.) + .5

   if(xcirc ge 0. and xcirc le 1.) then begin
      w = [w, xcirc, 1.]
   endif

   xcirc = -rp - (x-0.) + .5

   if(xcirc ge 0. and xcirc le 1.) then begin
      w = [w, xcirc, 1.]
   endif

endif

;intercept(s) of bottom
rp = sqrt(r^2-(0.-(y-0.5))^2)
if(rp ne 0.) then begin

   xcirc =  rp - (x-0.) + .5

   if(xcirc ge 0. and xcirc le 1.) then begin
      w = [w, xcirc, 0.]
   endif

   xcirc = -rp - (x-0.) + .5

   if(xcirc ge 0. and xcirc le 1.) then begin
      w = [w, xcirc, 0.]
   endif
endif

;intercept of left
rp = sqrt(r^2-(0.-(x+0.5))^2)
if(rp ne 0) then begin

   xcirc =  rp - (y-0.) + .5

   if(xcirc ge 0. and xcirc le 1.) then begin
      w = [w, 1., xcirc]
   endif

   xcirc = -rp - (y-0.) + .5

   if(xcirc ge 0. and xcirc le 1.) then begin
      w = [w, 1., xcirc]
   endif
endif


;intercept of right
rp = sqrt(r^2-(0.-(x-0.5))^2)
if(rp ne 0) then begin

   xcirc = rp - (y-0.) + .5

   if(xcirc ge 0. and xcirc le 1.) then begin
      w = [w, 0., xcirc]
   endif


   xcirc = -rp - (y-0.) + .5

   if(xcirc ge 0. and xcirc le 1.) then begin
      w = [w, 0., xcirc]
   endif

endif

if(n_elements(w) eq 1) then w = [w, -1.,-1.,-1.,-1.]
if(n_elements(w) eq 3) then w = [w, -1.,-1.]

return, w[1:*]


end
