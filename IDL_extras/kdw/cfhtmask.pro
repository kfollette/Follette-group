function cfhtmask,image

xcoord=INTARR(1024,1024)
ycoord=INTARR(1024,1024)

FOR x=0,1023 DO BEGIN
   xcoord[x,*]=x
   ycoord[*,x]=x
ENDFOR

image(WHERE(xcoord le 20))=!VALUES.F_NAN
image(WHERE(xcoord ge 1010))=!VALUES.F_NAN
image(WHERE(ycoord le 40))=!VALUES.F_NAN
image(WHERE(ycoord ge 1015))=!VALUES.F_NAN
image(WHERE(((xcoord ge 75) AND (xcoord le 160)) AND ((ycoord ge 400) AND (ycoord le 450))))=!VALUES.F_NAN
image(WHERE((xcoord le 70) AND (ycoord le 100)))=!VALUES.F_NAN
image(WHERE((ycoord le 80) AND (xcoord ge 280) AND (xcoord le 450)))=!VALUES.F_NAN 

RETURN, image

END
