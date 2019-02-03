; $Id: image_shift.pro, v 1.0 Aug 1999 e.d. $
;
;+
; NAME:
;IMAGE_SHIFT
;
; PURPOSE:
;Interpolate a 2D image in order to perform a fractional
;shift of the origin.
;
; CATEGORY:
;Mathematics. Interpolation.
;
; CALLING SEQUENCE:
;Result = IMAGE_SHIFT(Image, X_shift, Y_shift, Data)
;
; INPUTS:
;Image:2D array to be shifted
;
;X_shift, Y_shift:Components of the fractional shift
;
; OPTIONAL INPUTS:
;Data:See OPTIONAL OUTPUTS
;
; KEYWORD PARAMETERS:
;INTERP_TYPE:Use this keyword to choose an interpolation technique.
;The supported options are
;'F': interpolation by Fourier transform shift
;'S': interpolation by spline functions
;'I': interpolation by the IDL function INTERPOLATE (default)
;
; OPTIONAL OUTPUTS:
;Data:Structure containing useful information to be used on input
;if another shift of the same Image array must be performed.
;This allows the user to save some computation time
;
; RESTRICTIONS:
;1) Interpolation is suited to well-sampled data. For undersampled
;images,
;other techiques should be used. For this purpose, see the function
;STARS in the file 'stars.pro'.
;2) This routine may perform a fractional shift, i.e. abs(X_shift) <=
;0.5,
;abs(Y_shift) <= 0.5. When the fractional shift exceeds 0.5 pixels,
;strong
;edge effects might occur.
;
; MODIFICATION HISTORY:
; Written by:Emiliano Diolaiti, August 1999.
;-

FUNCTION image_shift, image, x_shift, y_shift, INTERP_TYPE = interp_type, data

  on_error, 2
  if  x_shift eq 0 and y_shift eq 0  then  return, image
  ; is AUX_DATA defined?
  no_data = n_tags(data) eq 0
  ; some operations if the auxiliary data are undefined
  if  no_data  then begin
        ; define interpolation type to use
        if  n_elements(interp_type) eq 0  then  interp_type = 'I'
           interp = strupcase(strmid(interp_type, 0, 1))
              if  interp ne 'F' and interp ne 'S' and interp ne 'I'  then  interp = 'I'
                 if  interp eq 'S' and strlen(interp_type) eq 7  then $
                          degree = fix(strmid(interp_type, 6, 1))  else  degree = 3
                    ; extend image to prevent edge effects
                    siz = size52(image, /DIM)  &  extsiz = siz + 2
                       imag = extend_array(image, extsiz[0], extsiz[1], OFFSET = offset)
                          if  interp eq 'F'  then begin
                                   imag = extend_array(imag, sx, sy, /POW2, OFFSET = add_off)
                                         extsiz = [sx, sy]  &  offset = offset + add_off
                                      endif
                             lo = offset  &  up = lo + siz - 1
                          endif else  interp = data.interp_type

  ; interpolate
  case  interp  of

     'S': begin; SPLINE
        if  no_data  then begin
              spline_coeff, imag, DEGREE = degree, $
                             coefficients, x_knots, y_knots, x, y
                 data = {coefficients: coefficients, $
                                    x_knots: x_knots, y_knots: y_knots, $
                                    x: x, y: y, degree: degree, lo: lo, up: up, $
                                    interp_type: interp}
              endif
        imag = spline_interp(data.coefficients, data.x_knots, $
                                  data.y_knots, DEGREE = data.degree, $
                                  data.x - x_shift, data.y - y_shift)
     end

     'F':begin; Fourier Transform
        if  no_data  then $
              data = {image_ft: fft(imag), extsiz: extsiz, $
                         lo: lo, up: up, interp_type: interp}
        phi = phase2d(x_shift, y_shift, data.extsiz[0], data.extsiz[1])
        imag = float(fft(data.image_ft * phi, /INVERSE))
     end

     'I': begin; IDL INTERPOLATE
        if  no_data  then $
              data = {image: imag, $
                         x: findgen(extsiz[0]), y: findgen(extsiz[1]), $
                         lo: lo, up: up, interp_type: interp}
        imag = interpolate(data.image, data.x - x_shift, data.y - y_shift, $
                              /GRID, CUBIC = -0.5, MISSING = 0.)
     end

  endcase

  return, imag[data.lo[0]:data.up[0],data.lo[1]:data.up[1]]
end
