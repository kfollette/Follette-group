; $Id: size52.pro, v 1.0 Aug 1999 e.d. $
;
;+
; NAME:
;SIZE52
;
; PURPOSE:
;Alias of IDL 5.2 intrinsic SIZE function for previous versions.
;
; CATEGORY:
;Array informational routines.
;
; CALLING SEQUENCE:
;Result = SIZE52(X)
;
; INPUTS:
;X:IDL variable
;
; KEYWORD PARAMETERS:
;N_DIMENSION:Set this keyword to a nonzero value to retrieve the
;number of dimensions of X
;
;DIMENSION:Set this keyword to a nonzero value to retrieve a
;long-integer vector containing the size of each dimension of X
;
;TYPE:Set this keyword to a nonzero value to retrieve the IDL
;type code of X
;
; OUTPUTS:
;Result:Same as output of SIZE if no keyword is set. Otherwise return
;the result specified by the KEYWORD.
;
; MODIFICATION HISTORY:
; Written by:Emiliano Diolaiti, August 1999.
;-

FUNCTION size52, x, N_DIMENSION = ndim, DIMENSION = dim, TYPE = type

  s = size(x)
  if  keyword_set(ndim)  then  s = s[0]  else $
     if  keyword_set(dim)  then  begin
        if  s[0] le 1  then  s = s[0]  else  s = s[1:s[0]]
     endif  else $
        if  keyword_set(type)  then  s = s[s[0] + 1]
  return, s
end
