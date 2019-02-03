;+
; NAME: psplot_stop
;
; PURPOSE:
;  Stop postscript plotting.
;
; DESCRIPTION:
;  Closes the file and sets the device to 'x'.  Also changes the banner to match the filename instead of
;  the idl graphics default.
;
; INPUTS:
;  none
;
; INPUT KEYWORDS:
;  none
;
; OUTPUTS:
;  none
;;
; HISTORY:
;  Written 2013 by Jared Males, jrmales@email.arizona.edu
;
;-
pro psplot_stop

;This is set by the call to psplot_start
common psplot, filename

device, /close_file

spawn,'cat '+filename+$                               ; replace irritating
     '| sed "s|Graphics produced by IDL|'+filename+$  ; IDL plot banner
     '|" >  idltemp.ps; mv idltemp.ps '+filename      ; with the file name 
  
  
set_plot, 'x'

end

