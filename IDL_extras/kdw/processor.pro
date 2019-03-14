
;script for easy running of processv2
;see processv2.pro for extra details
;
;just add lines or comment line out depending on what you want to do.
;
;NO INPUTS FOR PROCESSOR.PRO
;INPUTS FOR PROCESSV2 (for making a new line)
;
;start file number: found in headstrip .txt files in the data
;                   directories, /raid/DATE on pavo or /raid/chara/DATE. 
;end file number  : found in headstip or file containong data
;output text file : invent your own output text file name
;
;OPTIONAL INPUT FOR PROCESSV2.PRO: if files are of form pavoXXXXX.fits
;and NOT pavoXXXXX.fits.Z then add ,/z to the call to processv2
;    ---
;;%%%%%%%%%format
;processv2,'filename',[tels,*,*]startfileno,endfileno,output text file name, forgroundfile,darkfile 
;;%%%%%%%%%%%%%%

; include main subroutines in compiling
@processv2.pro
@ps_from_files.pro
@ppspec.pro
@slicedice.pro
@find_obs_status.pro
@find_rot_offset.pro
@loadobs.pro
@saveobs.pro

tstart = systime(1) 

; typical 2-tel PAVO@CHARA example
;processv2,'110814',[0,1,1],00000,42693, '110814.dat',foregroundtype=1,root_dir='/raid/chara/',outdir='110814/',/individual,/plot

; typical 3-tel PAVO@CHARA example
;;processv2, '090714', [1,1,1], 13163, 14580, '090714_nopol_test.txt', root_dir='/import/pendragon3/snert/chara/', /plot, /individual

;;;;;;;;ATLAS: CP only example;;;;;;; 
;processv2, '081117', [1,1,1], 00494, 00831, '081117_a.txt', /hann, root_dir='/import/pendragon3/snert/chara/', /plot, /indiv

; typical PAVO@SUSI example
;processv2,'100714',[1,1],1790,4565, 'susi100714.txt',root_dir='/import/pendragon1/aaron/',/plot,/susi, /individual


tend = systime(1)
print, "total time:", tend-tstart

end

