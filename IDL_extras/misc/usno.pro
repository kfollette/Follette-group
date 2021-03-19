pro usno, ramin,ramax, decmin,decmax, sepmin,sepmax,diffmin,diffmax,year


  ;;download new files from http://www.usno.navy.mil/USNO/astrometry/optical-IR-prod/wds/orb6
  ;;save as tab delimited text files with nothing (no space, comma, or other symbol) inside blank cells
  ;;note that some columns are removed in following format (e.g. ADS and HIP columns from raw table)
  ;year is current date in YYYY.YY form

  readcol, '~/Desktop/Science/IDL_extras/my_stuff/orb6orbits_new_match.txt', RA, Dec, WDS, DD, HD, Vp, Vp_flag,Vs, Vs_flag,  per, peru, pere, sma, smau, smae, $
    i, ie, omega, omegae, tperi, tperiu, tperie, $
    format='D,D,A,A,F,F,A,F,A,F,A,F,F,A,F,F,F,F,F,F,A,A', /preserve_null, skipline=2, delimiter=string(9b), /debug;  $
    
  readcol, '~/Desktop/Science/IDL_extras/my_stuff/orb6ephem_new.txt', WDS2, name, grade2, ref2, th2014, rho2014, th2015, rh2015, th2016, rh2016, th2017, $
    rh2017, th2018, rh2018, note, format='A,A,I,A,F,F,F,F,F,F,F,F,F,F,A', /preserve_null, skipline=2, delimiter=string(9b),  /debug
    
  match, WDS, WDS2, mtch1, mtch2
  good_list=strarr(n_elements(mtch1),11)
  k=0
  goods=fltarr(n_elements(mtch1))
  sep=fltarr(n_elements(mtch1), 4)
  PA=fltarr(n_elements(mtch1), 4)
  
  for i=0, n_elements(mtch1)-1 do begin
    if (RA[mtch1[i]] gt ramin) and (RA[mtch1[i]] lt ramax) and (Dec[mtch1[i]] gt decmin) and (Dec[mtch1[i]]) lt decmax $
      and (rh2016[mtch2[i]] gt sepmin) and (rh2016[mtch2[i]] lt sepmax) and (abs(Vs[mtch1[i]]-Vp[mtch1[i]]) lt diffmax) $
      and (abs(Vs[mtch1[i]]-Vp[mtch1[i]]) gt diffmin) and (Vp_flag[mtch1[i]] ne 'NaN') and (Vs_flag[mtch1[i]] ne 'NaN') then begin
      k=k+1
      print, k
      good_list[i,*]=[WDS[mtch1[i]],WDS2[mtch2[i]], string(HD[mtch1[i]]),string(RA[mtch1[i]]), string(Dec[mtch1[i]]), $
        string(Vp[mtch1[i]]), Vp_flag[mtch1[i]],string(Vs[mtch1[i]]),Vs_flag[mtch1[i]],$
        string(rh2016[mtch2[i]]),string(th2016[mtch2[i]])]
        print, WDS[mtch1[i]], '      ', WDS2[mtch2[i]]
      sep[i,0]=rh2015[mtch2[i]]
      sep[i,1]=rh2016[mtch2[i]]
      sep[i,2]=rh2017[mtch2[i]]
      sep[i,3]=rh2018[mtch2[i]]
      PA[i,0]=th2015[mtch2[i]]
      PA[i,1]=th2016[mtch2[i]]
      PA[i,2]=th2017[mtch2[i]]
      PA[i,3]=th2018[mtch2[i]]
      goods[i]=i
    endif else begin
      good_list[i,*]='reject'
    endelse
    
    good_index=where(goods ne 0)
    cut_list=strarr(n_elements(good_index)+1,11)
    cut_list[0,*]=['WDS', 'WDS check', 'HD number', 'RA', 'Dec', 'mag prim', '', 'magsec', '', '2016 rho', '2016 theta']
    dates=[2015, 2016, 2017, 2018]
    sep_use=fltarr(n_elements(good_index)+1)
    PA_use=fltarr(n_elements(good_index)+1)
    
    for j=1, k do begin
      cut_list[j,*]=good_list[good_index[j-1],*]
      sep_use[j]=(interpol(sep[good_index[j-1],*],dates,year))
      PA_use[j]=(interpol(PA[good_index[j-1],*],dates,year))
    endfor
   
    sep_print=string(sep_use)
    PA_print=string(PA_use)
    sep_print[0]='current rho'
    PA_print[0]='current theta'
    
    fmt='(A20,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)'
    forprint, f=fmt, $
      cut_list[*,0], cut_list[*,1], cut_list[*,2], cut_list[*,3], cut_list[*,4], cut_list[*,5], $
      cut_list[*,7], cut_list[*,9], cut_list[*,10], sep_print[*], PA_print[*], textout='USNO_parsed.txt', width=1000, /silent
      
  endfor
  
  
  stop
end