pro radprof,image,xcen,ycen,radius,oplot=oplot1,expand=expand1, $
  gaussvol=gaussvol,gaussparams=gaussparams,miniplot=miniplot,silent=silent, $
  noplot=noplot,guiderplot=guiderplot
;+
;NAME:
;      RADPROF
;PURPOSE:
;      Plot a radial profile of a (stellar) object and overplot a gaussian
;      fit.  Also prints out some fit parameters, e.g. Sky, FWHM, Max, Total.
;CALLING SEQUENCE:
;      RADPROF,image,xcen,ycen,radius,/oplot,expand=nn
;INPUTS:
;      IMAGE  - Input image
;      XCEN   - Centroid (exact!) X position of star (scalar!)
;      YCEN   - Centroid (exact!) Y position of star (scalar!)
;      RADIUS - Radius to examine (and fit) star (scalar)
;OPTIONAL INPUTS:
;      OPLOT  - Do an OPLOT instead of a PLOT
;      EXPAND - Does a cubic interpolative expansion of the object before
;                 doing the radial plot.  This has the advantage of generating
;                 more data points, but is not particularly honest.  It
;                 also doesn't gain you much, and probably should be avoided.
;OUTPUTS:
;      Plots a radial profile to graphics channel with overlay fit.
;      Prints in a row:
;        X,Y  - X,Y centroid of star (AS GIVEN!) This must be exact
;        Sky  - Sky value derived from fit.  may be too high if radius too smal
;        FWHM - Derived FWHM from fitted gaussian
;        Max  - Derived Maximum of fitted gaussian (not max pixel value!)
;        Total- Derived total volume of fitted gaussian 
;PROCEDURE:
;      Generate a plot of pixel value versus radius of the star.  Symmetrize
;      the plot and fit a gaussian using USERLIB GAUSSFIT.  Plot final fit
;      and print out some useful fitted parameters
;MODIFICATION HISTORY
;      06-JUN-94  Written by Eric W. Deutsch
;-

  if (n_params(0) lt 4) then begin
    print,'Call> RADPROF,image,xcen,ycen,radius,/oplot,expand=nn'
    print,'e.g.> RADPROF,img,xc,yc,6'
    return
    endif

  if (n_elements(oplot1) eq 0) then oplot1=0
  if (n_elements(expand1) eq 0) then expand1=1
  if (n_elements(silent) eq 0) then silent=0
  if (n_elements(noplot) eq 0) then noplot=0
  if (n_elements(guiderplot) eq 0) then guiderplot=0

  xc=fix(xcen) & yc=fix(ycen) & siz=fix(radius+2)
  tmp=extrac(image,xc-siz,yc-siz,siz*2,siz*2)*1.0d


  expand1=fix(expand1)
  if (expand1 gt 1) then begin
    tmp2=interpolate(tmp,findgen(siz*2*expand1)/expand1, $
      findgen(siz*2*expand1)/expand1,/grid,/cubic)
;   tmp2(where(tmp2 lt min(tmp)))=min(tmp)
    dist_circle,mask,siz*2*expand1,(xcen-xc+siz)*expand1, $
      (ycen-yc+siz)*expand1
  endif else begin
    tmp2=tmp
    dist_circle,mask,siz*2,xcen-xc+siz,ycen-yc+siz
    endelse


  s=size(mask)
  masktmp=dblarr(s(1)*2,s(2))
  masktmp(0:s(1)-1,*)=mask/expand1
  masktmp(s(1):s(1)*2-1,*)=-mask/expand1

  fluxtmp=dblarr(s(1)*2,s(2))
  fluxtmp(0:s(1)-1,*)=tmp2
  fluxtmp(s(1):s(1)*2-1,*)=tmp2

  srt=sort(masktmp)
  xarr=masktmp(srt)
  flux=fluxtmp(srt)

  fit=gaussfit(xarr,flux,c)

  xscl=1.0
  if (noplot eq 0) then begin
    if (n_elements(miniplot) gt 1) then begin
      plot,mask/expand1,tmp2,psym=4,xr=[0,radius], $
        yr=miniplot(4:5),ysty=5,xsty=5,/noerase,pos=miniplot(0:3)
      endif else $
    if (oplot1 eq 0) then begin
      if (guiderplot eq 1) then begin
        xscl=0.42
        plot,mask/expand1*xscl,tmp2,psym=4,xr=[0,radius*xscl], $
          yr=[c(3)*.95,(c(3)+c(0))*1.05],ysty=1,/noerase, $
          xtit='Radius (arcsec)',ytit='DN',pos=[0.15,0.06,0.98,0.415], $
          xcharsize=0.7,ycharsize=0.7
        endif $
      else begin
        ymin=min(tmp2) & ymax=max(tmp2)
        yrng=ymax-ymin
        ymax=ymax+yrng*0.05
        plot,mask/expand1,tmp2,psym=4,xr=[0,radius], $
;          yr=[c(3)*.95,(c(3)+c(0))*1.05],ysty=1, $
          yr=[ymin,ymax],ysty=1, $
          xtit='Radius (Pixels)',ytit='Counts'
        endelse
      endif $
    else oplot,mask/expand1,tmp2,psym=2
    endif

  x=dindgen(radius*100)/100
  if (noplot eq 0) then $
    oplot,x*xscl,c(0)*exp(-x^2/(2*c(2)^2)) + c(3) + c(4)*x + c(5)*x^2, $
    color=!d.table_size-2

  gaussvol=c(0)*2*!dpi*c(2)^2
  gaussparams=c			; [amplitude,center,sigma,sky level,linear term,
				;  x^2 term]

  if silent then return

  print,'    X       Y         Sky        FWHM       Max           Total'
  print,'--------  ------  ----------  --------  ----------  -------------'

  print,format='(2f8.2,f12.3,f10.3,f12.2,f15.2)', $
    xcen,ycen,c(3),c(2)*2.35,c(0),gaussvol


  return

end

