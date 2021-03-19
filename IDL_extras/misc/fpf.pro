pro fpf, im, res, loc

;+
;
;res = size of resolution element in pixels (lam/D for diffraction limited, FWHM for not)
;loc = 2 element vector with x and y coordinates of planet center
;
;-

npix=(size(im))[1]
xcen=(npix-1)/2.
ycen=xcen

r=dblarr(npix,npix)
rres=dblarr(npix,npix)

r=radmap(im, xcen, ycen)
rres=r/res

nsamples=round(2*!pi*r[loc(0),loc(1)]/res)
stamps=dblarr(res, res, nsamples)
cen=dblarr(nsamples,2)
cen[0,*]=loc
thetastep=2*!pi/nsamples
thetas=dblarr(nsamples)
thetas[0]=atan((loc[1]-ycen)/(loc[0]-xcen))
rad=r[loc[0],loc[1]]
inten=dblarr(nsamples)

for i=0, nsamples-1 do begin
  if i gt 0 then begin
    thetas[i]=thetas[i-1]+thetastep
    cen[i,0]=round(cos(thetas[i])*rad+xcen)
    cen[i,1]=round(sin(thetas[i])*rad+ycen)
  endif
  stamphldr=im[cen[i,0]-(res-1)/2:cen[i,0]+(res-1)/2,cen[i,1]-(res-1)/2:cen[i,1]+(res-1)/2]
  stamphldr(where(radmap(stamphldr) gt res/2.))='NaN'
  stamps[*,*,i]=stamphldr
  inten[i]=total(stamps[*,*,i], /nan)
endfor

x1=inten[0]
x2=mean(inten[1:nsamples-1])
s2=stdev(inten[1:nsamples-1])
tau=x1/sqrt(1+1/(nsamples-1))
nu=nsamples-1

;fpf=qromo('tdist', tau, /midexp)

stop

end