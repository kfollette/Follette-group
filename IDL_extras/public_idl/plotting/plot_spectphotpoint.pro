pro plot_spectphotpoint, filtlam, filtwidth, spectflux, measflux, measerr, $
    lthick=lthick, sthick=sthick, snum=snum, ssz=ssz, lcol=lcol, scol=scol, nohat=nohat

if (n_elements(lthick) eq 0) then lthick = 1
if (n_elements(sthick) eq 0) then sthick = 1
if (n_elements(snum) eq 0) then snum = 1
if (n_elements(ssz) eq 0) then ssz = 1
if (n_elements(scol) eq 0) then scol = 0
if (n_elements(lcol) eq 0) then lcol = 0

oplot, [filtlam-.5*filtwidth, filtlam+.5*filtwidth], [spectflux, spectflux], thick=lthick, color=lcol

oploterror, filtlam, measflux, measerr*measflux, psym=snum, symsize=ssz, thick=sthick, errthick=sthick, color=scol, errcol=scol, nohat=keyword_set(nohat)

end





