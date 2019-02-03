pro mxloadct, ctname, ncol=ncol


n = strpos(ctname, ':')

type = strmid(ctname, 0, n)
name = strmid(ctname, n+1)

if(type eq 'idl') then begin

   loadct, long(name), ncolors=ncol, /silent

   return
endif

if(type eq 'sao') then begin

   saoct, name

   return

endif

if(type eq 'matlab') then begin

   matlabct, name

   return
endif

end


