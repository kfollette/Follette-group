
; This code was written in 2007 by Nick Konidaris.
; ds9__define is now released under the GNU public license.
;
;	Communicate with the author:  npk@astro.caltech.edu
;
;	ds9 uses the XPA library to communicate with the excellent
;		ds9 display program.  This is a replacement for ATV
;
;	Requirements:
;		1) A working XPA library, binaries must be in the shell's
;			PATH.
;		2) A working DS9 binary.  
;
;	Code
;		1) Checks to see if ds9 has been started
;		1b) If not started, spawns ds9 and loops for 15sec until
;			ds9 stars or timeouts.  Does not handle the case
;			of ds9 not starting!!  (BUG)
;		2) Stores the old shared-memory segment into a temporary
;			variable.  Then, allocates a new shared-memory
;			segment.
;		3) Via XPA, tells ds9 the id number of the shared mem
;			segment.
;		4) Unmaps the old shared memory segment.
;

PRO ds9::set_frame, ix
	self->cmd, 'frame frameno ' + string(ix)
END 

PRO ds9::update_wcs, hdr
	s = strcompress(string(randomu(s,1)),/rem)
	openw, 1, '/tmp/' + s[0] + 'wcs.tmp'
	for i = 0, n_elements(hdr)-1 do begin
		printf, 1, hdr[i]
	endfor
	close,1
	self->cmd, 'wcs append file /tmp/' + s[0] + 'wcs.tmp'
	file_delete, '/tmp/' + s[0] + 'wcs.tmp'
END

PRO ds9::new_frame
	;self->cmd, 'frame new'
END 

PRO ds9::cleanup
	catch, error_status
	if error_status ne 0 then begin
		catch, /cancel
		return
	endif
	;self->cmd, 'quit'
	if self.segname ne '' then $
		shmunmap, self.segname
END

PRO ds9::delete_all_frames
	self->cmd, 'frame delete all'
END 

PRO ds9::zoom_to_fit
	self->cmd, 'zoom to fit'
END

PRO ds9::line, x0, y0, x1, y1, color=color, text=text
	
	if self.os eq -1 then return
	if keyword_set(text) eq 0 then text = ""

	if keyword_set(color) eq 0 then color = 'blue'
	cmd = 'physical; line ' + string(x0) + ' ' + string(y0)  + ' ' + string(x1) + $
		string(y1) + ' # color = ' + color + " text = {" + text + "}"

	self->pipecmd, cmd, 'regions'
end

PRO ds9::ellipse, x0, y0, rx, ry, theta=theta, color=color
	
	if keyword_set(color) eq 0 then color = 'blue'

	angle = (keyword_set(theta) ? string(theta) : '0')
	cmd = 'physical; ellipse' + string(x0) + ' ' + string(y0)  + ' ' + string(rx) + ' ' + $
		string(ry) + ' ' + angle + ' # color = ' + color

	self->pipecmd, cmd, 'regions'
end

; Take a bunch of segments  connect
; Take a bunch of segments  connect
PRO ds9::lineconnected, x, y, color=color
	
	if (size(self.os))[1] eq 0 then return

	if n_elements(x) ne n_elements(y) then begin
		message, 'The length of x and y vectors is not the same'
	endif

	if keyword_set(color) eq 0 then color = 'blue'

	for i = 0, n_elements(x) - 2 do begin
		x0 = x[i] & y0 = y[i]
		x1 = x[i+1] & y1 = y[i+1]
		ds9line, x0, y0, x1, y1, color=color
	endfor
end


PRO ds9::pipecmd, pipe, cmd
	spawn, 'echo "' + pipe+ '" | ' + self.path + 'xpaset ' + self.title + ' ' + cmd
END

PRO ds9::cmd,cmd 
	spawn, self.path + 'xpaset -p ' + self.title + ' ' + cmd
END

;Issue an xpaget command
;added by Jared Males
PRO ds9::get, cmd, fn
        fn = '/tmp/ds9idl.out'
        fullcmd = self.path + 'xpaget ' + self.title + ' ' + cmd + ' > ' + fn
        spawn, fullcmd
        
END
        
PRO ds9::spawn_ds9, geometry=geometry, _extra=extra
	p = self.path
	spawn, p + 'xpaaccess ' + self.title, res

	if keyword_set(geometry) eq 0 then geometry = '' $
	else geometry = ' -geometry ' + geometry
	if res[0] eq 'no' then begin
		print, 'Spawning ds9 with title: ' + self.title
		spawn, p+'xpans&'
		; -- 1b
		spawn, p + 'ds9 -title ' + self.title + ' -port 0' + geometry + '&'
		for i = 0, 130 do begin
			wait, .5
			spawn, p + 'xpaaccess ' + self.title, res
			if res[0] eq 'yes' then break
		endfor
	endif
END



FUNCTION ds9::INIT, _EXTRA=_extra
	p = getenv('DS9')
	;if p eq '' then self.path = '/usr/local/bin/' $
	;else self.path = p
        self.path = ''
        
	self.os = -1
	self.segname = 'none'
	;Use just plain ole ds9, instead of timestamped.  this means we can use any already opened ds9
	self.title = 'ds9' ; + strcompress(string(systime(/sec),format='(i13)'),/rem)

	; -- 1
	self->spawn_ds9, _extra=_extra
	
	self.statsw = 10L
        self.plotw = 20
        self.fwhm = 2.5
        self.xcen = -1.
        self.ycen = -1.
   
	return,1
END

PRO ds9::rgbimage, $
	red, $
	green, $
	blue

	self->cmd, 'rgb'
	self->image, red
	self->cmd, 'rgb green'
	self->image, green 
	self->cmd, 'rgb blue'
	self->image, blue 
END

pro ds9::setimage, img
;Store the image in shared memory
   oldseg = '0'
   if self.os ne -1 then begin
      oldseg = self.segname
   endif
   shmmap, template=float(img), get_name=segname,/sysv, get_os_handle=os
   self.segname = segname
   self.os = os
   var = shmvar(self.segname)

   if((size(img))[0] eq 2) then    var[0,0] = float(img)
   if((size(img))[0] eq 3) then    var[0,0,0] = float(img)
   
   ;Now unmap the old segment, if needed
   if oldseg ne '0' then begin
      shmunmap, oldseg
   endif
        
end

pro ds9::getimage, img
;Get a mapped image array
   img = shmvar(self.segname, /float)
end

PRO ds9::dispimage, $
            img, $	; The image to display 
               preserve=preserve, $	; Preserve scale during load
                  newframe=newframe, $ 	; Keyword set if image to display in new frame
                     frame=frame ;Keyword to load in a specific frame (e.g. frame=2)
                     
   self->spawn_ds9

   p = self.path

   self->setimage, img

   sz = size(img,/dim)

   ; -- 3
	if (n_elements(frame) eq 1) then begin
		; kwd/kf edit to force frame 1 each time (mar 2019)
		spawn, p + 'xpaset -p ' + self.title + ' frame ' + 'delete'
		spawn, p + 'xpaset -p ' + self.title + ' frame ' + string(frame)
	endif
	
	if (n_elements(frame) ne 1 and keyword_set(newframe)) then begin
                ;spawn, p + 'xpaset -p ' + self.title + ' frame new'
        endif
	
	if keyword_set(preserve) then begin
		;spawn, p + 'xpaset -p ' + self.title + ' preserve scale yes'
	endif else begin
		;spawn, p + 'xpaset -p ' + self.title + ' preserve scale no'
	endelse
	
	if((size(img))[0] eq 2) then begin
	spawn, p + 'xpaset -p ' +  self.title + ' shm array shmid ' + strtrim(self.os,2) + ' "[xdim='+$
		strtrim(sz[0],2)+',ydim='+strtrim(sz[1],2)+',bitpix=-32]"&'
        endif
        
        if((size(img))[0] eq 3) then begin
        spawn, p + 'xpaset -p ' +  self.title + ' shm array shmid ' + strtrim(self.os,2) + ' "[xdim='+$
                strtrim(sz[0],2)+',ydim='+strtrim(sz[1],2)+',zdim='+strtrim(sz[2],2)+',bitpix=-32]"&'
        endif
        

END

PRO ds9::basic_imexam, key, x, y, z
;+
; PURPOSE:
;   Get the key name and coordinate of a key press
;
; AUTHOR:
;   Added by Jared Males, April 2013
;
;-
   self->spawn_ds9
   cmd = 'imexam key coordinate image'
   self->get, cmd, fn
   
   readcol, fn, key, sx, sy, format='A,I,I', /silent
   
   ;Now convert to integer, 0 counted 
   x = strn(sx) - 1
   y = strn(sy) - 1
   
   if(arg_present(z)) then begin
      ;get the cube coordinate
      cmd = 'cube'
      self->get, cmd, fn
      
      readcol, fn, sz, format='I', /silent
      
      z = strn(sz) - 1 ;0 count
      
   endif
   
END

PRO ds9::coordinate, x, y, z
;+
; PURPOSE:
;   Get the coordinate of a mouse click
;
; AUTHOR:
;   Added by Jared Males, April 2013
;
;-

   self->spawn_ds9
   cmd = 'imexam coordinate image'
   self->get, cmd, fn
   
   readcol, fn, sx, sy, format='I,I', /silent
   
   ;Now convert to integer, 0 counted 
   x = strn(sx) - 1
   y = strn(sy) - 1
   
   if(arg_present(z)) then begin
      cmd = 'cube'
      self->get, cmd, fn
      
      readcol, fn, sz, format='I', /silent
      
      z = strn(sz) - 1 ;0 count
      
   endif
   
END

PRO ds9::data, xc, yc, data, x, y, w, h
;+
; PURPOSE:
;   Get data from ds9
;
; AUTHOR:
;   Added by Jared Males, April 2013
;
;-

   self->spawn_ds9
   cmd = 'data image ' + string(x) + ' ' + string(y) + ' ' + string(w) + ' ' + string(h) + ' no'
   self->get, cmd, fn
   
   readcol, fn, xc, yc, eqs, data, format='I,I,A,D', /silent
END

PRO ds9::dumpfile, fname
;+
; PURPOSE:
;   Dumps current frame to a fits file
; AUTHOR:
;   Jared Males
;-

   self->spawn_ds9
   cmd = 'fits'
   self->get, cmd, fname
END

pro ds9::set_statsw, sw
   self.statsw = sw
end

function ds9::statsw
   return, self.statsw
end

pro ds9::set_plotw, pw
   self.plotw = pw
end

function ds9::plotw
   return, self.plotw
end

pro ds9::set_fwhm, fw
   self.fwhm = fw
end

function ds9::fwhm
   return, self.fwhm
end

pro ds9::set_xcen, cen
   self.xcen = cen
end

function ds9::xcen
   return, self.xcen
end

pro ds9::set_ycen, cen
   self.ycen = cen
end

function ds9::ycen
   return, self.ycen
end

PRO ds9__define
	params = {DS9, $
		title: '', $
		os: -1L, $
		segname: '', $
		path: '', $
		statsw: 10L, $  ;the width of stats regions
		plotw:  20L, $  ;the width of imexam plots
		fwhm:   3.d, $  ;the full-width at half-maximum to assume for centroiding, etc.
		xcen:   -1.d,$  ;for imexam, the xcenter of the image.  if -1, then 0.5*(dim1-1)
		ycen:   -1.d $  ;for imexam, the ycenter of the image.  if -1, then 0.5*(dim2-1)
	}
END
