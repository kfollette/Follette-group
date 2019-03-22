; This code was written in 2007 by Nick Konidaris.
; ds9__define is now released under the GNU public license.
;
;   Communicate with the author:  npk@astro.caltech.edu
;
;   ds9 uses the XPA library to communicate with the excellent
;       ds9 display program.  This is a replacement for ATV
;
;   Requirements:
;       1) A working XPA library, binaries must be in the shell's
;           PATH.
;       2) A working DS9 binary.  
;
;   Code
;       1) Checks to see if ds9 has been started
;       1b) If not started, spawns ds9 and loops for 15sec until
;           ds9 stars or timeouts.  Does not handle the case
;           of ds9 not starting!!  (BUG)
;       2) Stores the old shared-memory segment into a temporary
;           variable.  Then, allocates a new shared-memory
;           segment.
;       3) Via XPA, tells ds9 the id number of the shared mem
;           segment.
;       4) Unmaps the old shared memory segment.
;
pro ds9::set_frame, ix
  self->cmd, 'frame frameno ' + string(ix)
end 
pro ds9::new_frame
  self->cmd, 'frame new'
end 
pro ds9::delete_all_frames
  self->cmd, 'frame delete all'
end 
pro ds9::update_wcs, hdr
  s = strcompress(string(randomu(s,1)),/rem)
  
  openw,1,'/tmp/'+s[0]+'wcs.tmp'
  for i = 0,n_elements(hdr)-1 do printf,1,hdr[i]
  close,1
  
  self->cmd, 'wcs append /tmp/' + s[0] + 'wcs.tmp'
  file_delete, '/tmp/' + s[0] + 'wcs.tmp'
end
pro ds9::zoom_to_fit
  self->cmd, 'zoom to fit'
end
pro ds9::zoom,val
  self->cmd, 'zoom to '+string(val)
end
pro ds9::line, x0, y0, x1, y1, color=color, text=text, fixed=fixed, $
               arrowb=arrowb,arrowe=arrowe,width=width
  if self.os eq -1 then return
  if keyword_set(text) eq 0 then text = ""
  if keyword_set(width) eq 0 then width = 1
  if keyword_set(fixed) eq 0 then move = 1 else move = 0
  if keyword_set(arrowb) eq 0 then arrowb = 0 else arrowb = 1
  if keyword_set(arrowe) eq 0 then arrowe = 0 else arrowe = 1
  ;; different coordinate system IDL/ds9
  xx0 = x0+1
  yy0 = y0+1
  xx1 = x1+1
  yy1 = y1+1
  
  if keyword_set(color) eq 0 then color = 'green'
  cmd = 'physical; line '+string(xx0)+' '+string(yy0)+' '+string(xx1)+ $
        string(yy1)+' # color = '+color+' width = ' + string(fix(width))+ $
        " text = {"+text+"}"+' '+'move='+string(move)+ $
        ' line='+string(arrowb)+' '+string(arrowe)
  self->pipecmd, cmd, 'regions'
end
pro ds9::circle, x0, y0, r, color=color, width=width, fixed=fixed
  if keyword_set(color) eq 0 then color = 'green'
  if keyword_set(width) eq 0 then width = 1
  if keyword_set(fixed) eq 0 then move = 1 else move = 0
  
  ;; different coordinate system IDL/ds9
  xx0 = x0+1
  yy0 = y0+1
  
  angle = (keyword_set(theta) ? string(theta) : '0')
  cmd = 'physical; circle' + string(xx0) + ' ' + string(yy0) + ' ' + string(r) + ' ' + $
        ' # color = ' + color + ' ' + 'width = ' + string(fix(width))+' '+'move='+string(move)
  
  self->pipecmd, cmd, 'regions'
end
pro ds9::ellipse, x0, y0, rx, ry, theta=theta, color=color, width=width , fixed=fixed
  if keyword_set(color) eq 0 then color = 'green'
  if keyword_set(fixed) eq 0 then move = 1 else move = 0
  
  ;; different coordinate system IDL/ds9
  xx0 = x0+1
  yy0 = y0+1  
  
  angle = (keyword_set(theta) ? string(theta) : '0')
  cmd = 'physical; ellipse' + string(xx0) + ' ' + string(yy0)  + ' ' + string(rx) + ' ' + $
        string(ry) + ' ' + angle + ' # color = ' + color+ ' ' + $
        'width = ' + string(fix(width))+' '+'move='+string(move)
  
  self->pipecmd, cmd, 'regions'
end
pro ds9::box, x0, y0, w, h, ang, color=color, width=width, fixed=fixed
  if keyword_set(color) eq 0 then color = 'green'
  if keyword_set(width) eq 0 then width = 1
  if keyword_set(fixed) eq 0 then move = 1 else move = 0
  
  ;; different coordinate system IDL/ds9
  xx0 = x0+1
  yy0 = y0+1  
  angle = (keyword_set(theta) ? string(theta) : '0')
  cmd = 'physical; box' + string(xx0) + ' ' + string(yy0) + ' ' + string(w) + ' ' + $
        string(h) + ' ' + string(ang)+ ' # color = ' + color + ' ' + 'width = ' + $
        string(fix(width))+' '+'move='+string(move)
  
  self->pipecmd, cmd, 'regions'
end
pro ds9::text, x0, y0, text, color=color, fixed=fixed
  if keyword_set(color) eq 0 then color = 'green'
  if keyword_set(width) eq 0 then width = 1
  if keyword_set(fixed) eq 0 then move = 1 else move = 0
  
  ;; different coordinate system IDL/ds9
  xx0 = x0+1
  yy0 = y0+1
  angle = (keyword_set(theta) ? string(theta) : '0')
  cmd = 'physical; text' + string(xx0) + ' ' + string(yy0) + ' ' + $
        ' # text={' + text + '}'+ ' '+'color = ' + color + ' ' + $
        'move='+string(move)
  
  self->pipecmd, cmd, 'regions'
end
pro ds9::point, x0, y0, style=style, color=color, width=width, fixed=fixed
  if keyword_set(color) eq 0 then color = 'green'
  if keyword_set(width) eq 0 then width = 1
  if keyword_set(fixed) eq 0 then move = 1 else move = 0
  if not keyword_set(style) then style = 'cross'
  
  ;; different coordinate system IDL/ds9
  xx0 = x0+1
  yy0 = y0+1
  
  angle = (keyword_set(theta) ? string(theta) : '0')
  cmd = 'physical; point' + string(xx0) + ' ' + string(yy0) + ' ' + $
        ' # point = '+style+' color = ' + color + ' ' + 'width = ' + $
        string(fix(width))+' '+'move='+string(move)
  
  self->pipecmd, cmd, 'regions'
end
pro ds9::imexam,x,y,wcs=wcs
  if (Keyword_set(wcs) eq 1) then coord = 'wcs fk5 degrees' $
  else coord = 'image'
  cmd = 'imexam coordinate '+coord
  
  spawn,'xpaget '+self.title+' '+cmd,output
  if (output ne ' ') then begin
     sep = strsplit(output,' ',/extract)
     ;; different coordinate system IDL/ds9
     x = float(sep[0])-1
     y = float(sep[1])-1
  endif else begin
     x = !values.f_nan
     y = !values.f_nan
  endelse
end
pro ds9::centroid,im,par,extend=extend,plot=plot,_extra=_extra
  if not keyword_set(extend) then extend = 20
  dim = 2*extend+1
  
  self->imexam,x,y
  if ((finite(x) ne 1) or (finite(y) ne 1)) then par = -1.
  
  cx = round(x)
  cy = round(y)
  sub = im[cx-extend:cx+extend,cy-extend:cy+extend]
  res = mpfit2dpeak(sub,par,_extra=_extra)
  if keyword_set(plot) then begin
     window,/free
     nx = par[4]
     ny = par[5]
     dx = nx-floor(nx)
     dy = ny-floor(ny)
     
     xx = (fltarr(dim)+1) ## findgen(dim)
     yy = findgen(dim) ## (fltarr(dim)+1)
     
     nsub = interpolate(sub,xx+dx,yy+dy,cubic=-0.5,missing=0.0)
     
     plot,findgen(dim)-extend-dx,nsub[floor(nx),*],xs=1,xr=[-extend,extend], $
          title='Centroid',/nodata
     plots,0,!y.crange
     oplot,findgen(dim)-extend-dx,nsub[floor(nx),*],linestyle=0
     oplot,findgen(dim)-extend-dy,nsub[*,floor(ny)],linestyle=1
     legend,['x','y'],linestyle=[0,1],/top,/left
  endif
end
pro ds9::scale,linear=linear,sqrt=sqrt,log=log,histeq=histeq,minmax=minmax,zscale=zscale, $
               limits=limits
  flag_scale = keyword_set(linear)+keyword_set(log)+keyword_set(sqrt)+keyword_set(histeq)  
  flag_mode  = keyword_set(minmax)+keyword_set(zscale)+keyword_set(limits)
  if (flag_scale eq 0) then begin
     linear = 1
  endif else if (flag_scale gt 1) then begin
     print,'You must choose only one scale: linear, sqrt, log, histeq!'
     return
  endif
  if (flag_mode eq 0) then begin
     minmax = 1
  endif else if (flag_mode gt 1) then begin
     print,'You must choose only one mode: minmax, zscale, limits!'
     return
  endif
  if ((keyword_set(limits) eq 1) and (n_elements(limits) ne 2)) then begin
     print,'Limits must have 2 elements!'
     return
  endif
  
  if keyword_set(linear) then self->cmd,'scale linear'
  if keyword_set(sqrt)   then self->cmd,'scale sqrt'
  if keyword_set(log)    then self->cmd,'scale log 100'  
  if keyword_set(histeq) then self->cmd,'scale histeq'
  if keyword_set(limits) then self->cmd,'scale limits '+string(limits[0])+' '+string(limits[1])
  if keyword_set(minmax) then self->cmd,'scale mode minmax'
  if keyword_set(zscale) then self->cmd,'scale mode zscale'
end
pro ds9::pipecmd, pipe, cmd
  spawn, 'echo "' + pipe+ '" | ' + self.path + 'xpaset ' + self.title + ' ' + cmd
end
pro ds9::cmd,cmd 
  spawn, self.path + 'xpaset -p ' + self.title + ' ' + cmd
end
pro ds9::spawn_ds9,geometry=geometry,_extra=extra
  p = self.path
  spawn, p + 'xpaaccess ' + self.title, res
  
  if keyword_set(geometry) eq 0 then geometry = '' $
  else geometry = ' -geometry ' + geometry
  if res[0] eq 'no' then begin
     print, 'Spawning ds9 with title: ' + self.title
     spawn, p+'xpans&'
     ;; -- 1b
     spawn, p + 'ds9 -title ' + self.title + ' -port 0' + geometry + '&'
     for i = 0, 130 do begin
        print,'waiting...'
	      wait, .5
        spawn, p + 'xpaaccess ' + self.title, res
        if res[0] eq 'yes' then begin
          print,'it worked!'
        endif
        if res[0] eq 'yes' then break
     endfor
  endif
end
pro ds9::cleanup
  catch, error_status
  if error_status ne 0 then begin
     catch, /cancel
     return
  endif
  self->cmd, 'quit'
  if self.segname ne '' then $
     shmunmap, self.segname
end
FUNCTION ds9::init, _EXTRA=_extra
  p = getenv('DS9')
  $if p eq '' then self.path = '~/local/xpaa/bin/' $ original
  if p eq '' then self.path = '/opt/local/bin/' $ kdw edit
  else self.path = p
  
  self.os = -1
  self.segname = 'none'
  self.title = 'ds9_' + strcompress(string(systime(/sec),format='(i13)'),/rem)
  
  ;; -- 1
  self->spawn_ds9, _extra=_extra
  return,1
end
pro ds9::rgbim, $
   red, $
   green, $
   blue
  
  self->cmd,'rgb'
  self->im,red
  self->cmd,'rgb green'
  self->im,green 
  self->cmd,'rgb blue'
  self->im,blue 
end
pro ds9::im, $
   img, $                       ; The image to display
   preserve=preserve, $         ; Preserve scale during load
   frame=frame, $               ; Keyword set if image to display in new frame
   hdr=hdr                      ; FITS header to update WCS informations
  
  self->spawn_ds9
  
  p = self.path
  ;; -- 2
  oldseg = '0'
  if self.os ne -1 then begin
     oldseg = self.segname
  endif
  shmmap, template=double(img), get_name=segname,/sysv, get_os_handle=os
  self.segname = segname
  self.os = os
  var = shmvar(self.segname)
  
  sz = size(img,/dim)
  if (n_elements(sz) eq 2) then var[0,0] = double(img) $
  else if (n_elements(sz) eq 3) then var[0,0,0] = double(img) $
  else begin
     print,'ds9 error: cannot handle data!'
     goto,fin
  endelse
  
  ;; -- 3
  if keyword_set(frame) then begin
     spawn, p + 'xpaset -p ' + self.title + ' frame new'
  endif
  ;if keyword_set(preserve) then begin
  ;   spawn, p + 'xpaset -p ' + self.title + ' preserve scale yes'
  ;endif else begin
  ;   spawn, p + 'xpaset -p ' + self.title + ' preserve scale no'
  ;endelse
  cube = (n_elements(sz) eq 3) ? ',zdim='+strtrim(sz[2],2) : ''  
  spawn, p + 'xpaset -p ' +  self.title + ' shm array shmid ' + strtrim(self.os,2) + ' "[xdim='+$
         strtrim(sz[0],2)+',ydim='+strtrim(sz[1],2)+cube+',bitpix=-64]"&'
  if keyword_set(hdr) then self->update_wcs,hdr
fin:
  ;; -- 4
  if oldseg ne '0' then begin
     shmunmap, oldseg
  endif  
end
pro ds9__define
  params = {DS9, $
            title: '', $
            os: -1L, $
            segname: '', $
            path: '' $
           }
end
pro ds9_slider_event,event
  ;; get GUI structure and data from top base
  widget_control,event.top,get_uvalue=state
  steps = ['10.0',' 5.0',' 1.0',' 0.5',' 0.1']
  
  widget_control,event.id,get_uvalue=uvalue
  case uvalue of
     'min':  state.data.min = event.value
     'max':  state.data.max = event.value
     'step': state.data.step = steps[event.index]
     'quit': begin
        !v->scale,/minmax
        widget_control,event.top,/destroy
        return
     end
     else:
  endcase   
  min     = state.data.min
  max     = state.data.max
  step    = state.data.step
  nlevels = (max-min)/step  
  
  if (uvalue eq 'slider') then begin
     state.data.slider = event.value*step+min
     !v->scale,limits=[state.data.slider,state.data.slider]
  end     
  slider = state.data.slider
  
  widget_control,state.id.slider,set_slider_min=min
  widget_control,state.id.slider,set_slider_max=min+nlevels
  widget_control,state.id.slider,set_value=(slider-min)/step
  widget_control,state.id.value,set_value=string(state.data.slider,format='(D8.2)')
  
  widget_control,event.top,set_uvalue=state
end
pro ds9::slider
  id = { $
       base   : 0, $
       min    : 0, $
       max    : 0, $
       step   : 0, $
       slider : 0, $
       value  : 0, $
       quit   : 0  $
       }
  data = { $
         min    : 0.,  $
         max    : 10., $
         step   : 1.0, $
         slider : 5.0  $
         }
  
  state = { $
          id   : id,  $
          data : data $
          }
  
  ;; root base
  root_base_id = widget_base(title='Slider',/col)
  
  widget_control,root_base_id,/managed
  state.id.base = widget_base(root_base_id,row=3)
  dummy = widget_base(state.id.base,col=3)
  state.id.min    = cw_field(dummy,title='Min: ',/row, $
                             value=state.data.min,/integer,uvalue='min', $
                             /all_events,/column,xsize=12)
  state.id.max    = cw_field(dummy,title='Max: ',/row, $
                             value=state.data.min,/integer,uvalue='max', $
                             /all_events,/column,xsize=12)
  steps = ['10.0',' 5.0',' 1.0',' 0.5',' 0.1']
  state.id.step   = widget_droplist(dummy,title='Step:',/flat,value=steps,uvalue='step')
  
  dummy = widget_base(state.id.base,col=2)
  state.id.slider = widget_slider(dummy,/drag, $
                                  value=state.data.slider,uvalue='slider', $ $
                                  minimum=state.data.min,maximum=state.data.max, $
                                  xsize=200,/suppress_value)
  state.id.value  = widget_label(dummy,value=string(state.data.slider,format='(D8.2)'),uvalue='value', $
                                 /align_left)
  dummy = widget_base(state.id.base)
  state.id.quit   = widget_button(dummy,/align_right,uvalue='quit',value=' Quit ')
  
  ;; update values
  widget_control,state.id.min,set_value=state.data.min
  widget_control,state.id.max,set_value=state.data.max
  widget_control,state.id.step,set_droplist_select=where(steps eq state.data.step)
  widget_control,state.id.slider,set_value=state.data.slider
  widget_control,state.id.slider,set_value=state.data.slider
  
  ;; save GUI structure and data in the top base uvalue
  widget_control,root_base_id,set_uvalue=state
  ;; draw GUI
  widget_control,root_base_id,/realize
  
  xmanager,'ds9_slider',root_base_id,/no_block
end
;;
;; instantiate a global ds9 object
;;
pro ds9
  defsysv,'!v',EXISTS = exists  
  if not exists then begin
     v = obj_new('ds9')
     defsysv, "!v",v,1
  endif else !v->spawn_ds9
end
