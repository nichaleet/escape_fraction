; =================
pro viewspec_event, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=value
    if (n_elements(value) eq 0) then value = ''
    name = strmid(tag_names(ev, /structure_name), 7, 4)

    case (name) of
        'BUTT': obj->handle_button, ev
        'TEXT': obj->handle_text, ev
        'TRAC': obj->handle_tracking, ev
        'COMB': obj->handle_combobox, ev
        'DRAW': begin
            if ev.type eq 0 then obj->handle_draw_click, ev
            if ev.type eq 2 then obj->handle_motion, ev
        end
        'DONE': widget_control, ev.top, /destroy
        else: begin
            widget_control, widget_info(ev.top, find_by_uname='spec'), /input_focus
        end
    endcase
end


; ================ BUTTONS ================
function line_event, ev
    self->redraw
end


pro viewspec::handle_button, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=uvalue
    
    case (uvalue) of
        'default_range': begin
            self->default_range
            self->redraw
        end
        'blue': self->lambdarange, /blue
        'red': self->lambdarange, /red
        else:
    endcase
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


pro viewspec::skyylim
    wsky = where(((*self.science).skyfit)[*,0,0] ne 0, csky)
    if csky gt 0 then begin
        skyfit = (*self.science).skyfit[wsky,*]
        yminmax = (*self.science).dlam
        if skyfit[0,0] ne -1 then yminmax = [yminmax, skyfit[wsky,1]+skyfit[wsky,2], skyfit[wsky,1]-skyfit[wsky,2]]    
        self.skyylim = minmax(yminmax*2.35)    
    endif else begin
        self.skyylim = [0.9, 1.7]
    endelse
 end


pro viewspec::newspec, noredraw=noredraw

    self->statusbox
    self.ylim = minmax((*self.science).spec,/nan)
    self->skyylim
    if ~keyword_set(noredraw) then self->draw
end


pro viewspec::default_range, update=update
    if ~keyword_set(update) then begin
        self.ylim = minmax((*self.science).spec,/nan)
        self.ylim[1] *= 1.1
        self.divylim = [-1.0, 2.5]
        self->skyylim
        self.lambdalim = (minmax((*self.science).lambda / (1d + (*self.science).zspec),/nan) < 1700) > 1100
        self.lambdalim[0] >= 1100.
        self.lambdalim[1] <= 1700.; / (1d + (*self.science).zspec)
     endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(g8.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(g8.2)'), /rem)
        end
        mode eq 2: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.skyylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.skyylim[1], format='(D5.2)'), /rem)
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.divylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.divylim[1], format='(D5.2)'), /rem)
        end
    endcase
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 3 ? (*self.science).zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
end

pro viewspec::lambdarange, red=red, blue=blue
    if keyword_set(red)+keyword_set(blue) ne 1 then message, 'You must specify red or blue.'
    lrange = self.lambdalim[1] - self.lambdalim[0]
    if keyword_set(blue) then begin
        if self.lambdalim[0] lt min((*self.science).lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is bluest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] - 0.6*lrange
        lhighnew = self.lambdalim[1] - 0.6*lrange
    endif
    if keyword_set(red) then begin
        if self.lambdalim[1] gt max((*self.science).lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is reddest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] + 0.6*lrange
        lhighnew = self.lambdalim[1] + 0.6*lrange
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 3 ? (*self.science).z : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
    self->lambdalow, llownew, /noredraw
    self->lambdahigh, lhighnew
end


; ============== TEXT BOXES =============
pro viewspec::handle_text, ev
    widget_control, ev.id, get_uvalue=uvalue
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_value=val

    case (uvalue) of
        'ylow': self->ylow, val
        'yhigh': self->yhigh, val
        'lambdalow': self->lambdalow, val
        'lambdahigh': self->lambdahigh, val
        else: 
    end
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

pro viewspec::lambdalow, lambdalow, noredraw=noredraw
    self.lambdalim[0] = lambdalow
    if ~keyword_set(noredraw) then self->draw
end


pro viewspec::lambdahigh, lambdahigh, noredraw=noredraw
    self.lambdalim[1] = lambdahigh
    if ~keyword_set(noredraw) then self->draw
end


pro viewspec::ylow, ylow, noredraw=noredraw
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: self.ylim[0] = ylow
        mode eq 2: self.skyylim[0] = ylow
        else: self.divylim[0] = ylow
    endcase
    if ~keyword_set(noredraw) then self->draw
end


pro viewspec::yhigh, yhigh
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: self.ylim[1] = yhigh
        mode eq 2: self.skyylim[1] = yhigh
        else: self.divylim[1] = yhigh
    endcase
    self->draw
end



; ============== COMBOBOX  =============
pro viewspec::handle_combobox, ev
    ;widget_control, ev.id, get_uvalue=uvalue
    ;widget_control, ev.top, get_uvalue=obj
    self->statusbox
    self.ylim = minmax((*self.science).spec,/nan)
    self->skyylim
    self.keystate = 0
    self->draw
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


; ============= DRAW CLICK =============
pro viewspec::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    if mode lt 3 then click_coords /= 1d + (*self.science).zspec
    case ev.modifiers of
        0: begin
            case ev.press of
                1: lrange = (self.lambdalim[1] - self.lambdalim[0]) / 2.
                4: lrange = (self.lambdalim[1] - self.lambdalim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.
            zl = mode lt 3 ? (*self.science).zspec : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        1: begin
            if ev.press ne 1 then begin
                widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                return
            endif
            lrange = self.lambdalim[1] - self.lambdalim[0]
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.        
            zl = mode lt 3 ? (*self.science).zspec : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin
            widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
            case mode of
                mode le 0: ylim = self.ylim
                mode eq 2: ylim = self.skyylim
                else: ylim = self.divylim
            endcase
            case ev.press of
                1: yrange = (ylim[1] - ylim[0]) / 2.
                4: yrange = (ylim[1] - ylim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            ylownew = click_coords[1] - yrange/2.
            yhighnew = click_coords[1] + yrange/2.
            case 1 of
                mode le 0: begin
                    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(g8.2)'), /rem)
                    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(g8.1)'), /rem)
                end
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(D5.2)'), /rem)
                    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(D5.1)'), /rem)
                end
            endcase
            self->ylow, ylownew, /noredraw
            self->yhigh, yhighnew
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
            return            
        end
     endcase
 end
; ============== REDRAW ===============
pro viewspec::draw

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Drawing ...'
    self->default_range, /update
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    widget_control, widget_info(self.base, find_by_uname='zmode'), get_value=zmode
    widget_control, widget_info(self.base, find_by_uname='spec'), get_value=index
    wset, index
    science = *self.science
    case zmode of
       0: z=science.zspec
       1: z=science.zsys_lya
       2: z=science.zsys_ism_ave
    endcase

    case mode of
        0: begin                ;continuum fit
           if n_elements(size(science.contmask,/dimensions)) gt 1 then contmask = science.contmask[*,0]
           t = round(-1*ts_diff(science.contmask, 1))
           wstart = where(t eq 1, cstart)+1
           wend = where(t eq -1, cend)
            if science.contmask[0] eq 1 then begin
                if cstart eq 0 then begin
                    wstart = 0
                endif else begin
                    wstart = [0, wstart]
                endelse
                cstart += 1
            endif
            if science.contmask[n_elements(t)-1] eq 1 then begin
                if cend eq 0 then begin
                    wend = n_elements(t)-1
                endif else begin
                    wend = [wend, n_elements(t)-1]
                endelse
                cend += 1
            endif
            if cstart ne cend then message, 'There are a different number of starting and ending continuum wavelengths.'
            if cstart eq 0 or cend eq 0 then message, 'There are no continuum regions.'

            plot, science.lambda, science.spec, xrange=self.lambdalim * (1d + z), yrange=self.ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
            for i=0,cstart-1 do begin
                x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]] > (self.lambdalim[0] * (1d + z))) < (self.lambdalim[1] * (1d + z))
                y = [self.ylim[0], self.ylim[1], self.ylim[1], self.ylim[0]]
                polyfill, x, y, color=fsc_color('light cyan')
            endfor
            oplot, science.lambda, science.spec, color=fsc_color('black')
            oplot, science.lambda, science.cont, color=fsc_color('green')
            plot, science.lambda, science.spec, xrange=self.lambdalim * (1d + z), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3', /nodata, /noerase
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] * (1d + z) le !X.CRANGE[0] or (*self.linewaves)[i] * (1d + z) ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]] * (1d + z), [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, ((*self.linewaves)[i] * (1d + z))+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
        end

        1: begin                ;continuum division
            plot, science.lambda, science.contdiv, xrange=self.lambdalim * (1d + z), yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata
            oplot, self.lambdalim, [1.0, 1.0], color=fsc_color('orange')
            oplot, science.lambda, science.contdiv, color=fsc_color('black')
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
        end

        2: begin        ;sky line fit
            wsky = where(((*self.science).skylinemask)[*,0] ne -1, csky)
            plot, science.lambda, 2.35*science.dlam, xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6sky line FWHM (!sA!r!u!9 %!6!n)!3', xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xrange=self.lambdalim * (1d + z), yrange=self.skyylim
            if csky gt 0 then begin
                fit = science.skyfit[wsky,*,0]
                if fit[0,0] ne -1 then begin
                    oploterror, fit[wsky,0], 2.35*fit[wsky,1], 2.35*fit[wsky,2], psym=1, color=fsc_color('black'), errcolor=fsc_color('black'), /nohat
                    wdel = where(science.skylinemask[wsky] eq 0, cdel)
                    if cdel gt 0 then plots, fit[wdel,0], 2.35*fit[wdel,1], color=fsc_color('red'), psym=7
                endif
             endif
        end

        3: begin        ;rest frame
            t = round(-1*ts_diff(science.fitmask, 1))
            wstart = where(t eq 1, cstart)+1
            wend = where(t eq -1, cend)
            if science.fitmask[0] eq 1 then begin
                if cstart eq 0 then begin
                    wstart = 0
                endif else begin
                    wstart = [0, wstart]
                endelse
                cstart += 1
            endif
            if science.fitmask[n_elements(t)-1] eq 1 then begin
                if cend eq 0 then begin
                    wend = n_elements(t)-1
                endif else begin
                    wend = [wend, n_elements(t)-1]
                endelse
                cend += 1
            endif
            if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
            
            plot, science.lambda/(1d + z), science.contdiv, xrange=self.lambdalim, yrange=self.divylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
           
            if cstart eq 0 or cend eq 0 then message, 'There are no fitmask regions.', /info else begin
                for i=0,cstart-1 do begin
                    x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]]/(1d + z) > self.lambdalim[0]) < self.lambdalim[1]
                    y = [self.divylim[0], self.divylim[1], self.divylim[1], self.divylim[0]]
                    polyfill, x, y, color=fsc_color('light cyan')
                endfor
            endelse
            oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
            oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')
            if z ne 0 and finite(z) then oplot, science.lambda/(1d + z), science.contdiv, color=fsc_color('black') else oplot, science.lambda/(1d + science.zspec), science.contdiv, color=fsc_color('black')

            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]] / (1d + z), 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] le !X.CRANGE[0] or (*self.linewaves)[i] ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, (*self.linewaves)[i]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
             endfor
            plot, science.lambda/(1d + z), science.contdiv, xrange=self.lambdalim, yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), /nodata,/noerase

         end
    endcase        
    zl = mode lt 3 ? (*self.science).zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
 end


pro viewspec::statusbox, science=science
    if ~keyword_set(science) then science = *self.science
    
    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(science.objname, 2)
    widget_control, widget_info(self.base, find_by_uname='curz'), set_value=strcompress(string(science.zspec, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzlya'), set_value=strcompress(string(science.zsys_lya, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzism'), set_value=strcompress(string(science.zsys_ism_ave, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=science.sn gt 0 ? strcompress(string(science.sn, format='(D10.1)'), /rem) : unknown
    if n_elements(science.exptime) gt 2 then totexptime = (total(science.exptime,2))[0]/3600. else totexptime = (science.exptime)[0]/3600.
    widget_control, widget_info(self.base, find_by_uname='expt'), set_value=strcompress(string(totexptime, format='(D5.3)'), /rem)

 end


pro viewspec::getscience
  common filename,sciencefile
  
  widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'

  science = mrdfits(sciencefile, 1, /silent)
  widget_control, widget_info(self.base, find_by_uname='filename'), set_value=sciencefile
  self.science = ptr_new(science)
  widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3

end

pro viewspec::initialize_directory
    common filename,sciencefile

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading file ...'

    self.sciencefile = sciencefile

    self->getscience
    science = *self.science
    self->statusbox, science=science
    self->default_range
    self.lambdalim = [1100, 1700]
    self->newspec
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus    
end


; =============== INIT ================
function viewspec::INIT
    common filename,sciencefile

    base = widget_base(/row, title='ViewSpec', uvalue=self, mbar=menu, tab_mode=0, units=1)
    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base
    ; ------ LEFT -------
    wplotmode = widget_base(wleft, /column, /align_center, /frame)
    wplotradio = cw_bgroup(wplotmode, ['continuum fit', 'continuum division', 'sky line fit', 'rest frame'], /column, /exclusive, set_value=3, uname='mode', uvalue='mode', /no_release)
    wzmode = widget_base(wleft, /column, /align_center, /frame)
    wzradio = cw_bgroup(wzmode, ['z input', 'z Lya', 'z ism'], /column, /exclusive, set_value=2, uname='zmode', uvalue='zmode', /no_release)
    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    widlabel = widget_label(widbase, value='object name:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=180)
    wsnbase = widget_base(wcurobj, /align_center, /row)
    wsnlabel = widget_label(wsnbase, value='SNR = ', /align_right, uname='snlabel', xsize=95)
    wcursn = widget_label(wsnbase, value='     ', /align_left, uname='cursn', uvalue='cursn', xsize=150)
    wexpbase = widget_base(wcurobj, /align_center, /row)
    wexplabel = widget_label(wexpbase, value='Exp Time = ', /align_right, uname='explabel', xsize=95)
    wcurexp = widget_label(wexpbase, value='     ', /align_left, uname='expt', uvalue='expt', xsize=150)
    wzbase = widget_base(wcurobj, /align_center, /row)
    wzlabel = widget_label(wzbase, value='z = ', /align_right, uname='zlabel', xsize=95)
    wcurz = widget_label(wzbase, value='     ', /align_left, uname='curz', uvalue='curz', xsize=150)
    wzlyabase = widget_base(wcurobj, /align_center, /row)
    wzlyalabel = widget_label(wzlyabase, value='z Lya = ', /align_right, uname='zlyalabel', xsize=95)
    wcurzlya = widget_label(wzlyabase, value='     ', /align_left, uname='curzlya', uvalue='curzlya', xsize=150)
    wzismbase = widget_base(wcurobj, /align_center, /row)
    wzismlabel = widget_label(wzismbase, value='z IS = ', /align_right, uname='zismlabel', xsize=95)
    wcurzism = widget_label(wzismbase, value='     ', /align_left, uname='curzism', uvalue='curzism', xsize=150)


    ; ------ RIGHT -------
    wfile = widget_base(wright, /frame, /row, /align_left, tab_mode=1)
    wfilelist = widget_text(wfile,xsize=80,value='              ', uname='filename', uvalue='filename',tab_mode=0)
    wstatus = widget_text(wfile, xsize=108, value='Initializing ...', uname='status', uvalue='status', tab_mode=0)

    wspec = widget_base(wright, /frame, /column)
    wspecplot = widget_draw(wspec, xsize=1600, ysize=700, uname='spec', /button_events, keyboard_events=1)

    wspeccontrol = widget_base(wright, /row, /align_center, tab_mode=1)
    wycontrol = widget_base(wspeccontrol, /frame, /row)
    wylow = widget_text(wycontrol, xsize=8, /editable, uname='ylow', uvalue='ylow')
    wylabel = widget_label(wycontrol, value=' < y < ', /align_center, uname='ylabel')
    wyhigh = widget_text(wycontrol, xsize=8, /editable, uname='yhigh', uvalue='yhigh')
    wlambdacontrol = widget_base(wspeccontrol, /frame, /row)
    wblue = widget_button(wlambdacontrol, value='<-', uname='blue', uvalue='blue', /align_center)
    wlambdalow = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdalow', uvalue='lambdalow')
    wlambdalabel = widget_label(wlambdacontrol, value=' < l < ', /align_center, uname='lambdalabel', font='-urw-standard symbols l-medium-r-normal--0-0-0-0-p-0-adobe-symbol')
    wlambdahigh = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdahigh', uvalue='lambdahigh')
    wred = widget_button(wlambdacontrol, value='->', uname='red', uvalue='red', /align_center)

    widget_control, base, /realize
    self.base = base
    xmanager, 'viewspec', self.base, /no_block, cleanup='viewspec_cleanup'

    readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
    wbands = [1,2,3,4,5]
    tellstart = tellstart[wbands]
    tellend = tellend[wbands]
    ptr_free, self.linewaves, self.linewaves, self.linecolors, self.tellstart, self.tellend, self.tellthick
    self.linewaves = ptr_new([1215.670,1238.821,1260.42,1264.73,1302.168,1304.37,1334.532,1335.708,1393.76,1402.77,1526.72,1533.45,1548.202,1550.772,1608.451])
    self.linenames = ptr_new(['Lya', 'NV', 'SiII', 'SiII*', 'OI', 'SiII', 'CII', 'CII*', 'SiIV', 'SiIV', 'SiII', 'SiII*', 'CIV', 'CIV', 'FeII'])
    self.linecolors = ptr_new(['blue', 'black', 'black', 'blue', 'black', 'black', 'black', 'blue','black', 'black', 'black', 'blue', 'black', 'black', 'black'])
    self.tellstart = ptr_new(tellstart)
    self.tellend = ptr_new(tellend)
    self.tellthick = ptr_new([5, 2, 5, 2, 2])

    common random, seed
    seed = systime(1)
    self->initialize_directory
    return, 1
end

pro viewspec__define 
    state = {viewspec, $
             base:0L, $
             sciencefile:'', $
             science:ptr_new(), $
             lambdalim:[-100d, 100d], $
             ylim:[-100d, 100d], $
             divylim:[-100d, 100d], $
             skyylim:[-100d, 100d], $
             linewaves:ptr_new(), $
             linenames:ptr_new(), $
             linecolors:ptr_new(), $
             tellstart:ptr_new(), $
             tellend:ptr_new(), $
             tellthick:ptr_new(), $
             linestart:ptr_new(), $
             lineend:ptr_new(), $
             linetype:ptr_new(), $
             keystate:0, $
             lambda1:0d}
 end

pro viewspec, fitsfile
    common filename,sciencefile
    sciencefile = '/scr2/nichal/workspace3/SCIENCE/'+fitsfile

    if ~file_test(sciencefile) then message, 'File not found.'
    n = obj_new('viewspec')
end
