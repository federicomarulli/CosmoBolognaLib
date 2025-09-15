; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2013  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.sourceforge.net
;
; -----------------------------------------------------------------------------
pro record_original_settings, orig_st, help=help, restore = restore
;+
;
; record_original_settings, orig_st, HELP=, RESTORE= 
; 
; will record (or restore) the !PLT settings and color table
;
;  orig_st: 
;       (output) structure containing current settings, 
;    if RESTORE is not set
;
;       (input) structure containing settings to be restored, 
;    if RESTORE is set
;
;
;
; example:
;           record_original_settings, st
;           [...]
;           record_original_settings, st, /restore
;
;
; v1.0, 2008-11-05
;-

syntax = 'Syntax: RECORD_ORIGINAL_SETTINGS, orig_st, HELP=, RESTORE= '
if (keyword_set(help)) then begin
    doc_library,'record_original_settings'
endif

if (n_params() ne 1) then begin
    print,syntax
    return
endif

if (keyword_set(restore)) then begin
    names = tag_names(orig_st)

    ; restore PLT settings, if stored in valid structure
    if (tag_names(orig_st.plt,/structure_name) eq '!PLT') then begin
        !p = orig_st.plt
    endif

    ; restore color table, if stored in valid structure
    j1 = where(names eq 'RED',   n1)
    j1 = where(names eq 'GREEN', n2)
    j1 = where(names eq 'BLUE',  n3)
    if ( (n1+n2+n3) eq 3 && (n1*n2*n3) eq 1) then begin
        tvlct, orig_st.red, orig_st.green, orig_st.blue
    endif

    ; restore original device if it was changed by this routine
    j1 = where(names eq 'DEVICE_IN',  n1)
    j1 = where(names eq 'DEVICE_OUT', n2)
    if ( (n1+n2) eq 2 && (n1*n2) eq 1) then begin
        if orig_st.device_in ne orig_st.device_out then begin
            set_plot, orig_st.device_in
            message_patch,/info,'Switching back to X windows display.', level=-1
        endif
    endif

endif else begin
    my_plt = !p
    tvlct, /get, my_red, my_green, my_blue

    my_device_in = !d.name
    my_device_out = my_device_in
    if (my_device_in eq 'X' && getenv('DISPLAY') eq '') then begin
        message_patch,/info,'X windows display requested but not properly set (DISPLAY=)',level=-1
        message_patch,/info,'Z buffer will be used instead (set_plot,"Z"), and X preview disabled.', level=-1
        my_device_out = 'Z'
        set_plot, my_device_out
    endif
    orig_st = {device_in:my_device_in, device_out:my_device_out, plt:my_plt, red:my_red, green:my_green, blue:my_blue}
endelse


return
end 
