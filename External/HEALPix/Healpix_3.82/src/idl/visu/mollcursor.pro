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

;=============================================================================

PRO mollcursor, cursor_type=cursor_type, file_out = file_out_user
;+
;=-----------------------------------------------------------------------------
; NAME:   
;     MOLLCURSOR
; PURPOSE:
;     get sky position (long, lat), number and value of a pixel
;     selected on a Healpix or QuadCube map in Mollview projection 
; CALLING SEQUENCE:
;       MOLLCURSOR, [CURSOR_TYPE=, FILE_OUT=]
; INPUTS:
;      None
; OPTIONAL INPUT:
; OPTIONAL KEYWORD
;      cursor_type
;      file_out
; OUTPUTS:
; OPTIONAL OUTPUTS : 
;      Long, Lat, value, pixel number
; SIDE EFFECT:
; RESTRICIONS:
;      runs on the X-window created by the latest execution of Mollview
; PROCEDURE CALLS:
;      Moll2pix
; HISTORY:
;       June 98 : copied by EH from 
;        coor_cursor (F.X.Desert, JUL-1992) and xcursor (R. Sterner, Nov-1993)
;       March 1999 : updated for improved Mollview
;       Sept 2000 : updated for polarisation
;=-----------------------------------------------------------------------------
;-

@viewcom ; define common
if (!d.window lt 0 or undefined(data_plot)) then begin
    print,'****************************************************'
    print,'  This should be run *after* GnomView/MollView/Orthview/Cartview !'
    print,'****************************************************'
    return
endif
do_polamplitude = (polar_type eq 1)
do_poldirection = (polar_type eq 2)
do_polvector    = (polar_type eq 3)

do_file = 0
if defined(file_out_user) then begin
    file_out_def = 'cursor_catalog.txt'
    file_out = (datatype(file_out_user) eq 'STR' and strtrim(file_out_user,2) ne '') ? file_out_user : file_out_def
    do_file = 1
endif

; Initialisation
N_POINTS= 1L

coordinate = decode_coord(coord_out,nlong=long_name,nlat=lat_name)
ilo = strlen(strtrim(long_name,2))
ila = strlen(strtrim(lat_name ,2))
if (ilo eq 2) then long_name = '  '+strtrim(long_name,2)
if (ilo eq 3) then long_name = ' '+strtrim(long_name,2)
if (ila eq 3) then  lat_name = ' '+strtrim( lat_name,2)

svalue = '         Pixel Value'
if (do_polamplitude) then svalue = '      Pol. amplitude'
if (do_poldirection) then begin
    svalue = '      Pol. angle'
    sunits = 'Deg'
endif

blank = '                    '
len1 = STRLEN(svalue)
if (sunits NE ' ') then begin
    sunits2 = STRTRIM(sunits,2)
    len2 = STRLEN(sunits2)
;     svalue = STRMID(blank,0,(len1-len2-8)>1)+'Value ('+sunits2+')'
;     svalue = strmid(svalue,len1+1-len2-20,len2-2+20)+'('+sunits2+')'
    svalue = strmid(svalue,len1+2-len2-20,len2-1+20)+'('+sunits2+')'
endif

pix_num0 = ' pix #'
if (pix_type eq 'R') then pix_num = '  Ring'+pix_num0
if (pix_type eq 'N') then pix_num = '  Nest'+pix_num0
if (pix_type eq 'Q') then pix_num = '  Q.C.'+pix_num0


string_instruct = '  Left: display position,   Middle: record position,   Right: Quit  '
if (do_file eq 1) then string_recfile  = '  Record in file : '+file_out
string_head     = '   # '+long_name+' (Deg)  '+lat_name+' (Deg)'+svalue+pix_num
string_coord    = coordinate + ' coordinates'

;--------------------------------------------------------------------
do_moll = 0
do_gnom = 0
do_cart = 0
do_orth = 0
if (strmid(projection,0,4) eq 'MOLL') then begin
    title = 'MollView: cursor information'
    do_moll = 1
endif
if (strmid(projection,0,4) eq 'GNOM') then begin
    title = 'GnomView: cursor information'
    do_gnom = 1
endif
if (strmid(projection,0,4) eq 'CART') then begin
    title = 'CartView: cursor information'
    do_cart = 1
endif
if (strmid(projection,0,4) eq 'ORTH') then begin
    title = 'OrthView: cursor information'
    do_orth = 1
endif

; set map window up-front
wset, w_num
wshow, w_num, 1, iconic=0 ; put the current window above all other

; create widget window
top = widget_base(/column,title=title)
id = widget_label(top,val=string_instruct)
if (do_file) then id = widget_label(top,val=string_recfile)
id = widget_label(top,val=string_coord)
id = widget_label(top,val=string_head)
id = widget_label(top,val=' ',/dynamic_resize)
widget_control, top, /real
widget_control, id, set_val=' ', show=1

PRINT
PRINT,'Instructions :'
PRINT,string_instruct
if (do_file) then print,string_recfile

if not keyword_set(cursor_type) then ct= 34 else ct= cursor_type
device, cursor_standard= ct
new_point:
;--------

cursor, xpos, ypos, /wait, /data
; button = !err ; until idl 3.6
button = !mouse.button
if (button eq 4) then goto, closing
if (button gt 3) then goto, new_point

if (do_moll) then moll2pix, xpos, ypos, id_pix, lon_deg, lat_deg, value
if (do_gnom) then gnom2pix, xpos, ypos, id_pix, lon_deg, lat_deg, value
if (do_cart) then cart2pix, xpos, ypos, id_pix, lon_deg, lat_deg, value
if (do_orth) then orth2pix, xpos, ypos, id_pix, lon_deg, lat_deg, value
if (id_pix ge 0) then begin
    if (do_poldirection) then begin
        value =  !RaDeg*value MOD 180. ; in [0,180]
        value = value - 180.*(value GT 90.)  ; in [-90,90]
        format = '(i4,f11.2,f11.2,g20.3,i12)'
    endif else begin
        format = '(i4,f11.2,f11.2,e20.5,i12)'
    endelse
    string_out = string(n_points, lon_deg, lat_deg, value, id_pix,format=format)
    if (button eq 2) then begin
        if n_points eq 1 then begin
            list_pix = id_pix 
            print
            print,'selected pixels ('+string_coord+')'
            print,string_head
        endif else begin
            list_pix = [list_pix, id_pix]
        endelse
        print,string_out
        if (do_file eq 1) then begin
            if (n_points eq 1) then begin
                openw,lunit,file_out,/get_lun
                printf,lunit,'# '
                printf,lunit,'# '+string_coord
                printf,lunit,'# Resolution Parameter : '+string(pix_param,form='(i5)')
                printf,lunit,'# '
                printf,lunit,'#'+string_head
            endif
            printf,lunit,string_out
        endif
        n_points= n_points + 1
    endif
    widget_control, id, set_val=string_out, show=1
endif else begin
    string_out = string(replicate('**',5),format='(a4,a11,a11,a20,a12)')
    widget_control, id, set_val=string_out, show=1
endelse

; Suspend going on for a while
WAIT, 0.5

GOTO, NEW_POINT

CLOSING:
;------
if (n_points gt 1 and do_file eq 1) then free_lun, lunit
widget_control, top, /dest
DEVICE, /CURSOR_ORIGINAL
;if (n_points gt 1) then print,list_pix
print

RETURN
END



