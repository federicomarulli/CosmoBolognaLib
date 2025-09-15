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
function my_asinh, x, type
;
; if type=1, returns asinh(x) [default]
;   which is ~x    when x<<1 and ~ln(2x)   when x >> 1
; if type=2, returns asinh(x/2)/ln(10) 
;  which is ~x*0.21 for x<<1 and ~log10(x) when x >> 1
;
mytype = defined(type) ? type : 1
case mytype of
    1: y = asinh(x)
    2: begin
        two = (size(x[0],/tname) eq 'DOUBLE') ? 2d : 2.0
        y = asinh(x/two)/alog(two*5)
    end
    else: message,'type must be 1 or 2'
endcase
return, y
end
;----------------------------------
function my_sub_interpolate, data, x
;
; interpolate function that keeps working for scalar data
;
if n_elements(data) eq 1 then data = [data]
y = interpolate(float(data), x)
return, y
end

;----------------------------------

function color_map, data, mindata, maxdata, Obs, $
                    color_bar = color_bar, mode = mode, minset = min_set, maxset = max_set, silent = silent
;+
; color_map
;
; color = color_map( data, mindata, maxdata, Obs=obs, color_bar =
; color_bar, mode = mode)
;
;
; INPUT
;    data : input data
;    mindata : min (can be changed by the routine)
;    maxdata : maxdata (can be changed by the routine)
;
;  Obs : list of valid pixels
;  if Obs is not defined, all pixels are assumed valid
;
; mode    0: nothing
;         1: histogram equalized
;         2: log
;         4: asinh
; asin can not be used together with log nor hist_eq
;
;
;
; Sep 2007: added /silent
; July 2008: added Asinh mode
; Oct 2008: make sure that MIN and MAX are properly taken into account in Asinh
; mode
; 2009-04-30: make sure that MINSET and MAXSET are taken into account in Log mode
;-

if undefined (mode) then mode = 0
do_hist  = (mode and  1)   ; lowest bit -> [0,1]
do_log   = (mode and  2)/2 ; 2nd bit    -> [0,1]
tasinh   = (mode and 12)/4 ; 3rd and 4th bit -> 0, 1 or 2

if ((tasinh gt 0) && (do_log || do_hist)) then begin
    message,'Asinh mode can NOT be used together with Log or Hist mode'
endif

N_Color = !D.n_colors < 256

sz = size(data)
npix = n_elements(data)
Color = MAKE_ARRAY(/BYTE, npix, Value = 2B)
color_bar = [2B]
if (sz(0) eq 2) then Color = reform(color,/over,sz(1),sz(2))
if (sz(0) eq 3) then Color = reform(color,/over,sz(1),sz(2),sz(3))
if defined(Obs) then begin
    if Obs[0] eq -1 then return, color
    N_obs = n_elements(Obs)
    N_no_Obs = npix - N_obs
endif else begin
    ;Obs = lindgen(npix)
    N_obs = npix
    N_no_obs = 0
endelse

; -------------------------------------------------------------
; sets MIN and MAX
; -------------------------------------------------------------
if (~keyword_set(silent)) then print,'plotted area original MIN and MAX: ',mindata, maxdata
IF DEFINED(min_set) THEN BEGIN
    if (min_set gt mindata) then begin
        IF (N_No_Obs eq 0) THEN data = data > min_set ELSE data(Obs)=data(Obs) > min_set
    endif
    mindata = min_set
    if (~keyword_set(silent)) then print,'new MIN : ',mindata
ENDIF

IF DEFINED(max_set) THEN BEGIN
    if (max_set lt maxdata) then begin
        IF (N_No_Obs eq 0) THEN data = data < max_set ELSE data(Obs)=data(Obs) < max_set
    endif
    maxdata = max_set
    if (~keyword_set(silent)) then print,'new MAX : ',maxdata
ENDIF

IF (do_log) THEN BEGIN
;   log
    if (N_No_Obs eq 0) then begin
;         data      = ALOG10(data + (0.0001 -mindata))
        data = Alog10(data > 1.e-6*abs(maxdata))
        mindata = MIN(data,MAX=maxdata)
    endif else begin 
;         data[Obs] = ALOG10(data[Obs] + (0.0001 -mindata))
        data[Obs] = ALOG10(data[Obs] > 1.e-6*abs(maxdata))
        mindata = MIN(data[Obs],MAX=maxdata)
    endelse
    ; use user-defined range if provided
    if defined(min_set) then begin
        if (min_set gt 0) then mindata = alog10(min_set)
    endif
    if defined(max_set) then begin
        if (max_set gt 0) then maxdata = alog10(max_set)
    endif
ENDIF

; IF (do_asinh) THEN BEGIN
; ;   Asinh
;     if (N_No_Obs eq 0) then begin
;         data = asinh(data)
;         mindata = MIN(data,MAX=maxdata)
;     endif else begin 
;         data[Obs] = asinh(data[Obs])
;         mindata = MIN(data[Obs],MAX=maxdata)
;     endelse
; ENDIF


; turn data into colors (in the range [3,N_color-1])
col_scl = FINDGEN(N_Color-3)/(N_Color-4)

if (do_hist) then begin
;   histogram equalised scaling
    Tmax = maxdata & Tmin = mindata
    Bs = (Tmax-Tmin)/5000. ;MIN( [(Tmax-Tmin)/5000.,2.] )
    if (Bs eq 0.) then Bs=1.d-30 ; robust to narrow (or zero) range maps
    if (N_No_Obs eq 0) then begin
        Phist = HISTOGRAM( data, Min=Tmin, Max=Tmax, Bin = Bs )
        Phist = total(Phist,/cumul)
        Phist[0] = 0.
        Junk = my_sub_interpolate( Phist, (data-Tmin)/Bs )
        Color = 3B + BYTSCL( Junk, Top=N_Color-4 )
    endif else begin 
        Phist = HISTOGRAM( data[Obs], Min=Tmin, Max=Tmax, Bin = Bs )
        Phist = total(Phist,/cumul)
        Phist[0] = 0.
        Junk = my_sub_interpolate( Phist, (data[Obs]-Tmin)/Bs )
        Color[Obs] = 3B + BYTSCL( Junk, Top=N_Color-4 )
    endelse

    junk2= my_sub_interpolate( Phist, col_scl*(Tmax-Tmin)/bs )
    color_bar = (3B + BYTSCL( junk2, TOP = N_Color-4 ))

endif else begin
;   linear scaling or Asinh mapping
    if (tasinh gt 0) then begin
        Tmax = maxdata*1.0 & Tmin = mindata*1.0
        if (N_No_Obs eq 0) then begin
            Color = 3B + BYTSCL( my_asinh(data, tasinh), MIN=my_asinh(Tmin, tasinh), MAX=my_asinh(Tmax, tasinh), Top=N_Color-4 )
        endif else begin 
            Color[Obs] = 3B + BYTSCL( my_asinh(data[Obs], tasinh), MIN=my_asinh(Tmin, tasinh), MAX=my_asinh(Tmax, tasinh), Top=N_Color-4 )
        endelse
        
        bs = 5 * N_color
        junk2= my_asinh(  dindgen(bs)/(bs-1)*(Tmax-Tmin) + Tmin, tasinh  )
        color_bar = (3B + BYTSCL( junk2, TOP = N_Color-4 ))

    endif else begin
        if (maxdata ne mindata && $
            ABS((maxdata+mindata)/FLOAT(maxdata-mindata)) lt 5.e-2) then begin
;       if Min and Max are symmetric
;       put data=0 at the center of the color scale
            Tmax = MAX(ABS([mindata,maxdata]))
            Tmin = -Tmax
        endif else begin 
            Tmax = maxdata & Tmin=mindata
        endelse
        if (N_No_Obs eq 0) then begin 
            color = 3B + BYTSCL(data, MIN=Tmin, MAX=Tmax, Top=N_Color-4 )
        endif else begin
            color[Obs] = 3B + BYTSCL(data[Obs], MIN=Tmin, MAX=Tmax, Top=N_Color-4 )
        endelse
        color_bar = (3B + BYTSCL( col_scl, TOP = N_Color-4 ))
    endelse
endelse

mindata = Tmin
maxdata = Tmax

return, color
end

