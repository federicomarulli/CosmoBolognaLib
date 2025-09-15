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
pro add_units_fits,info_header, units=units,colnum=colnum,error=error
;+
; add_units_fits
;
; add_units_fits,info_header, units=,colnum=error=
;
;  adds the UNITS keyword defined by 'units' in
;  the fits header 'info_header'
;
;  EH 2000-01
;  2005-02-08 : remove creation of header primer if no info_header provided
;-

if (datatype(units) ne 'STR' and datatype(units) ne 'UND') then begin
    print,'invalid UNITS'
    error=1
    return
endif
error=0

; info_header present and not empty
header_flag = 0
if (n_params() eq 1) then begin
    if defined(info_header) then header_flag=1
endif

if undefined(colnum) then colnum = [1]
scolnum = strtrim(string(colnum,form='(i4)'),2)
ncolnum = n_elements(colnum)

; units given in the header
units_flag = 0
if (header_flag) then units_fits = strtrim(sxpar(info_header,'TUNIT'+scolnum[0],count=units_flag),2)
if units_flag eq 0 then begin
    units_fits = ' '
    if undefined(units) then units='unknown'
endif else begin
    if undefined(units) then return
endelse

; no header, open one
;;; if (header_flag eq 0) then info_header = strarr(1) 

; change TUNIT*
if (units_fits ne units ) then begin
    if (strtrim(units_fits) ne '') then print,' value of UNITS changed from '+units_fits+' to '+units
    for i=0,ncolnum-1 do begin
        sxaddpar,info_header,'TUNIT'+scolnum[i],units
    endfor
endif


return
end

;--------------------------------
