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
pro add_bad_data_fits, info_header, bad_data=bad_data
;+
; add_bad_data_fits, info_header [, Bad_data= ]
;
;  deals with BAD_DATA FITS keyword
;
; if user provided bad_data is NaN, remove 'BAD_DATA' keyword from header.
; otherwise use user provided value, when applicable
;   if not, use !healpix.bad_value
;
;
; created 2010-05-11, EH
;
;-


kwd_flag = 0
kwd = 'BAD_DATA'

kwd_value = sxpar(info_header, kwd, count=kwd_flag)

if keyword_set(bad_data) then begin
    if finite(/NaN, bad_data) then begin
        sxdelpar, info_header, kwd
        return
    endif
endif else begin
    init_healpix
    bad_data = !healpix.bad_value
endelse

sxaddpar, info_header, kwd, bad_data, ' Sentinel value given to bad pixels'


return
end

