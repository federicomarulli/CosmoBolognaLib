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
function today_fits, utc=utc, time=time
;+
; today_fits
;  return the current date in the FITS Y2K compliant format
;  ccyy-mm-dd
; OR
;  ccyy-mm-ddThh:mm:ss
;
; result = today_fits(UTC=, TIME=)
;
;  if UTC is set, will return UTC time instead of current time zone time
;
;  if TIME is set, will return date and time, instead of time alone
;
; mainly copied from today() from IDL 5.2
;
; updated to work with IDL 5.5 (str_sep is discontinued)
;
; 2017-10-23: added UTC and TIME keywords
;
;-

ascii_time  = systime(0, utc=utc)         ;Current time
a_time      = strcompress(ascii_time)  ; compress all those spaces to 1
idl_version = float(!version.release)
if (idl_version le 5.4) then begin
    s = str_sep(a_time, ' ')    ;Separate fields
    t = str_sep(s[3], ':')      ;Time fields separated by colon
endif else begin
    s = strsplit(a_time, ' ',/extract)    ;Separate fields
    t = strsplit(s[3], ':',/extract)      ;Time fields separated by colon
endelse
m = where(strupcase(s[1]) eq $  ; = month  - 1
 ['JAN','FEB','MAR','APR', 'MAY', 'JUN', 'JUL', 'AUG','SEP','OCT','NOV','DEC'])

fdate = STRING([s[4],m[0]+1,s[2]],form='(i4.4,''-'',i2.2,''-'',i2.2)')

if keyword_set(time) then begin
    fdate += 'T'+s[3]
endif


return,fdate
end

