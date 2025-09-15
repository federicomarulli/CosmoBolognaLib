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
function add_prefix, file_in, prefix, except=except, suffix=suffix
;+
; file_out = add_prefix(file_in, prefix, Except=, Suffix=)
;
;  file_in : string scalar
;  prefix : string scalar
;
;  except : string vector
;
;  suffix: flag
;
;  prepend prefix to str_in (prefix+str_in) except if 
;    - the prefix is already in place OR
;    - file_in does not belong to Except
;
; if Suffix is set, the result is (str_in+prefix) except if
;     - the suffix is already in place OR
;     - file_in is the same as Except
;
; June 2007
;-
str_in = strtrim(file_in,2)

file_out = str_in

; is input string among the 'except' list ?
nfound = 0
if defined(except) then junk = where(strtrim(except,2) eq strtrim(str_in,2), nfound)

; if not, then proceed
if (nfound eq 0) then begin
    if (keyword_set(suffix)) then begin
        suff = strtrim(prefix,2)
        if (strpos(str_in, suff,/reverse_search) ne strlen(str_in)-strlen(suff)) then file_out = str_in+suff
    endif else begin
        pref = strtrim(prefix,2)
        if (strpos(str_in, pref) ne 0) then file_out = pref+str_in
    endelse
endif

return, file_out
end

