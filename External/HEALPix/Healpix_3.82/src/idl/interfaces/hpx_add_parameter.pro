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
;+
; line = hpx_add_parameter(keyword, value, ...)
;
; 2021-06-28: added a workaround for GDL, in which expand_path does not behave as IDL
;
;-
; ----------------------------------------------
function hpx_add_parameter, keyword, value, default=default, expand_path=expand, skip_if_not_set = skip, ifempty=ifempty, force_output=force_output


NoFile = " '' "

if (defined(value) || defined(default)) then begin
; value is defined or a default is provided
    if (defined(default)) then begin
        stvalue = string(set_parameter(value, default, ifempty=ifempty))
    endif else begin
        stvalue = string(value)
    endelse
    if keyword_set(expand) then begin
        if (~is_gdl() || stvalue ne NoFile) then stvalue  = expand_path(stvalue)
    endif
endif else begin
; no value, no default
    stvalue = ''
endelse

if (keyword_set(force_output)) then begin
    stvalue = add_prefix(stvalue,"!", except = [NoFile,''])
endif

line = string(keyword)+' = '+stvalue

if (keyword_set(skip) && stvalue eq '') then line = '# '+line

return, line
end
