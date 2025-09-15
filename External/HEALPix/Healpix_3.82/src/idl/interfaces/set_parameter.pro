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
function set_parameter, param_in, default, ifempty=ifempty, ifzero=ifzero
;+
;  value = set_parameter(Input_Value, Default, IFEMPTY=, IFZERO=)
;
;  if input_value in undefined 
;   (or is an empty string when IFEMPTY is set, 
;    or 0 when IFZERO is set) 
;  then value=default,
;  in all other cases, value=input_value
;-

if n_params() ne 2 then begin
    print,'value = set_parameter(Input_Value, Default, IFEMPTY=, IFZERO=)'
    message,'Abort'
endif

param_out = defined(param_in) ? param_in : default

if keyword_set(ifempty) then begin
    if (size(param_out,/tname) eq 'STRING' && strtrim(param_out,2) eq '') then param_out = default
endif

if keyword_set(ifzero) then begin
    if (n_elements(param_out) eq 1 && size(param_out,/tname) ne 'STRING' && param_out eq 0) then param_out = default
endif


return, param_out
end


