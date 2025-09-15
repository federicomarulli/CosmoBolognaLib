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
; NAME:
;   DECODE_ORDERING
;
;
; PURPOSE:
;   returns valid string (or number) identifying the RING or NESTED ordering scheme
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;   output = decode_ordering(in_string [,Full_Name=, Short_Name=, Number=,
;   Error=, Silent=, ToFull=, ToNumber])
;
; 
; INPUTS:
;   in_string: user provided string (case un-sensitive)
;      examples of valid input: 'ring' 'RiNg', 'ringgg', ...
;                              'nested' 'NEst', 'nesttt', ...
;
;      if non valid string is input, an error message is issued (unless /Silent
;      is set)
;
;      if a non-string or a undefined variable is provided an error message will
;      be issued too
;
; OPTIONAL INPUTS:
;    none
;
; KEYWORD PARAMETERS:
;    Silent: is set, the error flag is set to 1 in case of error, but no
;            error message is issued
;    ToFull: a full length string ('RING' or 'NESTED') is issued instead of
;    'RING' and 'NEST'
;    ToNumber: a number (1 or 2) is issued instead of 'RING' or 'NEST' respectively
;
; OUTPUTS:
;    output: either 'RING' or 'NEST' (or 'RING' or 'NESTED' ; or 1 or 2)
;
; OPTIONAL OUTPUTS:
;    Full_Name: contains full length string on output(RING/NESTED/Error)
;    Short_Name: (RING/NEST/ERRO)
;    Number:     (1/2/-1)
;    Error: set to 1 if input string unvalid; 0 if valid input
;
; COMMON BLOCKS:
;    none
;
; SIDE EFFECTS:
;    none
;
; RESTRICTIONS:
;    none
;
; PROCEDURE:
;    string editing and parsing
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;   2006-11-20: EH, IAP, 1st version
;-

function decode_ordering, ordin, full_name=full_name, short_name=short_name, number=number, error=error, silent=silent, tofull=tofull, tonumber=tonumber

if undefined(ordin) then begin
    message,'Undefined ORDERING string'
endif

if datatype(ordin) ne 'STR' then begin
    message,'Expected a string for ORDERING, got:'+datatype(ordin)
endif

short_name = strtrim( ordin, 2)  ; remove blanks
short_name = strupcase( strmid( short_name, 0, 4)) ; first 4 letters, upper case

error = 0
case short_name of
    'RING': begin
        full_name = 'RING'
        number = 1
    end
    'NEST': begin
        full_name = 'NESTED'
        number = 2
    end
    else: begin
        short_name = 'ERROR'
        full_name = 'Error'
        number = -1
        error = 1
    end
endcase

if (error eq 1 and keyword_set(silent) eq 0) then begin
    message,'Invalid ORDERING string: '+ordin
endif

; output short_name by default
out = short_name

; or full name
if keyword_set(tofull) then out = full_name

; or number
if keyword_set(tonumber) then out = number

return, out
end
