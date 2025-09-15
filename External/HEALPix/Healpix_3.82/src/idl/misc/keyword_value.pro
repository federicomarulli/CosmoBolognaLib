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
function keyword_value, keyword, default=default, HELP=help, minvalue=minvalue, maxvalue=maxvalue
;+
;
;  value = keyword_value(Keyword, DEFAULT=, HELP=, MINVALUE=, MAXVALUE=)
;
;  returns Keyword value,
;    or DEFAULT if Keyword is not defined, 
;    or 0 otherwise
;  The result will be bounded by MINVALUE and/or MAXVALUE if one or both are
;  provided
;
;  the consistency of MINVALUE and MAXVALUE is *not* checked
;
;  2013-04-11
;
;
;-

routine = 'keyword_value'
syntax  = 'value = '+routine+'(Keyword, DEFAULT=, HELP=, MINVALUE=, MAXVALUE=)'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if n_params() eq 0 then begin
    print,syntax
    return,-1
endif

mydef = keyword_set(default)? default : 0

value = defined(keyword)    ? keyword : mydef

if defined(minvalue) then value >= minvalue
if defined(maxvalue) then value <= maxvalue

return,value
end
