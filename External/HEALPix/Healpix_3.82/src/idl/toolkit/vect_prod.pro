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
function vect_prod, v1, v2, help=help
;+
;=======================================================================
;     v3=vect_prod( v1, v2  [,/help] )
;     returns in v3 the vectorial product of the two 3-dimensional 
;     vectors v1 and v2
;=======================================================================
;-

routine = 'vect_prod'
syntax = 'v3 = '+routine+'( v1, v2 [, /HELP] )'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if n_params() ne 2 then begin
    print,syntax
    return,-1
endif

v3 = dblarr(3)
v3[0] = v1[1] * v2[2] - v1[2] * v2[1]
v3[1] = v1[2] * v2[0] - v1[0] * v2[2]
v3[2] = v1[0] * v2[1] - v1[1] * v2[0]
    
return, v3
end
