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
;===============================================================
function fudge_query_radius, nside, radius_in, quadratic=quadratic
;+
; radius_out = fudge_query_radius( nside, [radius_in, QUADRATIC=])
;
; with 
;   radius_out = radius_in + fudge            (default)
; or
;   radius_out = sqrt(radius_in^2 + fudge^2)  (quadratic)
;
; if absent, radius_in = 0
; where
;   fudge = factor(nside) * Pi / (4 *Nside)
; with factor = sqrt( 5/9 + (8/Pi)^2 / 5 * (1-1/(2*Nside)) )
;
; an upper bound of the actual largest radius 
; (determined by pixel located at z=2/3 and its North corner).
; 
;
; 2011-10-?, EH, v1.0
;   
;- 

radius_out = 0.0d0
if (n_params() eq 2) then radius_out = radius_in

;factor = 1.362d0
factor = sqrt( 0.55555555555555555555d0 + 1.29691115062192347448165712908d0*(1.d0-0.5d0/nside) )
fudge  = factor * !DPI / (4.0d0 * nside) ; increase effective radius
if (keyword_set(quadratic)) then begin
    radius_out = sqrt(radius_out^2 + fudge^2)
endif else begin
    radius_out += fudge
endelse
radius_out <= !DPI

return, radius_out
end
;===============================================================
