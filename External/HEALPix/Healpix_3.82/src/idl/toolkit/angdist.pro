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
pro angdist, v1, v2, dist

;+
; angdist, v1, v2, dist
; computes the angular distance dist (in rad) between 2 vectors v1 and v2
; in general dist = acos ( v1 . v2 )
; except if the 2 vectors are almost aligned.
;
;
; HISTORY : Oct 2003,
;    EH, corrected bug on scalar product boundary
;    2010-03-01: works for non-normalized input vectors
;    2011-04-11: bug corrected in special cases test
;    2011-04-12: replaced with angulardistance function
;-

dist = -100000.d0
if (n_params() ge 2) then dist = angulardistance(v1, v2)
message,/info,' The routine ANGDIST is now obsolete and has been replaced by the ANGULARDISTANCE function,'
message,/info,'which accepts list of vectors.'
message,/info,'Replace'
message,/info,'           angdist, vec1, vec2, distance '
message,/info,'with'
message,/info,'           distance = angulardistance(vec1, vec2) '
; if (n_params() ne 3) then begin
;     print,'Syntax : ANGDIST, vec1, vec2, angular_dist'
;     return
; endif


; ; normalize both vectors
; r1 = v1 / sqrt(total(v1*v1))
; r2 = v2 / sqrt(total(v2*v2))

; ; scalar product
; ;sprod = total(v1*v2)
; sprod = total(r1*r2) ; bug corrected 2010-03-01

; if (sprod gt 0.999d0) then begin ; replaced > with gt
; ; almost colinear vectors
;     vdiff = r1 - r2
;     diff = sqrt(total(vdiff*vdiff)) ; norm of difference
;     dist = 2.0d0 * asin(diff * 0.5d0)

; endif else if (sprod lt (-0.999d0)) then begin ; replaced < with lt
; ; almost anti-colinear vectors
;     vdiff = r1 + r2
;     diff = sqrt(total(vdiff*vdiff)) ; norm of sum
;     dist = !DPI - 2.0d0 * asin(diff * 0.5d0)

; endif else begin
; ; other cases
;     dist = acos( sprod )
; endelse

return
end


