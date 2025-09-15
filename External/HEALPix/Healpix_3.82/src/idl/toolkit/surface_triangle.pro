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
function surface_triangle, vec1, vec2, vec3
;+
;=======================================================================
; surface = surface_triangle( vec1, vec2, vec3)
;
; returns the surface in steradians 
;  of the spherical triangle with vertices vec1, vec2, vec3
;
; algorithm : finds triangle sides and uses l'Huilier formula to compute
; "spherical excess" = surface area of triangle on a sphere of radius one
; see, eg Bronshtein, Semendyayev Eq 2.86
;
;
;  Original code:?
;  2011-04-12: replace angdist with angulardistance
;=======================================================================
;-

;   half perimeter  
;   hp = 0.5_dp * (side1 + side2 + side3)

;   ; l'Huilier formula
;   x0 = tan( hp          * 0.5_dp)
;   x1 = tan((hp - side1) * 0.5_dp)
;   x2 = tan((hp - side2) * 0.5_dp)
;   x3 = tan((hp - side3) * 0.5_dp)
;  surface = 4 * atan( sqrt(x0 x1 x2 x3) )

; find triangle sides, divided by 4
side1 = (angulardistance(vec2, vec3))[0] * 0.25d0
side2 = (angulardistance(vec3, vec1))[0] * 0.25d0
side3 = (angulardistance(vec1, vec2))[0] * 0.25d0

; l'Huilier formula
x0 = tan( side1 + side2 + side3 )
x1 = tan(-side1 + side2 + side3 )
x2 = tan( side1 - side2 + side3 )
x3 = tan( side1 + side2 - side3 )
surface = 4.d0 * atan( sqrt(x0 * x1 * x2 * x3) )

return, surface
end
