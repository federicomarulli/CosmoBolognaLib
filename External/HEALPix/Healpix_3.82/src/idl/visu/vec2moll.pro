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
pro vec2moll, vec, u, v, flip=flip
;+
; VEC2MOLL
;  returns u,v position on Mollweide map corresponding to 3D vector 'vec'
; 
; VEC2MOLL, vec, u, v
;
; 2010-04-01: use double precision variables
;             make sure output of single precision INTERPOL is within valid range
;-

if keyword_set(flip) then flipconv = +1 else flipconv = -1 ; longitude increase leftward by default (astro convention)
halfpi = !dpi / 2.d0

common vec2mollcom, ntab, vtab, ztab
if undefined(ntab) then begin
    ntab = 180L * 10
    vtab = 2.d0*DINDGEN(ntab)/(ntab-1.d0) - 1.d0 ; in [-1,1]
    ztab = (ASIN(vtab) + vtab*SQRT(1.d0-vtab^2))/halfpi ; in [-1,1]
endif


np1 = n_elements(vec)/3
vec = reform(vec,np1,3,/Over)

phi = ATAN(vec(*,1), vec(*,0))/ halfpi ; in [-2, 2]

v = INTERPOL(vtab,ztab,vec(*,2)/SQRT(TOTAL(vec^2,2)))
v <= 1.d0  ; interpol is single precision -> |v| may exceed 1.
v >= (-1.d0)
u = flipconv* phi * SQRT(1.d0-v^2) ; minus sign : astro convention

return
end
