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
function ring_num, nside, z, shift=ishift, withpoles=withpoles
;+
; ring=ring_num(nside, z [, shift=, withpoles=])
;
; gives the ring number corresponding to z for the resolution nside
;
; usually returns the ring closest to the z provided
; if shift < 0, returns the ring immediatly north (of smaller index) of z
; if shift > 0, returns the ring immediatly south (of larger index) of z
;
; 2008-03-28: accepts scalar and vector z
;             added shift
;
; 2009-04-30: returns long integers (32 bits) instead of 16-bit integers
; 2011-03-10: added withpoles= to return rings in [0,4Nside] instead of [1, 4Nside-1]
; 2018-07-10: activates withpoles
;-
twothird = 2.d0 /3.d0
long = 1

shift = 0.d0
if (keyword_set(ishift)) then begin
    if (ishift lt 0) then shift = -0.5d0 ; Northward shift
    if (ishift gt 0) then shift =  0.5d0 ; Southward shift
endif
lower = 1
if keyword_set(withpoles) then lower = 0

;     ----- equatorial regime ---------
iring = NINT( nside*(2.d0-1.500d0*z) + shift, long=long)

;     ----- north cap ------
kn = where(z gt twothird, nkn)
if (nkn gt 0) then begin
    my_iring = NINT( nside* SQRT(3.d0*(1.d0-z[kn]))  + shift, long=long )
    iring[kn] = my_iring > lower
endif

;     ----- south cap -----
ks = where (z lt -twothird, nks)
if (nks gt 0) then begin
    ; beware that we do a -shift in the south cap
    my_iring = NINT( nside* SQRT(3.d0*(1.d0+z[ks]))  - shift, long=long)
    my_iring = my_iring > lower
    iring[ks] = 4L*nside - my_iring
endif

if (n_elements(z) eq 1) then begin
    return, iring[0]
endif else begin
    return, iring
endelse

end
