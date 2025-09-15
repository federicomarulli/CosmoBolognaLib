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
function ring2z, nside, iz
;+
; ring = ring2z(nside, iz)
;
; for a given Nside, returns the z=cos(theta) of a ring index (in [1,4*nside-1])
;
;  v1.0: ?
;  v1.1: 2011-04-07: accepts vector iz
;
;-

routine = 'ring2z'
syntax = 'ring = '+routine+'(nside, iz)'

if n_params() ne 2 then begin
    print,syntax
    return,-1
endif


range = hpx_minmax(iz)
if (range[0] lt 0 || range[1] gt 4*nside) then begin
    print,syntax
    message,' iz must be in [0,4*Nside]'
endif
if (n_elements(nside) gt 1) then begin
    print,syntax
    message,' Nside must be a scalar integer'
endif


dth1 = 1.0d0 / (3.0d0*double(nside[0])^2)
dth2 = 2.0d0 / (3.0d0*double(nside[0]))

if (n_elements(iz) eq 1) then begin ; scalar iz (original code)

    if (iz le nside-1) then begin ; north polar cap
        zring = 1.0d0  - double(iz)^2 * dth1
    endif else if (iz le 3L*nside) then begin ; tropical band + equat.
        zring = double(2L*nside-iz) * dth2
    endif else begin
        zring = - 1.0d0 + double(4L*nside-iz)^2 * dth1
    endelse

endif else begin ; vector iz

    npix = nside2npix(nside, err=err)
    if (err ne 0) then begin
        print, syntax
        message,'expected valid Nside'
    endif
    zring = dblarr(n_elements(iz))

    count = histogram(iz, min = -nside, max=5*nside, bin=2*nside, rev=rev)

    i = 0
    if (count[i] gt 0) then begin ; north polar cap
        k = rev[rev[i]:rev[i+1]-1]
        zring[k] = 1.0d0  - double(iz[k])^2 * dth1
    endif

    i = 1
    if (count[i] gt 0) then begin ; tropical band + equat.
        k = rev[rev[i]:rev[i+1]-1]
        zring[k] = double(2L*nside-iz[k]) * dth2
    endif

    i = 2
    if (count[i] gt 0) then begin ; south polar cap
        k = rev[rev[i]:rev[i+1]-1]
        zring[k] = - 1.0d0 + double(4L*nside-iz[k])^2 * dth1
    endif

endelse

return, zring
end
