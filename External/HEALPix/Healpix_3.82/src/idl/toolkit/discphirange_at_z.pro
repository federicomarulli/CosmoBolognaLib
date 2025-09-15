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
;===============================================================
function discphirange_at_z, vcenter, radius, z, phi0=phi0
;+
; dphi = discphirange_at_z(vcenter, radius, z, phi0=phi0)
;
;  for a disc centered on  vcenter and given radius,
;  and for location z=cos(theta)
;  returns dphi such that the disc has   abs(phi-phi0) <= dphi
;
; solving disc equation on sphere:
;  sin(theta) sin(theta0) cos(phi-phi0) + cos(theta) cos(theta0) >= cos(radius)
;
;
; INPUTS
;    vcenter: disc center (3D vector)
;
;    radius: DP scalar, disc radius in RADIANS
;
;    z: DP scalar or vector in [-1,1]
;
; OUTPUTS:
;    dphi in RADIANS in [0,Pi], same size as z
;    will be negative if z is outsize the disc
;
; OPTIONAL OUTPUTS
;    phi0: longitude of disc center in RADIANS in [-Pi,Pi]
;
;
;
; v1.0: adapted from phi_range_at_z on 2011-04
;
;-

cosang = cos(radius[0]*1.d0)

norm = sqrt(total(vcenter[0:2]^2, /double))
x0 = vcenter[0] / norm
y0 = vcenter[1] / norm
z0 = vcenter[2] / norm

phi0=0.d0
if ((x0 ne 0.d0) || (y0 ne 0.d0)) then phi0 = ATAN(y0, x0)  
;;;cosphi0 = cos(phi0)

st0 = x0*x0 + y0*y0 ; sin(theta0)^2
diff = cosang - z*z0 ; cos(rad) - cos(theta)cos(theta0)
st = 1.d0 - z*z  > 1.d-12; sin(theta)^2

dphi = replicate(-1000.d0, n_elements(z)) ; out of disc by default

if (st0 eq 0.d0) then begin ; polar cap
    full = where(diff le 0.d0, nfull) ; fully in cap
    if (nfull gt 0) then dphi[full] = !DPI

endif else begin
    cosdphi = diff/sqrt(st0*st)

    full    = where(cosdphi      lt (-1.0d0), nfull)   ; fully in disc
    if (nfull    gt 0) then dphi[full]    = !DPI

    partial = where(abs(cosdphi) le 1.0d0,    npartial) ; partially in disc
    if (npartial gt 0) then dphi[partial] = acos(cosdphi[partial]) ; in [0,Pi]

endelse

return, dphi
end
