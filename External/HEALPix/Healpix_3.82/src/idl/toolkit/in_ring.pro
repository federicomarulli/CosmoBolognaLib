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
function in_ring, nside, iz, phi0, dphi, nir, nested=nested, conservative=conservative
;+
;
; result = in_ring(nside, iz, phi0, dphi, nir)
; gives the list of pixels contained in [phi0-dphi,phi0+dphi]
; on the ring iz for the resolution nside
; nir is the number of pixels found
; if no pixel is found, on exit nir =0 and result = -1
;
;  bug corrected on 1999-12-07 (turned ip_low and ip_hi to long, add
;  nir on output)
;
;  correction for bug at low phi, EH, 2004-05-28, 2004-06-01
;  get ready for nside > 8192, exit if dphi<0, 2008-03-29
;-

if (dphi lt 0) then begin
    nir = 0
    return, -1
endif

twopi = 2.d0*!DPI
take_all = 0 ; false
to_top = 0   ; false
;conservative = 0
do_nest = (keyword_set(nested))

npix = nside2npix(nside)
if (nside gt 8192) then begin
    one = 1LL
endif else begin
    one = 1L
endelse
ncap  = 2*nside*(nside-one) ; number of pixels in the north polar cap
listir = -1
nir = 0

phi_low = (phi0 - dphi) MOD twopi
if (phi_low lt 0) then phi_low = phi_low + twopi
phi_hi  = (phi0 + dphi) MOD twopi
if (phi_hi lt 0)  then phi_hi  = phi_hi  + twopi
take_all = (ABS(dphi-!DPI) lt 1.d-6)

;     ------------ identifies ring number --------------
if ((iz ge nside) && (iz le 3*nside)) then begin ; equatorial region
    ir = iz - nside + 1L  ; in {1, 2*nside + 1}
    ipix1 = ncap + 4L*nside*(ir-one) ;  lowest pixel number in the ring
    ipix2 = ipix1 + 4L*nside - one   ; highest pixel number in the ring
    kshift = ir MOD 2
    nr = nside*4L
endif else begin
    if (iz lt nside) then begin       ;    north pole
        ir = iz
        ipix1 = 2L*ir*(ir-one)        ;  lowest pixel number in the ring
        ipix2 = ipix1 + 4L*ir - one   ; highest pixel number in the ring
    endif else begin                         ;    south pole
        ir = 4*nside - iz
        ipix1 = npix - 2L*ir*(ir+1L) ;  lowest pixel number in the ring
        ipix2 = ipix1 + 4L*ir - 1   ; highest pixel number in the ring
    endelse
    nr = ir*4L
    kshift = 1
endelse

;     ----------- constructs the pixel list --------------
if (take_all) then begin
    nir    = ipix2 - ipix1 + 1
    listir = lindgen(nir)+ipix1
    if (do_nest) then ring2nest, nside, listir, listir
    return, listir
endif

shift = kshift * .5d0
if keyword_set(conservative) then begin
; conservative : include every intersected pixels, 
; even if pixel CENTER is not in the range [phi_low, phi_hi]
    ip_low = NINT (nr * phi_low / twopi - shift)
    ip_hi  = NINT (nr * phi_hi  / twopi - shift)
    ip_low = ip_low MOD nr      ; in {0,nr-1}
    ip_hi  = ip_hi  MOD nr      ; in {0,nr-1}
endif else begin
; strict : include only pixels whose CENTER is in [phi_low, phi_hi]
    ip_low = ceil (nr * phi_low / TWOPI - shift) 
    ip_hi  = floor(nr * phi_hi  / TWOPI - shift) 
;     if ((ip_low - ip_hi eq 1) and (dphi*nr lt !DPI)) then begin  ; EH, 2004-06-01
    diff = (ip_low - ip_hi) mod nr      ; in {-nr+1,nr-1}
;    if (diff < 0) then diff = diff + nr ; in {0,nr-1} ; EH, 2008-03-29
    if (diff lt 0) then diff = diff + nr ; in {0,nr-1}
    if ((diff eq 1) && (dphi*nr lt !DPI)) then begin
; the interval is so small (and away from pixel center)
; that no pixel is included in it
        nir = 0
        return,-1
    endif
;     ip_low = ip_low < (nr-1) ; EH, 2004-05-28
;     ip_hi  = ip_hi  > 0   
    if (ip_low ge nr) then ip_low = ip_low - nr
    if (ip_hi  lt 0 ) then ip_hi  = ip_hi  + nr
endelse
if (ip_low gt ip_hi) then to_top = 1
ip_low = ip_low + ipix1
ip_hi  = ip_hi  + ipix1

if (to_top) then begin
    nir1 = ipix2 - ip_low + 1
    nir2 = ip_hi - ipix1  + 1
    nir  = nir1 + nir2
    if ((nir1 gt 0) && (nir2 gt 0)) then begin
        listir   = [lindgen(nir2)+ipix1, lindgen(nir1)+ip_low]
    endif else begin
        if nir1 eq 0 then listir   = [lindgen(nir2)+ipix1]
        if nir2 eq 0 then listir   = [lindgen(nir1)+ip_low]
    endelse
    if (do_nest) then ring2nest, nside, listir, listir
endif else begin
    nir = ip_hi - ip_low + 1
    listir = lindgen(nir)+ip_low
    if (do_nest) then ring2nest, nside, listir, listir
endelse

; if min(listir) le 0 then help,nr,ip_low,ip_hi,phi_low,phi_hi,iz
; if min(listir) lt 0 then begin
;     help,take_all,to_top,nir1,nir2,nir,ipix1,ip_low
;     print,listir
;     message,'abort'
; endif
return, listir
end

