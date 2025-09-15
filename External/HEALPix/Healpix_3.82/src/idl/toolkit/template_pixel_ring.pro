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
pro template_pixel_ring, nside, ipix, itplt, reflexion
;+
; NAME:
;         template_pixel_ring
;
; PURPOSE:
;         Returns template-pixel index corresponding to each (RING ordered) pixel provided
;         The reflexion (E-W, or N-W or both) to apply to the pixel to match the
;         template is optionally returned
;
; CATEGORY:
;
; CALLING SEQUENCE:
;         TEMPLATE_PIXEL_RING, Nside, Ipix, Template [,Reflexion]
;
; INPUTS:
;         Nside : scalar,           integer, Healpix resolution parameter
;         Ipix  : vector or scalar, integer, list of ring ordered pixel indices
;
; OUTPUTS:
;         Template : vector or scalar, integer, template-pixel index
;         corresponding to Ipix (see Notes below). Same dimension as Ipix.
;
; OPTIONAL OUTPUTS:
;         Reflexion : vector or scalar, integer in {0,3} (see Notes below). 
;         Same dimension as Ipix.
;
; COMMON BLOCKS:
;         none
;
; SIDE EFFECTS:
;         none
;
; RESTRICTIONS:
;         ??
;
; PROCEDURE:
;
; EXAMPLE:
;
;
; RELATED ROUTINES:
;         NEST2RING, RING2NEST
;         SAME_SHAPE_PIXELS_NEST, SAME_SHAPE_PIXELS_RING
;         TEMPLATE_PIXEL_NEST,    TEMPLATE_PIXEL_RING
;         NSIDE2NTEMPLATES
;
; MODIFICATION HISTORY:
;         2004-09, EH, v1.0
;    Dec 2007, EH,  IAP, enabled nside > 8192
;    Aug 2008, EH, IAP: issues warning if ipix is not of integer type
;
; Notes on pixel template
;
;         For a given resolution Nside, there are 12*Nside*Nside pixels on the
;         sphere, whose shape belong to one of the (Nside+6)*Nside/4 templates.
;         (2 in the case Nside=1).
;         Any pixel can be matched to a single of these templates by a combination of
;          - a rotation around the polar axis with 
;          - reflexion(s) around a meridian and/or the equator.
;         The reflexion number returned by this routine for each pixel is
;         Reflexion         Operation to do to match template
;           0                 rotation only
;           1                 rotation + East-West swap
;           2                 rotation + North-South swap
;           3                 rotation + East-West + North-South swaps
;
;         The template pixels are all located in the Northern Hemisphere, or on the
;         Equator.
;         They are chosen to have
;             z >= 2/3,      0< phi <= Pi/2               [Nside*(Nside+2)/4]
;             2/3 > z >= 0,  phi=0, or phi=Pi/(4*Nside)   [Nside]
;         They are numbered from 0, starting at the North Pole, with the index
;         increasing in phi, and then increasing for decreasing z.
;
;-

routine = 'template_pixel_ring'
nparams = n_params()
if nparams lt 3 or nparams gt 4 then begin
    print,' syntax: '+routine+', Nside, Ipix, Template [,Reflexion]'
    message,'Abort'
    return
endif
do_sym = (nparams eq 4)

if (N_ELEMENTS(nside) GT 1) then message,'Nside should be a scalar in '+routine
npix = nside2npix(nside, error = error)
if (error ne 0) then message,'Invalid Nside '+string(nside)

assert_pixindex_type, ipix[0], /warning ; warning if ipix is not integer
min_pix = MIN(ipix, MAX = max_pix)
IF (min_pix LT 0) THEN BEGIN
    PRINT,'pixel index : ',min_pix,FORMAT='(A,I18)'
    PRINT,'is out of range : ',0,npix-1,FORMAT='(A,I2,I18)'
    message,'Abort'
ENDIF
IF (max_pix GT npix-1) THEN BEGIN
    PRINT,'pixel index : ',max_pix,FORMAT='(A,I18)'
    PRINT,'is out of range : ',0,npix-1,FORMAT='(A,I2,I18)'
    message,'Abort'
ENDIF

lnside = long(nside[0]) ; 2018-06-06
np = n_elements(ipix)

if (lnside gt 8192) then begin
    ncap = 2LL*lnside*(lnside+1LL)
    ns4  = 4LL*lnside
    n0   = (lnside*(lnside+2LL))/4
    one  = 1LL
    itplt = lon64arr(np)
    l64 = 1
endif else begin
    ncap = 2*lnside*(lnside+1)
    ns4  = 4*lnside
    n0   = (lnside*(lnside+2))/4 > 1
    one  = 1L
    itplt = lonarr(np)
    l64 = 0
endelse

if (do_sym) then begin
    reflexion = intarr(np)
endif

; North polar cap
knorth = where(ipix lt ncap, nnorth)
if (nnorth gt 0) then begin
    ip = ROUND(ipix[knorth], L64=l64) + one
    iring = LONG( SQRT( ip/2.d0 - SQRT(ip/2) ) ) + one ; counted from NORTH pole, starting at 1
    ifi  = (ip - 1L - 2L*iring*(iring-one)) MOD iring ; in [0,iring-1], upwards
    ifd  = iring - 1L - ifi                          ; in [0,iring-1], downwards
    iff  = ifi < ifd                                ; in [0,(iring-1)/2]
    itplt[knorth] = (iring*iring)/4 + iff
    if (do_sym) then begin
        reflexion[knorth] = (ifd lt ifi)             ; East-West swap
    endif
    ip = 0 & iring = 0 & ifi = 0 & ifd = 0 & iff = 0 ; free memory
endif
knorth = 0 ; free memory

; North-Equatorial region
kneq   = where(ipix ge ncap and ipix lt (npix+ns4)/2, nneq)
if (nneq gt 0) then begin
    ip = ROUND(ipix[kneq], L64=l64)
    itplt[kneq] = n0  + long((ip - ncap)/ns4)
    if (do_sym) then begin
        reflexion[kneq] = 0      ; no swap
    endif
    ip = 0
endif
kneq = 0 ; free memory

; South-Equatorial region
kseq   = where(ipix ge (npix+ns4)/2 and ipix lt (npix-ncap), nseq)
if (nseq gt 0) then begin
    ip = ROUND(ipix[kseq], L64=l64)
    itplt[kseq] = n0 + long((npix - one - ip - ncap)/ns4)
    if (do_sym) then begin
        reflexion[kseq] = 2      ; North-South swap
    endif
    ip = 0
endif
kseq = 0 ; free memory

; South polar cap
ksouth = where(ipix ge (npix-ncap),  nsouth)
if (nsouth gt 0) then begin
    ip =  npix - ROUND(ipix[ksouth], L64=l64)
    iring = LONG( SQRT( ip/2.d0 - SQRT(ip/2) ) ) + one ; counted from SOUTH pole, starting at 1
    ifi  = (2L*iring*(iring+one) - ip) MOD iring     ; in [0,iring-1], upwards
    ifd  = iring - 1L - ifi                          ; in [0,iring-1], downwards
    iff  = ifi < ifd                                ; in [0,(iring-1)/2]
    itplt[ksouth] = (iring*iring)/4 + iff
    if (do_sym) then begin
        reflexion[ksouth] = (ifd lt ifi) + 2       ; East-West swap + North-South swap
    endif
    ip = 0 & iring = 0 & ifi = 0 & ifd = 0 & iff = 0 ; free memory
endif
ksouth = 0 ; free memory


return
end
