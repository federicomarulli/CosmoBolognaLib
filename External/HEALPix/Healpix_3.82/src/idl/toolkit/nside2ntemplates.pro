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
function nside2ntemplates, nside, error=error
;+
; ntemplates = nside2ntemplates(nside, error=error)
;
; returns the number of template shapes for a given Nside
;           = Nside * (Nside + 6) /4   for Nside > 1
;           = 2                        for Nside = 1
;
; if nside is not a power of 2 <= 2^29,
; -1 is returned and the error flag is set to 1
;
; MODIFICATION HISTORY:
;
;     v1.0, EH, Caltech, ?
;     v1.1, EH, IAP,     2000-03-18: enabled nside > 8192
;     v1.2, EH, IAP,     2018-06-06: accept non-scalar Nside
;-

npix = nside2npix(nside,err=error)

if (error ne 0) then begin
    bad = where(npix lt 0, nbad)
    if nbad eq 0 then bad = 0
    message,' WARNING: invalid nside: '+string(nside[bad]),/info
    return, -1
endif

if (max(nside) gt 8192) then begin
    lnside = long64(nside)
    nt = 1LL + lnside * (lnside + 6LL)
    ntemplate = long64( nt / 4LL )
endif else begin
    lnside = long(nside)
    nt = 1L + lnside * (lnside + 6L)
    ntemplate = long( nt / 4L )
endelse


return, ntemplate
end
