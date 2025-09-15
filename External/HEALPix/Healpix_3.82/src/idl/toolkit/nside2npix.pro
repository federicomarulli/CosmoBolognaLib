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
function nside2npix, nside, error=error, help=help
;+
; npix = nside2npix(nside, ERROR=, HELP=)
;
; returns npix = 12*nside*nside
; number of pixels on a Healpix map of resolution nside
;
; if nside is not a power of 2 <= 2^29,
; -1 is returned and the error flag is set to 1
;
; MODIFICATION HISTORY:
;
;     v1.0, EH, Caltech, 2000-02-11
;     v1.1, EH, Caltech, 2002-08-16 : uses !Healpix structure
;     v2.1, EH, IAP,     2006-10-16 : enabled nside > 8192
;     v2.2, EH, IAP,     2011-04-07: returns syntax in case of improper use
;     v3.0, EH, IAP,     2018-05-22 : accepts non-scalar nside, added /HELP
;-

routine = 'nside2npix'
syntax = 'Npix = '+routine+' (Nside [,error=])'

if keyword_set(help) then begin
    error = 0
    doc_library,routine
    return, -1
endif

error = 1
if N_params() ne 1 then begin
    print,syntax
    return,-1
endif

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

ni        = n_elements(nside)
nside_max = max(nside)
tag       = replicate(1, ni)
for i=0,ni-1 do begin
    junk = where(nside[i] eq !healpix.nside, count)
    ;if count ne 1 then return,-1
    if count ne 1 then tag[i] = -1
endfor

; compute npix
if (nside_max gt 8192) then begin
    npix = 12LL* long64(nside)^2
endif else begin
    npix = 12L * long(nside)^2
endelse

; flag bad nside (and npix)
if (min(tag) le 1) then begin
    npix = npix*tag > (-1)
endif

; if input is scalar, return a scalar
if size(nside,/n_dim) eq 0 then npix = npix[0]

if (min(npix) gt 0) then error=0

return, npix
end

