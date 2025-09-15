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
function npix2nside, npix, error=error
;+
; NAME:
;   npix2nside
;
; PURPOSE:
;   returns the Nside resolution parameter corresponding to the pixel number Npix
;
; CATEGORY:
;   healpix toolkit
;
; CALLING SEQUENCE:
;   nside = npix2nside(npix [, err=err] )
;
; INPUTS:
;   npix : integer, number of healpix pixels over the full sky
;
; OUTPUTS:
;   nside : integer, healpix resolution parameter,
;           nside = sqrt(npix/12)
;
; OPTIONAL OUTPUTS:
;   error : is set to 0 only if the number of pixels does correpond to a Healpix
;     tesselation of the full sphere
;     ie, npix = 12*nside^2 with nside 
;     a power of 2 <= 2^29
;   in any other cases, error = 1
;
; PROCEDURE:
;   trivial
;
; EXAMPLE:
;    nside = npix2nside(786432,err=err)
;    print,nside
;           will return : 256 
;
; MODIFICATION HISTORY:
;
;     v1.0, EH, Caltech, 2000-02-11
;     v1.1, EH, Caltech, 2002-08-16 : uses !Healpix structure
;     v2.1, EH, IAP, 2008-03-17; enabled nside > 8192
;
;-

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

error = 1
fnside = sqrt(npix/12.)

; closest integer
nside = long(round(fnside))

; is npix = 12*nside^2 ?
nnpix = (12LL*nside)*nside
if abs(nnpix-npix) gt .1 then return,-1  ;fnside

; is nside a power of 2 ?
junk = where(nside eq !healpix.nside, count)
if count ne 1 then return,-1  ;nside

; npix and nside are fine
error=0
return,nside
end

