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
function gaussbeam, fwhm, lmax, dim, help=help
;+
; NAME:
;   gaussbeam
;
; PURPOSE:
;    returns the spherical transform of a gaussian beam
;    from l=0 to lmax with a fwhm in arcmin
;
; CATEGORY:
;
; CALLING SEQUENCE:
;   result = gaussbeam(fwhm, lmax, [dim, HELP=])
; 
; INPUTS:
;   fwhm : scalar, FWHM in arcmin
;   lmax : scalar integer, lmax
;
; OPTIONAL INPUT:
;   dim  : scalar integer, dimension of the output
;     if absent or set to 0 or 1,
;          the output has size (lmax+1) and is the temperature beam
;     if set to >= 2 <= 4,
;          the output has size (lmax+1,dim)
;          and contains in that order the TEMPERATURE beam,
;          the GRAD/ELECTRIC polarisation beam
;          the CURL/MAGNETIC polarisation beam
;          the TEMPERATURE*GRAD beam
;
;  OPTIONAL KEYWORD:
;    HELP:  if set, a documentation header is printed out
;
; OUTPUTS:
;   result contains the window for l in [0,lmax]
;   corresponding to a gaussian beam of FWHM 'fwhm'
;
; PROCEDURE:
;    result = exp(-0.5*l*(l+1)*sigma^2)
;    with sigma = fwhm/60.*!DtoR / sqrt(8.*alog(2.))
;    and l = lindgen(lmax+1)
;
;    for Grad and Curl, there is an extra factor exp(2*sigma^2)
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;     version 1.0, EH, Caltech, 11-1999
;     version 1.1, EH, Caltech, 03-2000 : double precision
;     version 1.2, EH, Caltech, 08-2002 : added polarisation
;     version 1.3, EH, IAP,     01-2013: double precision Deg to Rad conversion
;     version 1.4, EH, IAP,     08-2016: added HELP
;-

routine = 'gaussbeam'
Syntax = 'Syntax : beam = '+routine+'(fwhm, lmax [, dim, HELP=])'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if (n_params() lt 2 or n_params() gt 3) then begin
    message,syntax,/info
    message,'For extended help, type  x = gaussbeam(/help)'
endif

; translates fwhm (arcmin) into sigma = rms (radians)
;sigma = fwhm/60.d0*!DtoR / sqrt(8.d0*alog(2.d0))
sigma = fwhm/60.d0*!DPI/180.d0 / sqrt(8.d0*alog(2.d0))

; l is in [0,lmax]
l = dindgen(lmax+1)

; computes the gaussian window
sig2 = sigma*sigma
g = exp(-0.5d0*l*(l+1.d0)*sig2)

; compute polarisation beam if necessary
; (cf, eg. Challinor et al, 2000)
if (n_params() eq 3) then begin
    ndim = nint(dim)

    if (ndim lt 0 or ndim gt 4) then begin
        message,syntax,/info
        message,' Invalid dim : '+string(dim)+' (should be in [0, 4])'
    endif

    if (ndim eq 0) then begin
        gout = g
    endif else begin
        factor_pol = exp( [0.d0, 2.d0*sig2, 2.d0*sig2, sig2]  )
        gout = g # factor_pol[0:dim-1]
    endelse

endif else begin
    gout = g
endelse

return,gout

end

