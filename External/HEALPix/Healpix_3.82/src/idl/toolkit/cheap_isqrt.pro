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
function cheap_isqrt, in
;+
; NAME:
;        CHEAP_ISQRT
;
;
; PURPOSE:
;       Returns exact Floor(sqrt(x)) where x is a (64 bit) integer.
;       The double precision floating point operation is not accurate enough
;       when dealing with 64 bit integers, especially in the vicinity of 
;       perfect squares. 
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;     y = cheap_isqrt(x)
;
;
; INPUTS:
;      x: integer number
;
;
; OUTPUTS:
;      y: 64 bit integer, such that
;         y^2 <= x < (y+1)^2         (1)
;
; PROCEDURE:
;     applies double precision SQRT and then enforces Eq (1) using 
;     integer numbers
;
;
; RESTRICTION:
;     the input must be a positive (signed) 64-bit integer
;
;
; EXAMPLE:
;    in = 2LL^54-1
;    print, cheap_isqrt(in), floor(sqrt(in))
;
;  will return 134217727   1324217728
;  while the exact answer is 2LL^27-1 = 134217727
;
; MODIFICATION HISTORY:
;   Aug 2011; v1.0
;-

din = double(in) ; input integer number has ~19 significant digits (in base 10)
dout = sqrt(din) ; accurate to ~15 digits (base 10) , for ~10 needed
out  = floor(dout, /l64) ; limited accuracy creates round-off error

; integer arithmetics solves round-off error
diff = in - out*out
toohi = where(diff lt 0,     nhi)
toolo = where(diff gt 2*out, nlo)
if (nhi gt 0) then out[toohi] -= 1
if (nlo gt 0) then out[toolo] += 1

;print, out
return, out
end
