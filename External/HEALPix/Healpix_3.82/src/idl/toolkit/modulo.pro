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

;===============================================================
function modulo, arg1, arg2
;+
;  x = modulo(x1, x2)
;  returns x = x1 - floor(x1/x2) * x2
;
;  while x = x1 MOD x2
;  returns x = x1 - int(x1/x2) * x2
;
;  (they differ when x1 in ]-x2, 0[  )
;
;
;-

n2 = n_elements(arg2)
out = arg1 mod arg2
neg = where(out lt 0, nneg)
if (nneg gt 0) then out[neg] += arg2[neg < (n2-1)]

return, out
end
