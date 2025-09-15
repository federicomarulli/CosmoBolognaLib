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
pro sixpack2vector, sixpack, vector
;+
; sixpack2vector, sixpack, vector
;
; turns a sixpack raster into a data vector ordered by pixel number
;
;
; uses : pix2xy from CGIS
;
; written by EH, oct-2001
;-


np = n_elements(sixpack)
nxf = sqrt(np/6)             ; side of each face
resq = alog(nxf)/alog(2.)+1  ; resolution parameter

pixel = lindgen(np)
pix2xy, pixel, ix, iy, res = resq, /sixpack
vector = sixpack(ix,iy)



return
end
