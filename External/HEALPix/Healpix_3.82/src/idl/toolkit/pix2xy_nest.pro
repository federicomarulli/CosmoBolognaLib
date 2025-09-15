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
pro pix2xy_nest, nside, ipf, ix, iy
;+
; pix2xy_nest, nside, ipf,      ix, iy
;
; v1.0 June 2006, EH
; Dec 2007, EH, IAP, enabled Nside > 8192
;-

common pix2xy, pix2x, pix2y
if (n_elements(pix2x) eq 0) then init_pix2xy ; initiate pix2x and pix2y

;     finds the x,y on the face (starting from the lowest corner)
;     from the pixel number
ipf_tmp = ipf
scale = 1L
ix = 0L & iy = 0L
smax = (nside le 8192) ? 2 : 5
for i=0, smax-1 do begin
    ip_low = ipf_tmp and 1023   ; last 10 bits
    ix += scale * pix2x[ip_low]
    iy += scale * pix2y[ip_low]
    scale *= 32L
    ipf_tmp /= 1024             ; truncate out last 10 bits
endfor
ix += scale * pix2x[ipf_tmp]
iy += scale * pix2y[ipf_tmp]
;;ip_low = 0 & ipf_tmp = 0        ; free memory


return
end
