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
pro xy2pix_nest, nside, ix, iy, face_num, ipnest
;+
; xy2pix_nest, nside, ix, iy, face_num,    ipnest
;
; v1.0 June 2006, EH
; Dec 2007, EH, IAP, enabled Nside > 8192
;-

common xy2pix, x2pix, y2pix
if (n_elements(x2pix) eq 0) then init_xy2pix ; initiate x2pix and y2pix

syntax = 'xy2pix_nest, nside, ix, iy, face_num,    ipnest'

if n_elements(nside) gt 1 or n_elements(face_num) gt 1 then begin
    print,syntax
    message,'Nside and face_num should be scalar numbers'
    return
endif

if n_elements(ix) ne n_elements(iy) then begin
    print,syntax
    message,'ix and iy should have the same number of elements'
    return
endif

ix_tmp = ix
iy_tmp = iy
if (nside gt 8192) then begin
    smax = 4
    scale = 1LL
    ipnest = 0LL
    scale_factor = 128LL*128LL
    npface = long64(nside)^2
endif else begin
    smax = 1
    scale = 1L
    ipnest = 0L
    scale_factor = 128L*128L
    npface = long(nside)^2
endelse      

for i=0, smax-1 do begin
    ix_low = ix_tmp and 127     ; last 7 bits
    iy_low = iy_tmp and 127     ; last 7 bits
    ipnest += (x2pix[ix_low] + y2pix[iy_low]) * scale
    scale *= scale_factor
    ix_tmp /= 128               ; truncate out last 7 bits
    iy_tmp /= 128 
endfor
ipnest += (x2pix[ix_tmp] + y2pix[iy_tmp]) * scale
ipnest += face_num* npface      ; in {0, 12*nside**2 - 1}

; ix_low = 0 & ix_tmp = 0 & iy_low = 0 & iy_tmp = 0

return
end
