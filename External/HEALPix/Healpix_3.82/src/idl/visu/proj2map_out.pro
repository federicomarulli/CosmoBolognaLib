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
function proj2map_out, planmap, offmap=offmap, bad_data = bad_data
;+
; projected_map = proj2map_out(planmap [,Offmap=, Bad_Data=] )
;
; created 2010-05-11, EH
;
;-

if (n_params() eq 0) then begin
    print,'projected_map = proj2map_out(planmap [,Offmap=, Bad_Data=] )'
    return,-1
endif

map_out = planmap

out_map_switch = 1

init_healpix
my_bad_data = defined(bad_data) ? bad_data : !healpix.bad_value

case out_map_switch of
    0: begin
; --------------------------------------------
; on output, bad pixels take value NaN
; --------------------------------------------
;
; deal with pixels=my_bad_data
        bad_pixels = where(abs(map_out/my_bad_data - 1.) lt 1.e-5, nbad_pixels)
        if (nbad_pixels gt 0) then map_out[bad_pixels] = !values.f_nan

; deal with pixels of projected map not matching any sphere pixel
        if defined(offmap)    then map_out[offmap]   = !values.f_nan

    end

    1: begin
; --------------------------------------------
; on output, bad pixels take value my_bad_data
; --------------------------------------------
;
; deal with pixels=NaN
        bad_pixels = where(finite(map_out,/nan), nbad_pixels)
        if (nbad_pixels gt 0) then map_out[bad_pixels] = my_bad_data

; deal with pixels of projected map not matching any sphere pixel
        if defined(offmap)    then map_out[offmap]    = my_bad_data

    end
endcase

return, map_out
end
