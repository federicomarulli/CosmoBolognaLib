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
function vector_map, norm, Obs, vector_scale = vector_scale

;
; given an array of vector norms, 
; returns the average norm and divide norms by this average
;
;
sz = size(norm)
npix = n_elements(norm)
normout = make_array(/float, npix, value=-1.)

if defined(Obs) then begin
    if Obs[0] eq -1 then begin
        vector_scale = 0.
        return, normout
    endif
    N_obs = n_elements(Obs)
    N_no_Obs = npix - N_obs
endif else begin
    N_obs = npix
    N_no_obs = 0
endelse

if (N_NO_Obs eq 0) then begin
    if undefined(vector_scale) then begin
;     mean1 = total(norm)/N_obs
        mean1 = median(norm)
        vector_scale = mean1
    endif
    normout      = norm/vector_scale
endif else begin
    if undefined(vector_scale) then begin
;    mean1 = total(norm[Obs])/N_obs
        mean1 = median(norm[Obs])
        vector_scale = mean1
    endif
    normout[Obs] = norm[Obs]/vector_scale
endelse


return, normout
end

