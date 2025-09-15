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
pro find_min_max_valid, data, mindata, maxdata, valid=valid, bad_data=bad_data
;+
;
; find_min_max_valid, data, mindata, maxdata, [valid=, bad_data=]
;
;
;  returns minvalue, maxvalue of valid pixels among data
;
;  valid pixels are those with finite value and > bad_data
;
;  if any pixel is unvalid, then valid returns the list of valid pixels
;
;
; v1.0: 2009-10-02
;
;-


nd = n_elements(data)

threshold = keyword_set(bad_data) ? bad_data : !healpix.bad_value

mindata = MIN(data, MAX=maxdata)



if (mindata le threshold || total(finite(data)) ne nd) then begin
    ; some unvalid data: create list of valid pixels
    Valid    = where(data gt threshold and finite(data), N_Valid )
    if (N_valid gt 0) then mindata = MIN(data[Valid], MAX=maxdata)

endif else begin
 
    if defined(Valid) then begin
        Valid = -1
        junk = temporary(Valid)   ; Valid is not defined
    endif

endelse


return
end
