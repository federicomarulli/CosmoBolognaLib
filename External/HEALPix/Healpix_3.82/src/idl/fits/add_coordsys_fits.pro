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
pro add_coordsys_fits, info_header, coordsys = coordsys, error = error
;+
; add_coordsys_fits, info_header, coordsys = coordsys, error = error
;
;-

if datatype(coordsys) ne 'STR' then begin
    print,'invalid coordinate system'
    return
endif

coord = strmid(strupcase(strtrim(coordsys,2)),0,1)

if (coord eq 'C' or coord eq 'G' or coord eq 'E') then begin
    sxaddpar, info_header,'COORDSYS',coord, ' Pixelization coordinate system'
    sxaddpar, info_header,'COMMENT',$
      '        G = Galactic, E = ecliptic, C = celestial = equatorial',after='COORDSYS' 
endif else begin
    print,'coordinate system unknown'
    return
endelse

return
end

