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
pro tbfree, stc
;+
;  free heap pointers created by TBINFO
;
;
;  tbfree, structure
;
;  EH, Jan 2008
;-

to_remove = ['TSCAL','TZERO']

tags = tag_names(stc)
for i=0, n_elements(to_remove)-1 do begin
    ; make sure that tag exist and is a pointer before freeing it
    id = where(tags eq to_remove[i], nfound)
    if (nfound gt 0) then begin
        if size(stc.(id), /tname) eq 'POINTER' then ptr_free, stc.(id)
    endif
endfor

return
end
