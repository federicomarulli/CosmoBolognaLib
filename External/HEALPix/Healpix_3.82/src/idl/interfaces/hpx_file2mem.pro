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
pro hpx_file2mem, tmpfile, variable, cl=cl, map=map, show_cl=show_cl, silent=silent

;; common hpx_xface_com, tmp_head, to_remove

; 2010-02-22: added SILENT keyword

do_map = keyword_set(map)
do_cl  = keyword_set(cl)

if ((do_map + do_cl) ne 1) then begin
    message,'Choose 1 and only 1 among, /cl and /map'
endif


if (size(variable, /tname) ne 'STRING') then begin
                                ; valid array: read file into data
    case 1 of
        do_cl:  fits2cl, variable, tmpfile, show=show_cl, silent=silent
        do_map: begin
            npix = getsize_fits(tmpfile, nmaps=nmaps)
            if (nmaps eq 1) then read_fits_map, tmpfile, variable, silent=silent else read_tqu, tmpfile, variable
        end
    endcase
endif else if (keyword_set(show_cl) && do_cl) then begin
    if (strtrim(tmpfile,2) ne '' && strtrim(tmpfile,2) ne "''") then fits2cl, junk, tmpfile, show=1, silent=silent
endif


return
end
