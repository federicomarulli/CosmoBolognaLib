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
function hpx_mem2file, variable, input=input, output=output, beam=beam, cl=cl, map=map, ring=ring, nested=nested, ordering=ordering
;+
; HPX_MEM2FILE
;   read/write data from/to disk to/from memory
;
;
;
;  version 1.0: ??
;  2008-10-14: explicit message if ordering of input map is not provided
;
;-
common hpx_xface_com, tmp_head, to_remove, uid

if undefined(uid) then uid = 0
if (undefined(tmp_head)) then begin
    message,'tmp directory not defined'
endif

do_in  = keyword_set(input)
do_out = keyword_set(output)
do_map = keyword_set(map)
do_cl  = keyword_set(cl)
do_bl  = keyword_set(beam)

if (do_in eq do_out) then begin
    message,'Choose either /input or /output'
endif
if (do_in && (do_map + do_cl + do_bl) ne 1) then begin
    message,'Choose 1 and only 1 among, /beam, /cl and /map'
endif

autofile = string(uid,form='(i2.2)')+'tmp.fits'
uid = (uid + 1) mod 100

if (do_out) then begin
    ; output data
    if (size(variable, /tname) ne 'STRING') then begin
        ; valid array : define file to receive data
        tmpfile = tmp_head + autofile
        to_remove = [to_remove, tmpfile] ; add this file to those to remove
    endif else begin
        ; file: leave it as is
        tmpfile = variable
    endelse
endif else begin
    ; input data
    if (size(variable, /tname) ne 'STRING' && n_elements(variable) gt 1) then begin
        ; valid array -> write into file
        tmpfile = tmp_head + autofile
        to_remove = [to_remove, tmpfile] ; add this file to those to remove
        case 1 of
            do_map: begin
                nmaps = n_elements(variable)/n_elements(variable[*,0])
                if (undefined(ring) && undefined(nested) && undefined(ordering)) then begin
                    message_patch,'Set either /RING or /NESTED or ORDERING=  when dealing with online maps',level=-1
                endif
                error = 0
                if (nmaps eq 1) then begin
                    write_fits_map, tmpfile, variable, ring=ring, nested=nested, ordering=ordering, error=error
                endif else begin
                    write_tqu, tmpfile, variable, ring=ring, nested=nested, ordering=ordering, error=error
                endelse
                if (error ne 0) then message,'Map not written.'
            end
            do_cl:  cl2fits, variable, tmpfile
            do_bl:  bl2fits, variable, tmpfile
            else: message,'Invalid usage.'
        endcase
    endif else begin
        ; file : leave it as is
        tmpfile = variable
    endelse
endelse

return, tmpfile
end
