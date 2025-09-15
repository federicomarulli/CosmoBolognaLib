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
pro add_ordering_fits,info_header, nested=nested, ring=ring, ordering=ordering,error=error
;+
; add_ordering_fits
;
; add_ordering_fits,info_header, nested=, ring=, ordering=,error=
;
;  adds the ORDERING keyword defined by 'nested' or 'ring' or 'ordering' in
;  the fits header 'info_header'
;  it also adds the card PIXTYPE = 'HEALPIX'
;
;  EH April 1999
;  2005-02-08 : remove creation of header primer if no info_header provided
;  2005-05-25: will add POLCCONV with value COSMO if none is present in info_header, otherwise
;        will leave it unchanged
;-

error=0
; info_header present and not empty
header_flag = 0
if (n_params() eq 1) then begin
    if defined(info_header) then header_flag=1
endif
; ordering given in the header
order_flag = 0
if (header_flag) then order_fits = strtrim(sxpar(info_header,'ORDERING',count=order_flag),2)
; user defined ordering
ring_set = keyword_set(ring)
nest_set = keyword_set(nested)
if defined(ordering) then begin
    case strmid(strupcase(strtrim(ordering,2)),0,4) of
        'RING' : ring_set = 1
        'NEST' : nest_set = 1
        else : begin
            print,'ordering has to be either ''RING'' or ''NESTED'' '
            print,'current ordering value: '''+strtrim(ordering,2)+''''
            error = 1
            return
        end
    endcase
endif
; conflict
if (ring_set*nest_set eq 1 or (order_flag+ring_set+nest_set) eq 0) then begin
    print,' Choose either RING or NESTED in write_fits_map'
    print,' ******* File not written ********'
    error=1
    return
endif

if (ring_set) then ordering='RING'
if (nest_set) then ordering='NESTED'
;; if (header_flag eq 0) then info_header = strarr(1) ; no header, open one

sxaddpar,info_header,'PIXTYPE','HEALPIX ',' HEALPIX pixelisation '
if (order_flag eq 0) then begin ; no Ordering in the header, define it
    sxaddpar,info_header,'ORDERING',ordering,' Pixel ordering scheme, either RING or NESTED'
endif else begin ; one ordering in the header, change it if necessary
    if undefined(ordering) then ordering = order_fits
    if (order_fits ne ordering) then begin
        print,' value of ORDERING changed from '+order_fits+' to '+ordering
        sxaddpar,info_header,'ORDERING',ordering
    endif
endelse

polar = sxpar(info_header,'POLAR')
if (polar) then begin
    polcconv = sxpar(info_header,'POLCCONV',count=npolcc)
    if (npolcc eq 0) then begin
        sxaddpar,info_header,'POLCCONV','COSMO',' Coord. convention for polarisation (COSMO/IAU)',after='POLAR'
;;        sxaddpar,info_header,'COMMENT',' either IAU or COSMO',after='POLCCONV'
    endif
endif

return
end

;--------------------------------
