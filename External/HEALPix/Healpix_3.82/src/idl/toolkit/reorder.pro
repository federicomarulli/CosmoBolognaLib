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
function reorder, map_in, in=in, out= out, n2r= n2r, r2n=r2n, help=help
;+
; NAME:
;   reorder
;
; PURPOSE:
;   change the ordering of a full sky Healpix map from Nested/Ring to
;   Ring/Nested
; CATEGORY:
;   Healpix pixel toolkit
;
; CALLING SEQUENCE:
;   map_out = reorder(map_in, [in=in, out=out, /n2r, /r2n, /help])
;
; INPUTS:
;   map_in : a full sky Healpix map, can be of any type 
;  (float, integer, double) and any dimension
;
; OPTIONAL INPUTS:
;   none
;
; KEYWORD PARAMETERS:
;   help- if set, print out this Help header and exits
;   in- is either 'RING' or 'NESTED' : input map ordering
;   out- is either 'RING' or 'NESTED' : output map ordering
;   r2n- if set, equivalent to in='RING', out='NESTED'
;   n2r- if set, equivalent to in='NESTED', out='RING'
;
; OUTPUTS:
;   map_out : a reordered full sky Healpix map,
;   same size and same type of map_in
;
; OPTIONAL OUTPUTS:
;   none
;
; COMMON BLOCKS:
;   none
;
; SIDE EFFECTS:
;   none
;
; RESTRICTIONS:
;   none
;
; PROCEDURE:
;
; EXAMPLES:
;   map_nest = reorder(map_ring,in='ring',out='nest')
;   map_nest = reorder(map_ring,/r2n)
;
; MODIFICATION HISTORY:
;    April 1999, EH, Caltech
;    Jan   2000, EH,  improved documentation header 
;    Oct   2004, EH,  added extra dimension, added R2N  and N2R
;    Mars  2012, EH, IAP: added HELP keyword
;    May   2017, EH, IAP: replaced obsolete DATATYPE() with IDL's SIZE(/TYPE)
;-

 ok = (defined(in) && defined(out)) + keyword_set(r2n) + keyword_set(n2r)
 if keyword_set(help) then begin
     doc_library,'reorder'
     return,-1
 endif

 if (N_params() ne 1) then begin
     print,' map_out = reorder(map_in, [In=in, Out=out,] [/R2N, /N2R], [/HELP])'
     return,0
 endif

 if (ok ne 1) then begin
     message,/info,'should be either /R2N *or* /N2R *or* (IN=order_in, OUT=order_out)'
     message,'Abort'
 endif

 if (defined(in)) then begin
     kin  = strmid(strupcase(strtrim(in ,2)),0,4)
     kout = strmid(strupcase(strtrim(out,2)),0,4)
 endif
 if (keyword_set(r2n)) then begin
     kin = 'RING' & kout = 'NEST'
 endif
 if (keyword_set(n2r)) then begin
     kout = 'RING' & kin = 'NEST'
 endif


 if (kin ne 'RING' and kin ne 'NEST') then begin
     print,' In has to be either ''RING'' or ''NESTED'' '
     print,' map_out = reorder(map_in, In=in, Out=out)'
     return,0
 endif

 if (kout ne 'RING' and kout ne 'NEST') then begin
     print,' Out has to be either ''RING'' or ''NESTED'' '
     print,' map_out = reorder(map_in, In=in, Out=out)'
     return,0
 endif

 if (kin eq kout) then return, map_in

 np   = n_elements(map_in[*,0])
 ncol = n_elements(map_in[0,*])
 nside = npix2nside(np,err=err)

 if (err ne 0) then begin
     print,' map_in should be a full sky Healpix map'
     print,' pixel #  = ',np
     return,0
 endif

 npf = np/12
 map_out = make_array(np, ncol, type=size(/type,map_in))
 if (kin eq 'RING') then begin  ; ring -> nest
     for j=0L, 11 do begin
         ipn = lindgen(npf) + j*npf
         nest2ring, nside, ipn, ipr
         for ic=0,ncol-1 do map_out[ipn,ic] = map_in[ipr,ic]
     endfor
 endif
 if (kin eq 'NEST') then begin  ; nest -> ring
     for j=0L, 11 do begin
         ipn = lindgen(npf) + j*npf
         nest2ring, nside, ipn, ipr
         for ic=0,ncol-1 do map_out[ipr,ic] = map_in[ipn,ic]
     endfor
 endif


return, map_out
end

