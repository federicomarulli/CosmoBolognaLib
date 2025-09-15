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

;===========================================================================

pro query_strip, nside, theta1, theta2, listpix, nlist, help=help, nested=nested, inclusive=inclusive
;+
;=======================================================================
;
;    query_strip, Nside, theta1, theta2, Listpix, Nlist, HELP=, NESTED=, INCLUSIVE=
;    --------------
;     nside       = resolution parameter (a power of 2)
;     theta1,theta2  = colatitude in Rad  in [0,Pi]
;      if theta1<theta2 then returns pixels with
;              theta1<colat<theta2
;      if theta2<theta1 then returns
;        0<= colatitude < theta2 or theta1 < colatitude < Pi
;
;     list_pix    = list of pixel lying in the strip
;     nlist       = number of pixels in the list
;     nested  (OPT), :0 by default, the output list is in RING scheme
;                  if set to 1, the output list is in NESTED scheme
;     inclusive (OPT) , :0 by default, only the pixels whose center 
;                       lie in the strip are listed on output
;                  if set to 1, all pixels overlapping the strip are output
;     help (OPT):   prints this documentation header and exits
;
;
;
; v1.0, EH, IAP, 2008-03-28 : adapted from F90 code
; 2012-01-14: systematically returns Listpix=[-1], Nlist=0 in case of problem
;=======================================================================
;-

routine = 'query_strip'
syntax = routine+', Nside, theta1, theta2, Listpix [, Nlist, HELP=, NESTED=, INCLUSIVE=]'
nlist = 0 & listpix = [-1]

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if (n_params() lt 4 or n_params() gt 5) then begin
    print,syntax
    return
endif


npix = nside2npix(nside)
if (npix lt 0) then begin 
    message,"Invalid Nside = "+string(nside)
endif
lnside = long(nside)

if (theta1 lt 0.d0 or theta1 gt !PI or theta2 lt 0.d0 or theta2 gt !PI) then begin
    message,/info," the colatitudes are in RADIAN "
    message,/info," and should lie in [0,Pi] "
    message,/info," current value = ", theta1, theta2
    message," program abort "
endif

if (theta1 lt theta2) then begin
    nstrip = 1
    colrange = [ theta1, theta2 ]
endif else begin
    nstrip = 2
    colrange = [ 0.0d0, theta2, theta1, !DPI ]
endelse

do_inclusive = keyword_set(inclusive)


nlist=-1
for is =0, nstrip-1 do begin
    zu = cos(colrange[2*is  ])  ; upper (=Northern) bound in z
    zd = cos(colrange[2*is+1])  ; lower (=Southern) bound in z
    irmin = ring_num(lnside, zu, shift=-do_inclusive) ; shift up north for inclusive
    irmax = ring_num(lnside, zd, shift= do_inclusive) ; shift down south for inclusive


                                ;     ------------- loop on ring number ---------------------
    for iz = irmin, irmax do begin
        
                                ; z of ring being considered
        zring = ring2z(lnside, iz)

        if ((zring ge zd and zring le zu) or do_inclusive) then begin
                                ;        ------- finds pixels in the ring ---------
            listir=in_ring(lnside, iz, 0.d0, !DPI, nir, nested=nested)
             
            if nlist le 0 then begin
                listpix = listir 
                nlist = nir
            endif else begin
                listpix = [listpix,listir]
                nlist = nlist + nir
            endelse

        endif

    endfor

endfor


return
end

