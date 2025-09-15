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
;+
; bit_manipulation:
;
; defines the function
;  iout = swapLSBMSB(iin, [L64=])
;  iout = invLSBMSB(iin, [L64=])
;  iout = invswapLSBMSB(iin, [L64=])
;  iout = invLSB(iin, [L64=])
;  iout = invMSB(iin, [L64=])
;
;  2006-2007: first version
;  2007-12-28: enabled nside > 8192
;-
function swapLSBMSB, i, l64=l64
; Returns i with even and odd bit positions interchanged.

if keyword_set(l64) then begin
; nside <= 2^28
    oddbits = 96076792050570581LL   ; 2^0 + 2^2 + 2^4+..+2^56
    evenbits=192153584101141162LL   ; 2^1 + 2^3 + 2^4+..+2^57

    li = long64(i)
    swapLSBMSB = (li AND evenbits)/2 + (li AND oddbits)*2

endif else begin
; nside <= 2^13
    oddbits = 89478485L         ; 2^0 + 2^2 + 2^4+..+2^26
    evenbits=178956970L         ; 2^1 + 2^3 + 2^4+..+2^27

    li = long(i)
    swapLSBMSB = (li AND evenbits)/2 + (li AND oddbits)*2
endelse

return, swapLSBMSB
end

; -------

function invLSBMSB, i, l64=l64
; Returns NOT(i)

if keyword_set(l64) then begin
    invLSBMSB = NOT (long64(i))
endif else begin
    invLSBMSB = NOT (long(i))
endelse

return, invLSBMSB
end

; -------

function invswapLSBMSB, i, l64=l64
; Returns NOT(i) with even and odd bit positions interchanged.

invswapLSBMSB = NOT swapLSBMSB(i, l64=l64)

return, invswapLSBMSB
end

; -------

function invLSB, i, l64=l64
; Returns i with odd (1,3,5,...) bits inverted.

if keyword_set(l64) then begin
    oddbits = 96076792050570581LL ; 2^0 + 2^2 + 2^4+..+2^56
    invLSB = (long64(i) XOR oddbits)
endif else begin
    oddbits = 89478485L
    invLSB = (long(i) XOR oddbits)
endelse

return, invLSB
end

; -------

function invMSB, i, l64=l64
; Returns i with even (0,2,4,...) bits inverted.

if keyword_set(l64) then begin
    evenbits=192153584101141162LL ; 2^1 + 2^3 + 2^4+..+2^57
    invMSB = (long64(i) XOR evenbits)
endif else begin
    evenbits=178956970L
    invMSB = (long(i) XOR evenbits)
endelse

return, invMSB
end

; ;--------------------

; pro bit_manipulation

; x =  swapLSBMSB(1)
; x =  invLSBMSB(1)
; x =  invswapLSBMSB(1)
; x =  invLSB(1)
; x =  invMSB(1)

; return
; end

;=====================================
pro neighbours_nest, nside, ipix, list, nneigh
;+
; NAME: 
;          neighbours_nest
;
; PURPOSE: 
;          find nearest neighbors pixels in Nested scheme
;
; CATEGORY:
;
; CALLING SEQUENCE: 
;          neighbours_nest, nside, ipix, list [, nneigh]
;
; INPUTS:
;          nside: Healpix resolution parameter
;          ipix: pixel index in [0, 12*Nside*Nside-1]
;
; OUTPUTS:
;          list: list of neighbours of size nneigh (see below)
;
;       The neighbours are ordered in the following way:
;       First pixel is the one to the south (the one west of the south
;        direction is taken
;       for the pixels which don't have a southern neighbour). From
;       then on the neighbours are ordered in the clockwise direction 
;       (from outside the sphere)
;       about the pixel with number ipix.
;
; OPTIONAL OUTPUTS:
;          nneigh: number of neighbours: usually 8, sometimes 7 (for 8 pixels)
;             or 6 (for Nside=1)
;
; COMMON BLOCKS:
;          xy2pix
;
; SIDE EFFECTS:
;         runs init_xy2pix (defining common block xy2pix)
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;         June 2006: v1.0 adapted from B. Wandelt neighbours_nest F90 routine by E.H.@IAP
;         Dec 2007: enabled Nside > 8192
;         July 2008: merge bit manipulation functions in same file
;-

routine = 'NEIGHBOURS_NEST'
syntax = 'syntax: '+routine+', Nside, Ipix, List [, Nneigh]'

if n_params() lt 3 then begin
    print,syntax
    list = [-1]
    nneigh = -1
    return
endif

if (n_elements(nside) gt 1) then message,'Nside should be scalar in '+routine
npix = nside2npix(nside, error = error)
if (error ne 0) then begin
    message, /info, syntax
    message, 'Invalid Nside: '+string(nside)
endif

if (ipix ge npix or ipix lt 0 or n_elements(ipix) gt 1) then begin
    message,/info,syntax
    message,/info,'Ipix should be scalar in [0, ',npix-1,']'
endif
ipix0 = ipix[0]

if (nside eq 1) then begin
    nneigh = 6
    case ipix0 of
        0: list =   [  8, 4, 3, 2, 1, 5 ]
        1: list =   [  9, 5, 0, 3, 2, 6 ]
        2: list =   [ 10, 6, 1, 0, 3, 7 ]
        3: list =   [ 11, 7, 2, 1, 0, 4 ]
        4: list =   [ 11, 7, 3, 0, 5, 8 ]
        5: list =   [  8, 4, 0, 1, 6, 9 ]
        6: list =   [  9, 5, 1, 2, 7,10 ]
        7: list =   [ 10, 6, 2, 3, 4,11 ]
        8: list =   [ 10,11, 4, 0, 5, 9 ]
        9: list =   [ 11, 8, 5, 1, 6,10 ]
        10: list =  [  8, 9, 6, 2, 7,11 ]
        11: list =  [  9,10, 7, 3, 4, 8 ]
    endcase
    return
endif

common xy2pix, x2pix, y2pix
if (n_elements(x2pix) eq 0) then init_xy2pix ; initiate x2pix and y2pix
;resolve_routine,'bit_manipulation',/compile_full_file ; initiate bit manipulation routines (swap*, inv*, ...)
;bit_manipulation

if (nside gt 8192) then begin
    l64 = 1
    nsidesq = npix / 12LL
    face_num=long64(ipix0)/nsidesq
endif else begin
    l64 = 0
    nsidesq = npix / 12L
    face_num=long(ipix0)/nsidesq
endelse
local_magic1=(nsidesq-1)/3
local_magic2=2*local_magic1

ipf = ipix0 mod nsidesq

pix2xy_nest, nside,ipf,ix,iy
ixm=ix-1
ixp=ix+1
iym=iy-1
iyp=iy+1

nneigh=8                  ;Except in special cases below
icase = 0 ;general case
                                ;     Exclude corners
if (ipf eq local_magic2)     then icase = 5 ;WestCorner
if (ipf eq (nsidesq-1))      then icase = 6 ;NorthCorner
if (ipf eq 0)                then icase = 7 ;SouthCorner
if (ipf eq local_magic1)     then icase = 8 ;EastCorner
if (icase ne 0) then goto, skip
                                ;     Detect edges
if ((ipf and local_magic1) eq local_magic1) then icase = 1 ;NorthEast
if ((ipf and local_magic1) eq 0)            then icase = 2 ;SouthWest
if ((ipf and local_magic2) eq local_magic2) then icase = 3 ;NorthWest
if ((ipf and local_magic2) eq 0)            then icase = 4 ;SouthEast

if (icase eq 0) then begin
    ; inside a face
    vx = [ixm, ixm, ixm, ix , ixp, ixp, ixp, ix]
    vy = [iym, iy , iyp, iyp, iyp, iy , iym, iym]
    xy2pix_nest, nside, vx, vy, face_num, list
    return
endif

skip:
ia= face_num/4                  ; in {0,2}
ib= face_num mod 4          ; in {0,3}
ibp=(ib+1)   mod 4
ibm=(ib+4-1) mod 4
ib2=(ib+2)   mod 4

if (ia eq 0) then begin          ;North Pole region
    case icase of
        1: begin                ;NorthEast edge
            other_face=0+ibp
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ixm, iym, face_num, n_1
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ixm, iyp, face_num, n_3
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            ipo=(swapLSBMSB(ipf, l64=l64) mod nsidesq) ;East-West flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo+1 , iyo, other_face, n_5
            n_6=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo-1, iyo, other_face, n_7
        end
        2:begin                 ;SouthWest edge
            other_face=4+ib
            ipo=(invLSB(ipf, l64=l64) mod nsidesq) ;SW-NE flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo, iyo-1, other_face, n_1
            n_2=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo, iyo+1, other_face, n_3
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            xy2pix_nest, nside, ixp, iym, face_num, n_7
            xy2pix_nest, nside, ixp, iy , face_num, n_6
            xy2pix_nest, nside, ixp, iyp, face_num, n_5
        end
        3:begin                 ;NorthWest edge
            other_face=0+ibm
            ipo=(swapLSBMSB(ipf, l64=l64) mod nsidesq) ;East-West flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo, iyo-1, other_face, n_3
            n_4=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo, iyo+1, other_face, n_5
            xy2pix_nest, nside, ixm, iym, face_num, n_1
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ixp, iym, face_num, n_7
            xy2pix_nest, nside, ixp, iy , face_num, n_6
        end
        4:begin                 ;SouthEast edge
            other_face=4+ibp
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ixm, iyp, face_num, n_3
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            xy2pix_nest, nside, ixp, iyp, face_num, n_5
            xy2pix_nest, nside, ixp, iy , face_num, n_6
            ipo=(invMSB(ipf, l64=l64) mod nsidesq) ;SE-NW flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo+1, iyo, other_face, n_7
            n_8=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo-1, iyo, other_face, n_1
        end
        5:begin                 ;West corner
            nneigh=7
            other_face=4+ib
            n_2=other_face*nsidesq+nsidesq-1
            n_1=n_2-2
            other_face=0+ibm
            n_3=other_face*nsidesq+local_magic1
            n_4=n_3+2
            n_5=ipix0+1
            n_6=ipix0-1
            n_7=ipix0-2
        end
        6:begin                 ;North corner
            n_1=ipix0-3
            n_2=ipix0-1
            n_8=ipix0-2
            other_face=0+ibm
            n_4=other_face*nsidesq+nsidesq-1
            n_3=n_4-2
            other_face=0+ib2
            n_5=other_face*nsidesq+nsidesq-1
            other_face=0+ibp
            n_6=other_face*nsidesq+nsidesq-1
            n_7=n_6-1
        end
        7:begin                 ;South corner
            other_face=8+ib
            n_1=other_face*nsidesq+nsidesq-1
            other_face=4+ib
            n_2=other_face*nsidesq+local_magic1
            n_3=n_2+2
            n_4=ipix0+2
            n_5=ipix0+3
            n_6=ipix0+1
            other_face=4+ibp
            n_8=other_face*nsidesq+local_magic2
            n_7=n_8+1
        end
        8:begin                 ;East corner
            nneigh=7
            n_2=ipix0-1
            n_3=ipix0+1
            n_4=ipix0+2
            other_face=0+ibp
            n_6=other_face*nsidesq+local_magic2
            n_5=n_6+1
            other_face=4+ibp
            n_7=other_face*nsidesq+nsidesq-1
            n_1=n_7-1
        end
    endcase                     ; north
endif
if (ia eq 1) then begin          ;Equatorial region
    case(icase) of
        1:begin                 ;NorthEast edge
            other_face=0+ib
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ixm, iym, face_num, n_1
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ixm, iyp, face_num, n_3
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            ipo=(invLSB(ipf, l64=l64) mod nsidesq) ;NE-SW flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo , iyo+1, other_face, n_5
            n_6=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo, iyo-1, other_face, n_7
        end
        2:begin                 ;SouthWest edge
            other_face=8+ibm
            ipo=(invLSB(ipf, l64=l64) mod nsidesq) ;SW-NE flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo, iyo-1, other_face, n_1
            n_2=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo, iyo+1, other_face, n_3
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            xy2pix_nest, nside, ixp, iym, face_num, n_7
            xy2pix_nest, nside, ixp, iy , face_num, n_6
            xy2pix_nest, nside, ixp, iyp, face_num, n_5
        end
        3:begin                 ;NorthWest edge
            other_face=0+ibm
            ipo=(invMSB(ipf, l64=l64) mod nsidesq) ;NW-SE flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo-1, iyo, other_face, n_3
            n_4=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo+1, iyo, other_face, n_5
            xy2pix_nest, nside, ixm, iym, face_num, n_1
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ixp, iym, face_num, n_7
            xy2pix_nest, nside, ixp, iy , face_num, n_6
        end
        4:begin                 ;SouthEast edge
            other_face=8+ib
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ixm, iyp, face_num, n_3
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            xy2pix_nest, nside, ixp, iyp, face_num, n_5
            xy2pix_nest, nside, ixp, iy , face_num, n_6
            ipo=(invMSB(ipf, l64=l64) mod nsidesq) ;SE-NW flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo+1, iyo, other_face, n_7
            n_8=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo-1, iyo, other_face, n_1
        end
        5:begin                 ;West corner
            other_face=8+ibm
            n_2=other_face*nsidesq+nsidesq-1
            n_1=n_2-2
            other_face=4+ibm
            n_3=other_face*nsidesq+local_magic1
            other_face=0+ibm
            n_4=other_face*nsidesq
            n_5=n_4+1
            n_6=ipix0+1
            n_7=ipix0-1
            n_8=ipix0-2
        end
        6:begin                 ;North corner
            nneigh=7
            n_1=ipix0-3
            n_2=ipix0-1
            other_face=0+ibm
            n_4=other_face*nsidesq+local_magic1
            n_3=n_4-1
            other_face=0+ib
            n_5=other_face*nsidesq+local_magic2
            n_6=n_5-2
            n_7=ipix0-2
        end
        7:begin                 ;South corner
            nneigh=7
            other_face=8+ibm
            n_1=other_face*nsidesq+local_magic1
            n_2=n_1+2
            n_3=ipix0+2
            n_4=ipix0+3
            n_5=ipix0+1
            other_face=8+ib
            n_7=other_face*nsidesq+local_magic2
            n_6=n_7+1
        end
        8:begin                 ;East corner
            other_face=8+ib
            n_8=other_face*nsidesq+nsidesq-1
            n_1=n_8-1
            n_2=ipix0-1
            n_3=ipix0+1
            n_4=ipix0+2
            other_face=0+ib
            n_6=other_face*nsidesq
            n_5=n_6+2
            other_face=4+ibp
            n_7=other_face*nsidesq+local_magic2
        end
    endcase                     ; equator
endif
if (ia eq 2) then begin         ;South Pole region
    case(icase) of
        1:begin                 ;NorthEast edge
            other_face=4+ibp
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ixm, iym, face_num, n_1
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ixm, iyp, face_num, n_3
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            ipo=(invLSB(ipf, l64=l64) mod nsidesq) ;NE-SW flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo , iyo+1, other_face, n_5
            n_6=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo, iyo-1, other_face, n_7
        end
        2:begin                 ;SouthWest edge
            other_face=8+ibm
            ipo=(swapLSBMSB(ipf, l64=l64) mod nsidesq) ;W-E flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo-1, iyo, other_face, n_1
            n_2=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo+1, iyo, other_face, n_3
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            xy2pix_nest, nside, ixp, iym, face_num, n_7
            xy2pix_nest, nside, ixp, iy , face_num, n_6
            xy2pix_nest, nside, ixp, iyp, face_num, n_5
        end
        3:begin                 ;NorthWest edge
            other_face=4+ib
            ipo=(invMSB(ipf, l64=l64) mod nsidesq) ;NW-SE flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo-1, iyo, other_face, n_3
            n_4=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo+1, iyo, other_face, n_5
            xy2pix_nest, nside, ixm, iym, face_num, n_1
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ix , iym, face_num, n_8
            xy2pix_nest, nside, ixp, iym, face_num, n_7
            xy2pix_nest, nside, ixp, iy , face_num, n_6
        end
        4:begin                 ;SouthEast edge
            other_face=8+ibp
            xy2pix_nest, nside, ixm, iy , face_num, n_2
            xy2pix_nest, nside, ixm, iyp, face_num, n_3
            xy2pix_nest, nside, ix , iyp, face_num, n_4
            xy2pix_nest, nside, ixp, iyp, face_num, n_5
            xy2pix_nest, nside, ixp, iy , face_num, n_6
            ipo=(swapLSBMSB(ipf, l64=l64) mod nsidesq) ;E-W flip
            pix2xy_nest, nside,ipo,ixo,iyo
            xy2pix_nest, nside, ixo, iyo+1, other_face, n_7
            n_8=other_face*nsidesq+ipo
            xy2pix_nest, nside, ixo, iyo-1, other_face, n_1
        end
        5:begin                 ;West corner
            nneigh=7
            other_face=8+ibm
            n_2=other_face*nsidesq+local_magic1
            n_1=n_2-1
            other_face=4+ib
            n_3=other_face*nsidesq
            n_4=n_3+1
            n_5=ipix0+1
            n_6=ipix0-1
            n_7=ipix0-2
        end
        6:begin                 ;North corner
            n_1=ipix0-3
            n_2=ipix0-1
            other_face=4+ib
            n_4=other_face*nsidesq+local_magic1
            n_3=n_4-1
            other_face=0+ib
            n_5=other_face*nsidesq
            other_face=4+ibp
            n_6=other_face*nsidesq+local_magic2
            n_7=n_6-2
            n_8=ipix0-2
        end
        7:begin                 ;South corner
            other_face=8+ib2
            n_1=other_face*nsidesq
            other_face=8+ibm
            n_2=other_face*nsidesq
            n_3=n_2+1
            n_4=ipix0+2
            n_5=ipix0+3
            n_6=ipix0+1
            other_face=8+ibp
            n_8=other_face*nsidesq
            n_7=n_8+2
        end
        8:begin                 ;East corner
            nneigh=7
            other_face=8+ibp
            n_7=other_face*nsidesq+local_magic2
            n_1=n_7-2
            n_2=ipix0-1
            n_3=ipix0+1
            n_4=ipix0+2
            other_face=4+ibp
            n_6=other_face*nsidesq
            n_5=n_6+2
        end
    endcase
endif                           ; south

list = (nneigh eq 8) ? [n_1, n_2, n_3, n_4, n_5, n_6, n_7, n_8] : [n_1, n_2, n_3, n_4, n_5, n_6, n_7]
; if (nneigh eq 8) then begin
;     list = [n_1, n_2, n_3, n_4, n_5, n_6, n_7, n_8]
; endif else begin
;     list = [n_1, n_2, n_3, n_4, n_5, n_6, n_7, n_8]
; endelse

return
end
