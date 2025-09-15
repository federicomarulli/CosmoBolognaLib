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
pro neighbours_ring, nside, ipix, list, nneigh
;+
; NAME: 
;          neighbours_ring
;
; PURPOSE: 
;          find nearest neighbors pixels in Ring scheme
;
; CATEGORY:
;
; CALLING SEQUENCE: 
;          neighbours_ring, nside, ipix, list [, nneigh]
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
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;        calls NEIGHBOURS_NEST and converts to ring scheme
;        using RING2NEST and NEST2RING (and NSIDE2NPIX)
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;         June 2006: v1.0 relies on neighbours_nest, E.H.@IAP
;-

syntax = 'syntax: NEIGHBOURS_RING, Nside, Ipix, List [, Nneigh]'

if n_params() lt 3 then begin
    print,syntax
    list = [-1]
    nneigh = -1
    return
endif

npix = nside2npix(nside, error = error)
if (error ne 0) then begin
    message, /info, syntax
    message, 'Invalid Nside: '+string(nside)
endif

if (ipix ge npix or ipix lt 0 or n_elements(ipix) gt 1) then begin
    message,/info,syntax
    message,/info,'Ipix should be scalar in [0, ',npix-1,']'
endif

; turns central pixel into Nested indexing
ring2nest, nside, ipix, ipnest
; get list of neighbours in Nested scheme
neighbours_nest, nside, ipnest, listnest, nneigh
; get back to Ring scheme
nest2ring, nside, listnest, list


return
end
