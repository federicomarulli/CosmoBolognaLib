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
pro identify_projection, projtype, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, orthographic=orthographic, diamonds=diamonds, azeq=azeq
;+
; identify_projection, projtype, projection=, mollweide=, $
;                                gnomic=, cartesian=, orthographic=, diamonds=, azeq=
;
;
;
;  OUTPUT : projtype
;    1 : mollweide
;    2 : gnomic
;    3 : cartesian
;    4 : orthographic
;    5 : diamonds
;    6 : azeq (azimuthal equidistant)
;
;-
proj = ' '
if keyword_set(projection) then proj = strupcase(strmid(projection,0,4))
do_moll = (keyword_set(mollweide)    || proj eq 'MOLL')
do_gnom = (keyword_set(gnomic)       || proj eq 'GNOM')
do_cart = (keyword_set(cartesian)    || proj eq 'CART')
do_orth = (keyword_set(orthographic) || proj eq 'ORTH')
do_diam = (keyword_set(diamonds)     || proj eq 'DIAM')
do_azeq = (keyword_set(azeq)         || proj eq 'AZEQ' || proj eq 'AZIM')

if ((do_moll+do_gnom+do_cart+do_orth+do_diam+do_azeq) ne 1) then begin
    message,'should be either Mollweide or Gnomic or Cartesian or Diamonds or AzEq'
endif
; if ((do_moll+do_gnom+do_cart+do_orth) ne 1) then begin
;     message,'should be either Mollweide or Gnomic or Cartesian'
; endif

molltype = 1
gnomtype = 2
carttype = 3
orthtype = 4
diamtype = 5
azeqtype = 6

projtype = do_moll * molltype    $
         + do_gnom * gnomtype    $
         + do_cart * carttype    $
         + do_orth * orthtype    $
         + do_diam * diamtype    $
         + do_azeq * azeqtype

return
end
