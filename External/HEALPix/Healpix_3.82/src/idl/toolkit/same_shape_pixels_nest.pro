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
pro same_shape_pixels_nest, nside, template, list, reflexion, nreplications=nrep
;+
; NAME:
;         same_shape_pixels_nest
;
; PURPOSE:
;         Returns the list of NESTED ordered pixels having the same shape as the
;         template provided
;
; CATEGORY:
;
; CALLING SEQUENCE:
;         SAME_SHAPE_PIXELS_NEST, Nside, Template, ListPix, [Reflexion, Nreplications=]
;
; INPUTS:
;         Nside    : scalar, integer, Healpix resolution parameter
;         Template : scalar, integer, index of the pixel template
;           should be in [0, (1 + Nside*(Nside+6))/4 - 1]
;           (see Nside2Ntemplates)
;
; KEYWORD PARAMETERS:
;         Nreplications: on output, number of pixels having the same shape as template
;          (depending on the template, Nreplications is either
;          8, 16, 4*Nside or 8*Nside)
;
; OUTPUTS:
;         ListPix : vector, integer, list of pixel (NEST numbered) having the
;          same shape as Template, see Notes below
;
;         Reflexion : vector, integer in [0,3], see Notes below.
;            (same size as ListPix)
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;         None
;
; SIDE EFFECTS:
;         none
;
; RESTRICTIONS:
;         ??
;
; PROCEDURE:
;        calls SAME_SHAPE_PIXELS_RING, RING2NEST and SORT
;
; EXAMPLE:
;
;
; RELATED ROUTINES:
;         NEST2RING, RING2NEST
;         SAME_SHAPE_PIXELS_NEST, SAME_SHAPE_PIXELS_RING
;         TEMPLATE_PIXEL_NEST,    TEMPLATE_PIXEL_RING
;         NSIDE2NTEMPLATES
;
; MODIFICATION HISTORY:
;         2004-09, EH, v1.0
;
;
; Notes: 
;         see Notes in TEMPLATE_PIXEL_RING
;
;
;-

routine = 'same_shape_pixels_nest'
syntax  = 'syntax: '+routine+', Nside, Template, ListPix, [Reflexion, Nreplications=]'
nparams = n_params()
if nparams lt 3 or nparams gt 4 then begin
    print,syntax
    message,'Abort'
    return
endif


if (N_ELEMENTS(nside) GT 1) then message,'Nside should be a scalar in '+routine
npix = nside2npix(nside, error = error)

if (error ne 0) then message,'Invalid Nside '+string(nside)
lnside = long(nside[0]) ; 2018-06-06

if (N_ELEMENTS(template) GT 1) then message,'template should be a scalar in '+routine
ntplt= (lnside*(lnside+6))/4 > 2

if (template lt 0 or template ge ntplt) then begin
    message,/info,'Error on template argument'
    print,syntax
    print,'Nside = ',nside,', Template = ',template
    print,'Template should be in [0, Max(Nside*(Nside+6)/4,2)-1=',ntplt-1,']'
    message,'Abort'
endif

same_shape_pixels_ring, nside, template, list_ring, reflexion_ring, nreplications=nrep
ring2nest, nside, list_ring, list

reshuf    = sort(list)
list      = list[reshuf]
reflexion = reflexion_ring[reshuf]

return
end
