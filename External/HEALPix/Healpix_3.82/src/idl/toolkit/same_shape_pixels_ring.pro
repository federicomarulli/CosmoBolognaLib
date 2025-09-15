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
pro same_shape_pixels_ring, nside, template, list, reflexion, nreplications=nrep
;+
; NAME:
;         same_shape_pixels_ring
;
; PURPOSE:
;         Returns the list of RING ordered pixels having the same shape as the
;         template provided
;
; CATEGORY:
;
; CALLING SEQUENCE:
;         SAME_SHAPE_PIXELS_RING, Nside, Template, ListPix, [Reflexion, Nreplications=]
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
;         ListPix : vector, integer, list of pixel (RING numbered) having the
;          same shape as Template, see Notes below
;
;         Reflexion : vector, integer in [0,3], see Notes below
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
;
; EXAMPLE:
;
;
; RELATED ROUTINES:
;         NEST2RING, RING2NEST
;         SAME_SHAPE_PIXELS_NEST, SAME_SHAPE_PIXELS_RING
;         TEMPLATE_PIXEL_NEST,    TEMPLATE_PIXEL_RING,
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

routine = 'same_shape_pixels_ring'
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

ncap = 2*lnside*(lnside+1)
ns4  = 4*lnside
n0   = (lnside*(lnside+2))/4 > 1

if (template ge n0) then begin
; Non Polar pixels (most frequent)
    ramp = lindgen(ns4)
    zero = replicate(0,ns4)
    
    in0 = ncap + (template - n0)*ns4 ; Number of pixels North of current template

    if (template lt ntplt-1) then begin
        ; Non- Equator
        is0 = npix - in0 - ns4
        nrep = ns4 * 2
        list = [in0+ramp, is0+ramp]
        reflexion = [zero, zero+2]
    endif else begin
        ; Equator
        nrep = ns4
        list = [in0+ramp]
        reflexion = zero
    endelse

endif else begin
    iring = long(sqrt(4*template+1)) ; in [1,Nside]
    ifi   = template - (iring*iring)/4L
    in0   = 2L*iring*(iring-1L) ; Number of pixels North of current template
    is0   = npix - in0 - 4*iring ; Number of pixels North of south image of current template

    ramp4 = [0,1,2,3]
    if (2L*ifi+1 eq iring) then begin
        ; pixel on phi = Pi/4
        nrep = 8
        list = [in0+ramp4*iring, is0+ramp4*iring] + ifi
        zero4 = replicate(0,4)
        reflexion = [zero4, zero4 + 2]
    endif else begin
        nrep = 16
        rr = [[ramp4*iring+ifi],[(ramp4+1)*iring-ifi-1]]
        rr = reform(transpose(rr),8) ; [ifi, iring-ifi-1,iring+ifi, 2*iring-ifi-1,2*iring+ifi,..]
        list = [in0+rr, is0+rr]
        alt8  = indgen(8) mod 2 ; [0,1,0,1,...]
        reflexion = [alt8, alt8+2]
    endelse
endelse

; map = replicate(-1,npix) & map[list]=reflexion
; mollview,/online,map,/flip ; ,min=-1,max=4


return
end
