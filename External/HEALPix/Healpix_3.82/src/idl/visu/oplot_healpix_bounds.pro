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
pro zphi2vec, z, phi, vec

x = cos(phi)
y = sin(phi)
sz = sqrt((1. - z)*(1. + z))

vec = [ [sz*x], [sz*y], [z] ] 

return
end



pro oplot_healpix_bounds, nside, eul_mat, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, orthographic=orthographic, azeq=azeq, flip = flip, _extra = oplot_kw, half_sky=half_sky, coordsys=coordsys
;+
; NAME:
;       OPLOT_HEALPIX_BOUNDS
;
; PURPOSE:
;       overplots Healpix pixels boundary on top
;       of existing map
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;       oplot_healpix_bounds, nside
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;      calls oplot_sphere
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;         Nov 2004
;   Dec 2010, EH: slight edition to avoid "Floating illegal operand" warning
;                 increased number of points (np) from 200 to 300
;-

identify_projection, projtype, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, azeq=azeq

if keyword_set(flip) then flipconv=1 else flipconv = -1  ; longitude increase leftward by default (astro convention)

; define variables
do_rot = (n_elements(eul_mat) eq 9)
do_ig = (n_elements(coordsys) eq 2) 
if (projtype eq 4) then begin  ; orthographic : deal with boundaries
    if keyword_set(half_sky) then begin
        nd = 1 ; number of half-sky disc
        c0 = 0
    endif else begin
        nd = 2
        c0 = 1
    endelse
endif

np = 300 > (2*nside)
ramp = dindgen(np)/(np-1.d0)
halfpi = !dpi/2.d0
fn  = double(nside)
fn2 = fn*fn
const = halfpi / fn / sqrt(3.d0) 
zero = replicate(0.d0, np)


for regime=0,1 do begin
    for kk=0L,nside-1 do begin
        if (regime eq 0) then begin ; equatorial regime
            z = 2.d0/3.d0 - ramp * (4.d0/3.d0)
            phi = halfpi*(ramp - kk/fn)
        endif else begin            ; polar regime
            ;z = (2.d0/3.d0 + ramp * (1.d0-kk^2/fn2)/3.d0) 
            r2 = (kk le (nside/2 +1)) ? ramp * (2.d0-ramp) : ramp; denser sampling in z close to poles, just like Healpix
            z = (2.d0/3.d0 + r2 * (1.d0-kk^2/fn2)/3.d0)
            phi = (kk gt 0) ? kk * const / sqrt(1.d0-z) : zero ; avoid 0/0 at pole
        endelse
        for quad=0,3 do begin ; quadrant
            for j=0,regime do begin
                sj = (-1)^j
                for i=0,1 do begin
                    si = (-1)^i
                    zphi2vec, si*z, (j+quad) *halfpi + sj*phi, vv
                    if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
                    if (do_rot) then vv = vv # transpose(eul_mat)
                    
                    case projtype of
                        1 : begin ; mollweide
                            vec2moll, vv, u, v
                            oplot_sphere, -flipconv * u, v, _extra = oplot_kw
                        end
                        2 : begin ; gnomic
                            keep = where(vv(*,0) gt 0, nk)
                            if (nk gt 0) then begin
                                u = vv(keep,1)/vv(keep,0) ; Y/X
                                v = vv(keep,2)/vv(keep,0) ; Z/X
                                oplot, flipconv * u, v, _extra = oplot_kw
                            endif
                        end
                        3 : begin ; cartesian
                            phic = atan(vv[*,1],vv[*,0]) ; longitude in [0,2pi]
                            theta = asin(vv[*,2]) ; latitude in [-pi/2,pi/2]
                            oplot_sphere, flipconv * phic, theta, _extra = oplot_kw
                        end
                        4 : begin ; orthographic
                            for sign = 1,1-nd,-2 do begin ; either (1,-1) or (1)
                                keep = where(vv[*,0]*sign ge 0, nk)
                                if (nk gt 0) then begin
                                    u = vv[keep,1] ; Y
                                    v = vv[keep,2] ; Z
                                    oplot_sphere, flipconv*(u+c0)*sign, v, _extra = oplot_kw
                                endif ; nk>0
                            endfor ; loop on sign
                        end
                        6: begin ; azimuthal equidistant
                            phia = atan(vv[*,2],vv[*,1]) ; in [0,2pi]
                            theta = acos(vv[*,0]) ; colatitude in [0,pi]
                            oplot_sphere, flipconv * theta * cos(phia), theta *sin(phia), _extra = oplot_kw
                        end
                        else: begin
                            message,projtype,' unsupported projection'
                        end
                    endcase
                endfor
            endfor
        endfor
    endfor
endfor


return
end
