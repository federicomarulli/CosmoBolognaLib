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
pro orth2pix, xpos, ypos, id_pix, lon_deg, lat_deg, value
;+
; convert a uv position on a Orthographic map into a pixel number and
; (long, lat)
; only for scalar input
;
;-

@viewcom ; define common
do_poldirection = (polar_type eq 2)
indlist = (n_elements(pixel_list) eq n_elements(data_plot))

if (do_flip) then flipconv = +1 else flipconv = -1 ; longitude increase leftward by default (astro convention)

if (do_poldirection) then begin
    ; compute new position of pixelisation North Pole in the plot coordinates
    north_pole = [0.,0.,1.]
    if (do_conv) then north_pole = SKYCONV(north_pole, inco= coord_in, outco=coord_out)
    if (do_rot) then north_pole = north_pole # transpose(eul_mat)
endif

id_pix  = -1
lon_deg = -1000.
lat_deg = -1000.


if (do_fullsky) then begin
    c0 =  1
    c1 = -1
    valid = ((xpos-c0)^2 + ypos^2 le 1.) or ((xpos-c1)^2 + ypos^2 le 1.)
endif else begin
    valid = (xpos)^2 + ypos^2 le 1.
endelse

if (valid) then begin

    nested = keyword_set(nest)
    v1 = ypos
    u1 = flipconv * xpos
    if (do_fullsky) then begin
        sign = (u1 ge 0.)*2 - 1 ; =1 for u1>0 , =-1 otherwise
        ys = abs(u1)-c0
        xs = sqrt(1.-ys*ys-v1*v1)
        vector = [[sign*xs],[ys],[v1]] ; normalized vector
    endif else begin
        xs = sqrt(1.-u1*u1-v1*v1)
        vector = [[xs],[u1],[v1]] ; normalized vector
    endelse
        ; ---------
        ; rotation
        ; ---------
    if (do_rot) then vector = vector # eul_mat
    vec2ang, vector, lat_deg, lon_deg, /astro
    if (do_conv) then vector = SKYCONV(vector, inco = coord_out, outco =  coord_in)
                                ; we go from the final Orthographic map (system coord_out) to
                                ; the original one (system coord_in)
        ; -------------------------------------------------------------
        ; converts the position on the sphere into pixel number
        ; and project the corresponding data value on the map
        ; -------------------------------------------------------------
    case pix_type of
        'R' : begin
            vec2pix_ring, pix_param, vector, id_pix ; Healpix ring
            pix2vec_ring, pix_param, id_pix, vectorp ; vector pointing to Healpix pixel center
        end
        'N' : begin
            vec2pix_nest, pix_param, vector, id_pix ; Healpix nest
            pix2vec_nest, pix_param, id_pix, vectorp ; vector pointing to Healpix pixel center
        end
        'Q' : begin
            id_pix = uv2pix(vector, pix_param) ; QuadCube (COBE cgis software)
            vectorp = pix2uv(id_pix, pix_param) ; vector pointing to QuadCube pixel center
        end
        else : print,'error on pix_type'
    endcase
        ; --------------------------------
        ; deal with polarisation direction
        ; --------------------------------
    if (do_poldirection) then begin
        phi = 0.
        if (do_rot or do_conv) then begin
                                ; rotate pixel center to final coordinate system
            if (do_conv) then vectorp = SKYCONV(vectorp, inco= coord_in, outco=coord_out)
            if (do_rot) then vectorp = vectorp # transpose(eul_mat)
                                ; compute rotation of local coordinates around each vector
            tmp_sin = north_pole[1] * vectorp[*,0] - north_pole[0] * vectorp[*,1]
            tmp_cos = north_pole[2] - vectorp[*,2] * (north_pole[0:2] ## vectorp)
            if (flipconv lt 0) then tmp_cos = flipconv * tmp_cos
            phi = ATAN(tmp_sin, tmp_cos) ; angle in radians
        endif
        value = (data_plot[id_pix] - phi + 4*!PI) MOD (2*!PI) ; in 0,2pi
    endif else begin
        if (indlist) then begin
            value = sample_sparse_array(data_plot, id_pix, in=pixel_list, default=!healpix.bad_value)
        endif else begin
            value = data_plot[id_pix]
        endelse
    endelse
    if (n_elements(xpos) eq 1) then begin
        id_pix = id_pix[0]
        value = value[0]
    endif
endif

return
end
