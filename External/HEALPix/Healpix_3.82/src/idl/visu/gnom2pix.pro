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
pro gnom2pix, xpos, ypos, id_pix, lon_deg, lat_deg, value
;+
; convert a uv position on a Gnomic map into a pixel number and
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

if (xpos ge !x.crange(0) and xpos le !x.crange(1) and ypos ge !y.crange(0) and ypos le !y.crange(1)) then begin

    vector = [1., flipconv*xpos, ypos] ; minus sign = astro convention
    ; ---------
    ; rotation
    ; ---------
    if (do_rot) then vector = vector # eul_mat
    vec2ang, vector, lat_deg, lon_deg, /astro
    if (do_conv) then vector = SKYCONV(vector, inco = coord_out, outco =  coord_in)
          ; we go from the final Gnomonic map (system coord_out) to
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
            vector = vector / (sqrt(total(vector^2, 2))#replicate(1,3)) ; normalize vector
            ; compute rotation of local coordinates around each vector
            tmp_sin = north_pole[1] * vector[*,0] - north_pole[0] * vector[*,1]
            tmp_cos = north_pole[2] - vector[*,2] * (north_pole[0:2] ## vector)
            if (flipconv lt 0) then tmp_cos = flipconv * tmp_cos
            phi = ATAN(tmp_sin, tmp_cos) ; angle in radians
        endif
        value = (data_plot[id_pix] - phi + 4*!PI) MOD (2*!PI) ; in 0,2pi
    endif else begin
;        if (n_elements(pixel_list) eq n_elements(data_plot)) then begin
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

