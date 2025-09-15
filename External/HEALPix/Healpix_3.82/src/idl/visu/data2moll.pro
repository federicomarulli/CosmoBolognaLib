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
pro data2moll, data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
               planmap, Tmax, Tmin, color_bar, planvec, vector_scale, $
               PXSIZE=pxsize, LOG=log, HIST_EQUAL=hist_equal, MAX=max_set, MIN=min_set, FLIP=flip,$
               NO_DIPOLE=no_dipole, NO_MONOPOLE=no_monopole, UNITS = units, DATA_PLOT = data_plot, $
               GAL_CUT=gal_cut, POLARIZATION=polarization, SILENT=silent, PIXEL_LIST=pixel_list, ASINH=asinh, $
               TRUECOLORS=truecolors, DATA_TC=data_tc, MAP_OUT = map_out, ROT=rot_ang, FITS=fits
;+
;==============================================================================================
;     DATA2MOLL
;
;     turns a Healpix or Quad-cube map into in Mollweide egg
;
;     DATA2MOLL,  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in,
;                     coord_out, eul_mat
;          planmap, Tmax, Tmin, color_bar, planvec, vector_scale,
;          pxsize=, log=, hist_equal=, max=, min=, flip=, no_dipole=,
;          no_monopole=, units=, data_plot=, gal_cut=, polarization=, silent=,
;          pixel_list=, asinh=, truecolors=, data_tc=, map_out=, rot= , fits=
;
; IN :
;      data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in,
;      coord_out, eul_mat, rot
; OUT :
;      planmap, Tmax, Tmin, color_bar, planvec, vector_scale, map_out
; KEYWORDS
;      pxsize, log, hist_equal, max, min, flip, no_dipole, no_monopole, units,
;      polarization
;
;  called by mollview
;
;  HISTORY
; Sep 2007: added /silent
; April 2008: added pixel_list
; July 2008: added asinh
; May 2009: can deal with maps without any valid pixel
; April 2010: added Map_Out
;==============================================================================================
;-

do_true = keyword_set(truecolors)
truetype = do_true ? truecolors : 0
du_dv = 2.    ; aspect ratio
fudge = 1.02  ; spare some space around the Mollweide egg
if keyword_set(flip) then flipconv=1 else flipconv = -1  ; longitude increase leftward by default (astro convention)
if undefined(polarization) then polarization=0
do_polamplitude = (polarization[0] eq 1)
do_poldirection = (polarization[0] eq 2)
do_polvector    = (polarization[0] eq 3)
do_map_out      = arg_present(map_out)
do_fits         = keyword_set(fits)


!P.BACKGROUND = 1               ; white background
!P.COLOR = 0                    ; black foreground

mode_col = keyword_set(hist_equal)
mode_col = mode_col + 2*keyword_set(log) + 4*keyword_value(asinh, default=0, min=0, max=2)

sz = size(data)
obs_npix = sz[1]
bad_data= !healpix.bad_value

if (do_poldirection or do_polvector) then begin
    ; compute new position of pixelisation North Pole in the plot coordinates
    north_pole = [0.,0.,1.]
    if (do_conv) then north_pole = SKYCONV(north_pole, inco= coord_in, outco=coord_out)
    if (do_rot) then north_pole = north_pole # transpose(eul_mat)
endif
;-----------------------------------
; mask out some data
;-----------------------------------
;--------------------------------------
; remove monopole and/or dipole (for temperature, not polarisation)
;----------------------------------------
if not (do_poldirection or do_polamplitude) then begin
    if undefined(gal_cut) then bcut = 0. else bcut = abs(gal_cut)
    if keyword_set(no_dipole) then remove_dipole, data, $
      nside=pix_param, ordering=pix_type, units=units, coord_in = coord_in, coord_out=coord_out, $
      bad_data=bad_data, gal_cut=bcut, pixel=pixel_list
    if keyword_set(no_monopole) then remove_dipole, data, $
      nside=pix_param, ordering=pix_type, units=units, coord_in = coord_in, coord_out=coord_out, $
      bad_data=bad_data, gal_cut=bcut,/only, pixel=pixel_list
endif
; -------------------------------------------------------------
; create the rectangular window
; -------------------------------------------------------------
if DEFINED(pxsize) then xsize= LONG(pxsize>200) else xsize = 800L
ysize = xsize/2L
zsize = (do_true) ? 3 : 1
n_uv = xsize*ysize
indlist = (n_elements(pixel_list) eq obs_npix)
small_file = (n_uv GT obs_npix)  && ~do_map_out &&~do_fits
;small_file = ((n_uv GT npix)  and not do_poldirection)
dtype = size(data,/type) eq 5 ? 5 : 4 ; double (5) or float (4) by default

if (small_file) then begin
    ; file smaller than final map, make costly operation on the file
    ; initial data is destroyed and replaced by color
    if (do_poldirection or do_polvector) then begin
        phi = 0.
        if (do_rot or do_conv) then begin
            ; position of each map pixel after rotation and coordinate changes
            id_pix = lindgen(obs_npix)
            case pix_type of
                'R' : PIX2VEC_RING, pix_param, id_pix, vector ; Healpix ring
                'N' : PIX2VEC_NEST, pix_param, id_pix, vector; Healpix nest
                'Q' : vector = PIX2UV(pix_param, id_pix) ; QuadCube (COBE cgis software)
                else : print,'error on pix_type'
            endcase
            id_pix = 0
            if (do_conv) then vector = SKYCONV(vector, inco= coord_in, outco=coord_out)
            if (do_rot) then vector = vector # transpose(eul_mat)
            ; compute rotation of local coordinates around each vector
            tmp_sin = north_pole[1] * vector[*,0] - north_pole[0] * vector[*,1]
            tmp_cos = north_pole[2] - vector[*,2] * (north_pole[0:2] ## vector)
            if (flipconv lt 0) then tmp_cos = flipconv * tmp_cos
            phi = ATAN(tmp_sin, tmp_cos) ; angle in radians
            tmp_sin = 0. & tmp_cos = 0 & vector = 0.
        endif
        data_plot = data
        if (do_poldirection) then begin
            data = (data - phi + 4*!PI) MOD (2*!PI) ; angle
            min_set = 0. & max_set = 2*!pi
        endif
        if (do_polvector) then begin
            pol_data[*,1] = (pol_data[*,1] - phi + 4*!PI) MOD (2*!PI) ; angle is rotated
        endif
    endif else begin ; temperature only or polarisation amplitude only
        data_plot = data
    endelse
    ; color observed pixels
    if (do_true) then begin
        if (truetype eq 2) then begin
            for i=0,2 do begin
                find_min_max_valid, data_tc[*,i], mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
                data_tc[0,i] = COLOR_MAP(data_tc[*,i], mindata, maxdata, Obs, $
                                    color_bar = color_bar, mode=mode_col, silent=silent )
            endfor
        endif else begin
            find_min_max_valid, data_tc, mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
            data_tc = COLOR_MAP(data_tc, mindata, maxdata, Obs, $
                                color_bar = color_bar, mode=mode_col, $
                                minset = min_set, maxset = max_set, silent=silent )
        endelse
    endif else begin
        find_min_max_valid, data[*,0], mindata, maxdata, valid=Obs, bad_data= 0.9 * bad_data
        data = COLOR_MAP(data, mindata, maxdata, Obs, $
                         color_bar = color_bar, mode=mode_col, $
                         minset = min_set, maxset = max_set, silent=silent )
    endelse
    if (do_polvector) then begin ; rescale polarisation vector in each valid pixel
        pol_data[0,0] = vector_map(pol_data[*,0], Obs, vector_scale = vector_scale)
    endif
    if defined(Obs) then Obs = 0
    Tmin = mindata & Tmax = maxdata
    planmap = MAKE_ARRAY(/BYTE, xsize, ysize, zsize, Value = !P.BACKGROUND) ; white
endif else begin ; large file
    planmap = MAKE_ARRAY(type=dtype, xsize, ysize, zsize, Value = bad_data) 
    plan_off = 0L
endelse
if do_polvector then planvec = MAKE_ARRAY(type=dtype,xsize,ysize, 2, Value = bad_data) 

; -------------------------------------------------
; make the projection
;  we split the projection to avoid dealing with to big an array
; -------------------------------------------------
if (~keyword_set(silent)) then print,'... making the projection ...'
; -------------------------------------------------
; generate the (u,v) position on the mollweide map
; -------------------------------------------------
xll= 0 & xur =  xsize-1
yll= 0 & yur =  ysize-1
xc = 0.5*(xll+xur) & dx = (xur - xc)
yc = 0.5*(yll+yur) & dy = (yur - yc)


yband = LONG(5.e5 / FLOAT(xsize))
for ystart = 0, ysize - 1, yband do begin 
    yend   = (ystart + yband - 1) < (ysize - 1)
    nband = yend - ystart + 1
    u = FINDGEN(xsize)     # REPLICATE(1,nband)
    v = REPLICATE(1,xsize) # (FINDGEN(nband) + ystart)
    u =  du_dv*(u - xc)/(dx/fudge)   ; in [-2,2]*fudge
    v =        (v - yc)/(dy/fudge)   ; in [-1,1] * fudge

    ; -------------------------------------------------------------
    ; for each point on the mollweide map 
    ; looks for the corresponding position vector on the sphere
    ; -------------------------------------------------------------
    ellipse  = WHERE( (u^2/4. + v^2) LE 1. , nellipse)
    if (~small_file) then begin
        off_ellipse = WHERE( (u^2/4. + v^2) GT 1. , noff_ell)
        if (noff_ell NE 0) then plan_off = [plan_off, ystart*xsize+off_ellipse]
    endif
    if (nellipse gt 0) then begin
        u1 =  u(ellipse)
        v1 =  v(ellipse)
        u = 0 & v = 0
        s1 =  SQRT( (1-v1)*(1+v1) )
        a1 =  ASIN(v1)

        z = 2./!PI * ( a1 + v1*s1)
        phi = (flipconv *!Pi/2.) * u1/s1 ; lon in [-pi,pi], the minus sign is here to fit astro convention
        sz = SQRT( (1. - z)*(1. + z) )
        vector = [[sz * COS(phi)], [sz * SIN(phi)], [z]]
        u1 = 0 & v1 = 0 & s1 = 0 & a1 = 0 & z = 0 & phi = 0 & sz = 0
        ; --------------------------------
        ; deal with polarisation direction
        ; --------------------------------
        if ((do_poldirection || do_polvector) && ~small_file) then begin
            phi = 0.
            if (do_rot or do_conv) then begin
                ; compute rotation of local coordinates around each vector
                tmp_sin = north_pole[1] * vector[*,0] - north_pole[0] * vector[*,1]
                tmp_cos = north_pole[2] - vector[*,2] * (north_pole[0:2] ## vector)
                if (flipconv lt 0) then tmp_cos = flipconv * tmp_cos
                phi = ATAN(tmp_sin, tmp_cos) ; angle in radians
                tmp_sin = 0. & tmp_cos = 0
            endif
        endif
        ; ---------
        ; rotation
        ; ---------
        if (do_rot) then vector = vector # eul_mat
        if (do_conv) then vector = SKYCONV(vector, inco = coord_out, outco =  coord_in)
                                ; we go from the final Mollweide map (system coord_out) to
                                ; the original one (system coord_in)
        ; -------------------------------------------------------------
        ; converts the position on the sphere into pixel number
        ; and project the corresponding data value on the map
        ; -------------------------------------------------------------
        case pix_type of
            'R' : VEC2PIX_RING, pix_param, vector, id_pix ; Healpix ring
            'N' : VEC2PIX_NEST, pix_param, vector, id_pix ; Healpix nest
            'Q' : id_pix = UV2PIX(vector, pix_param)    ; QuadCube (COBE cgis software)
            else : print,'error on pix_type'
        endcase
        if (small_file) then begin ; (data and data_pol are already rescaled and color coded)
            if (do_true) then begin
                for i=0,zsize-1 do planmap[(ystart*xsize+i*n_uv)+ellipse] = data_tc[id_pix,i]
            endif else begin
                if (~ (do_polvector || do_polamplitude || do_poldirection) ) then begin
                    planmap[ystart*xsize+ellipse] = sample_sparse_array(data,id_pix,in_pix=pixel_list,default=2B) ; temperature
                endif else begin
                    planmap[ystart*xsize+ellipse] = data[id_pix]
                endelse
                if (do_polvector) then begin
                    planvec[ystart*xsize+ellipse]         = pol_data[id_pix,0] ; amplitude
                    planvec[(ystart*xsize+n_uv)+ellipse]  = pol_data[id_pix,1] ; direction
                endif
            endelse
        endif else begin ; (large file : do the projection first)
            if (do_true) then begin
                for i=0,zsize-1 do planmap[(ystart*xsize+i*n_uv)+ellipse] = data_tc[id_pix,i]
            endif else begin
                if (do_poldirection) then begin
                    planmap[ystart*xsize+ellipse] = (data[id_pix] - phi + 4*!PI) MOD (2*!PI) ; in 0,2pi
                endif else if (do_polvector) then begin
                    planmap[ystart*xsize+ellipse]         = data[id_pix] ; temperature
                    planvec[ystart*xsize+ellipse]         = pol_data[id_pix,0] ; amplitude
                    planvec[(ystart*xsize+n_uv)+ellipse]  = (pol_data[id_pix,1] - phi + 4*!PI) MOD (2*!PI) ; angle
                endif else begin ; temperature only or amplitude only
                                ;planmap[ystart*xsize+ellipse]         = data[id_pix] ; temperature
                    planmap[ystart*xsize+ellipse]         = sample_sparse_array(data,id_pix,in_pix=pixel_list,default=!healpix.bad_value) ; temperature
                endelse
            endelse
        endelse
    endif
    ellipse = 0 & id_pix = 0
endfor


;-----------------------------------
; export in FITS and as an array the original mollweide map before alteration
;----------------------------------------------

; planmap -> IDL array
if (do_map_out) then map_out = proj2map_out(planmap, offmap=plan_off, bad_data=bad_data)

; planmap -> FITS file
if keyword_set(fits) then begin 
    reso_arcmin = 60.d0 * 360.d0/(xsize-1) * fudge
    reso_arcmin *=  sqrt(8.d0) / !dpi ; WCS convention, ellipse surface is 4Pi
    proj2fits, planmap, fits, $
               projection = 'MOLL', flip=flip, $
               rot = rot_ang, coord=coord_out, reso = reso_arcmin, unit = sunits, min=mindata, max = maxdata
endif



if (small_file) then begin
    data = 0 & pol_data = 0
endif else begin
; file larger than final map, make
; costly coloring operation on the Mollweide map
    data_plot = temporary(data)
    pol_data = 0
    if (do_poldirection) then begin
        min_set = 0.
        max_set = 2*!pi
    endif
    find_min_max_valid, planmap, mindata, maxdata, valid= Obs, bad_data = 0.9 * bad_data
    case truetype of
        2: begin
                                ; truecolors=2 map each field to its color independently
            color = bytarr(xsize, ysize, zsize)
            for i=0,zsize-1 do begin
                find_min_max_valid, planmap[*,*,i], mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
                color[0,0,i] = COLOR_MAP(planmap[*,*,i], mindata, maxdata, Obs, $
                                         color_bar = color_bar, mode=mode_col, silent=silent)
            endfor
            planmap = color
        end
        3: begin
            intensity = total(planmap,3)/3.
            find_min_max_valid, intensity, mindata, maxdata, valid= Obs, bad_data = 0.1 * bad_data
            bint = COLOR_MAP(intensity, mindata, maxdata, Obs, $
                                color_bar = color_bar, mode=mode_col, $
                                minset = min_set, maxset = max_set, silent=silent)
;             for i=0,2 do begin
;                 ioff = i*xsize*ysize
;                 color[Obs+ioff] = 3B + bytscl((planmap[Obs+ioff]>0)/intensity[Obs] * bint[Obs], min=0,top=252)
;             endfor
            ioff = xsize * ysize
            mat = planmap[[[Obs], [Obs+ioff], [Obs+2*ioff]]] > 0
            mat2 = (bint[Obs]/intensity[Obs]) # [1,1,1]
            color = 3B + bytscl(mat*mat2, min=0, top=252)
            planmap = MAKE_ARRAY(/BYTE, xsize, ysize, 3, Value = 2B)
            for i=0,2 do planmap[Obs + i*ioff] = color[*,i]
            
        end
        else: begin
                                ; same for truecolors=1 and false colors:
            planmap = COLOR_MAP(planmap, mindata, maxdata, Obs, $
                                color_bar = color_bar, mode=mode_col, $
                                minset = min_set, maxset = max_set, silent=silent)
        end
    endcase
    for i=0,zsize-1 do planmap[plan_off+i*n_uv] = !p.background ; white
    if (do_polvector) then begin ; rescale polarisation vector in each valid pixel
        planvec[*,*,0] = vector_map(planvec[*,*,0], Obs, vector_scale = vector_scale)
        planvec[plan_off] = -1
    endif
    Obs = 0 & plan_off = 0
    Tmin = mindata & Tmax = maxdata
endelse


return
end

