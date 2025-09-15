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
pro data2orth, data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
               planmap, Tmax, Tmin, color_bar, planvec, vector_scale, $
               PXSIZE=pxsize, LOG=log, HIST_EQUAL=hist_equal, MAX=max_set, MIN=min_set, FLIP=flip,$
               NO_DIPOLE=no_dipole, NO_MONOPOLE=no_monopole, UNITS = units, DATA_PLOT = data_plot, $
               GAL_CUT=gal_cut, POLARIZATION=polarization, HALF_SKY=half_sky, SILENT=silent, PIXEL_LIST=pixel_list, $
               ASINH=asinh, $
               DO_SHADE=shade, SHADEMAP = shademap, TRUECOLORS=truecolors, DATA_TC=data_tc, $
               MAP_OUT = map_out, ROT=rot_ang, FITS=fits, STAGGER=stagger
;+
;==============================================================================================
;     DATA2ORTH
;
;     turns a Healpix or Quad-cube map into in Orthographic disc(s)
;
;     DATA2ORTH,  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in,
;                     coord_out, eul_mat
;          planmap, Tmax, Tmin, color_bar, planvec, vector_scale,
;          pxsize=, log=, hist_equal=, max=, min=, flip=, no_dipole=,
;          no_monopole=, units=, data_plot=, gal_cut=,
;          polarization=, half_sky, silent=, pixel_list, asinh=,
;          do_shade=, shademap=, truecolors=, data_tc=, map_out=, rot=, fits=, stagger=
;
; IN :
;      data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in,
;      coord_out, eul_mat, rot
; OUT :
;      planmap, Tmax, Tmin, color_bar, planvec, vector_scale, shademap, map_out
; KEYWORDS
;      pxsize, log, hist_equal, max, min, flip, no_dipole, no_monopole, units,
;      polarization, half_sky, do_shade, fits
;
;  called by mollview
;  HISTORY
; Sep 2007: added /silent
; April 2008: added pixel_list
; July 2008: added asinh
; May 2009: added do_shade, shademap
;           can deal with map without any valid pixel
; May 2010: added Map_Out and FITS
; Dec 2011: STAGGER implemented
; Jan 2012: added a twist to STAGGER
;==============================================================================================
;-

do_shade = keyword_set(shade)
do_true = keyword_set(truecolors)
truetype = do_true ? truecolors : 0
do_stagger =  keyword_set(stagger)
do_fullsky = ~ (keyword_set(half_sky) || do_stagger)
nspheres = 1
if (do_fullsky) then nspheres = 2
if (do_stagger) then begin
    nspheres = 3
    if (stagger[0] le 0. || stagger[0] gt 2.) then begin
        message_patch, level=-1, 'Stagger must be in ]0,2]'
    endif
endif

if (nspheres gt 1) then begin
    du_dv = 2. 
    xsize_min = 200
    xsize_default = 800L
endif else begin
    du_dv = 1.                  ; aspect ratio
    xsize_min = 200
    xsize_default = 600L
endelse
fudge = 1.02  ; spare some space around the Orthographic disc(s)
if keyword_set(flip) then flipconv=1 else flipconv = -1  ; longitude increase leftward by default (astro convention)
if undefined(polarization) then polarization=0
do_polamplitude = (polarization[0] eq 1)
do_poldirection = (polarization[0] eq 2)
do_polvector    = (polarization[0] eq 3)
do_map_out      = arg_present(map_out)
do_fits         = keyword_set(fits)
;vec_shine = [1,1,1]/sqrt(3.)
vec_shine = [0.5,-1*flipconv,1] ; light source position (do_shade option)
vec_shine = vec_shine / sqrt(total(vec_shine^2))


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
      nside=pix_param, ordering=pix_type, units=units, coord_in = coord_in, coord_out=coord_out, bad_data=bad_data, gal_cut=bcut
    if keyword_set(no_monopole) then remove_dipole, data, $
      nside=pix_param, ordering=pix_type, units=units, coord_in = coord_in, coord_out=coord_out, bad_data=bad_data, gal_cut=bcut,/only
endif
; -------------------------------------------------------------
; create the rectangular window
; -------------------------------------------------------------
if DEFINED(pxsize) then xsize= LONG(pxsize>xsize_min) else xsize = xsize_default
ysize = xsize/du_dv
zsize = (do_true) ? 3 : 1
n_uv = xsize*ysize
indlist = (n_elements(pixel_list) eq obs_npix)
small_file = (n_uv GT obs_npix)   && ~do_map_out &&~do_fits
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
        find_min_max_valid, data, mindata, maxdata, valid=Obs, bad_data= 0.9*bad_data
        data    = color_map(data, mindata, maxdata, Obs, $
                         color_bar = color_bar, mode=mode_col, $
                         minset = min_set, maxset = max_set, silent=silent)
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
if do_polvector then planvec  = MAKE_ARRAY(type=dtype,xsize,ysize, 2, Value = bad_data) 
if do_shade     then shademap = MAKE_ARRAY(type=dtype,xsize,ysize,    Value = 1.0) 

; -------------------------------------------------
; make the projection
;  we split the projection to avoid dealing with too big an array
; -------------------------------------------------
if (~keyword_set(silent)) then print,'... making the projection ...'
; -------------------------------------------------
; generate the (u,v) position on the orthographic map
; -------------------------------------------------
xll= 0 & xur =  xsize-1
yll= 0 & yur =  ysize-1
xc = 0.5*(xll+xur) & dx = (xur - xc)
yc = 0.5*(yll+yur) & dy = (yur - yc)

if (nspheres eq 2) then begin
    c0 =  1
    c1 = -1
endif
if (nspheres eq 3) then begin
    cstag =   abs(stagger[0])*[1,0,-1]
    ; astag = (n_elements(stagger) gt 1) ? double(stagger[1]) : 0.0_dp
    do_twist = (n_elements(stagger) gt 1)
    if (do_twist) then begin
        astag = stagger[1] * !DtoR ; twist angle converted to Radians
        rotstag = [[cos(astag), 0, sin(astag)],[0,1,0],[-sin(astag),0,cos(astag)]]
    endif
endif

yband = LONG(5.e5 / FLOAT(xsize))
for ystart = 0, ysize - 1, yband do begin 
    yend   = (ystart + yband - 1) < (ysize - 1)
    nband = yend - ystart + 1
    u = FINDGEN(xsize)     # REPLICATE(1,nband)
    v = REPLICATE(1,xsize) # (FINDGEN(nband) + ystart)
    u =  flipconv*du_dv*(u - xc)/(dx/fudge)   ; in [-2,2]*fudge or [-1,1]*fudge
    v =                 (v - yc)/(dy/fudge)   ; in [-1,1] * fudge
    ; -------------------------------------------------------------
    ; for each point on the orthographic map 
    ; looks for the corresponding position vector on the sphere
    ; -------------------------------------------------------------
    case nspheres of
        3: disc  = WHERE( (((u-cstag[0])^2 + v^2) LE 1.) OR $
                          (((u-cstag[1])^2 + v^2) LE 1.) OR $
                          (((u-cstag[2])^2 + v^2) LE 1.), ndisc, $
                       complement=off_disc, ncomplement=noff_ell)
        2: disc  = WHERE( (((u-c0)^2 + v^2) LE 1.) OR (((u-c1)^2 + v^2) LE 1.), ndisc, $
                       complement=off_disc, ncomplement=noff_ell)
        1: disc  = WHERE( ((u^2 + v^2) LE 1.), ndisc, $
                       complement=off_disc, ncomplement=noff_ell)
    endcase
    if (~small_file && noff_ell NE 0) then plan_off = [plan_off, ystart*xsize+off_disc]
    if (ndisc gt 0) then begin
        u1 =  u[disc]
        v1 =  v[disc]
        u = 0 & v = 0
        case nspheres of
            3: begin
                trigger = (u1/cstag[0]) < 1
                trigger >= (-1)
                trigger = nint(trigger) ; -1, 0 or 1
                ys = u1 - trigger*cstag[0]
            end
            2: begin
                sign = (u1 ge 0.)*2 - 1 ; =1 for u1>0 , =-1 otherwise
                ys = abs(u1)-c0 
            end
            1: begin
                ys = u1
            end
        endcase
        xs = sqrt(1.-ys*ys-v1*v1)
        if (nspheres eq 2) then xs *= sign
        vector = [[xs],[ys],[v1]] ; normalized vector
        if (do_shade) then begin
            ; sphere shading (Phong model)
            ; ambiant light = 50%, specular light = 90%, diffuse light = ( 1 - ambiant)
            ambiant = 0.5 ;0.5
            spec_frac = 0.9
            ; specular: (r.v) = (2 (n.i) n - i). v
            ; i = incident light direction
            ; r = reflected light direction, |r| = 1
            ; n = normal to surface
            ; v = observer (here at (+/-Inf, 0, 0) )
            cos_in = vec_shine ## vector
            spec = (2 * cos_in * vector[*,0] - vec_shine[0])
            spec = spec * ( cos_in gt 0) ; no specular highlight accros the sphere
            shade = (ambiant + (1.0 - ambiant) * cos_in + spec_frac* abs(spec)^16) > 0.
            spec = 0 & cos_in = 0.
        endif
        if (nspheres eq 3 && do_twist) then begin ; extra rotation for stagger with twist
            km = where(trigger lt 0, nkm)
            kp = where(trigger gt 0, nkp)
            if nkm gt 0 then vector[km,0:2] = vector[km,0:2] # rotstag
            if nkp gt 0 then vector[kp,0:2] = vector[kp,0:2] # transpose(rotstag)
        endif

        u1 = 0 & v1 = 0 & sign = 0 & xs = 0 & ys = 0
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
                                ; we go from the final Orthographic map (system coord_out) to
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
        if (small_file) then begin ; (data, data_pol and data_tc are already rescaled and color coded)
            if (do_true) then begin
                for i=0,zsize-1 do planmap[(ystart*xsize+i*n_uv)+disc] = data_tc[id_pix,i]
            endif else begin
                if (~ (do_polvector || do_polamplitude || do_poldirection) ) then begin
                    planmap[ystart*xsize+disc] = sample_sparse_array(data,id_pix,in_pix=pixel_list,default=2B) ; temperature
                endif else begin
                    planmap[ystart*xsize+disc] = data[id_pix]
                endelse
                if (do_polvector) then begin
                    planvec[ystart*xsize+disc]         = pol_data[id_pix,0] ; amplitude
                    planvec[(ystart*xsize+n_uv)+disc]  = pol_data[id_pix,1] ; direction
                endif
            endelse
        endif else begin ; (large file : do the projection first)
            if (do_true) then begin
                for i=0,zsize-1 do planmap[(ystart*xsize+i*n_uv)+disc] = data_tc[id_pix,i]
            endif else begin
                if (do_poldirection) then begin
                    planmap[ystart*xsize+disc] = (data[id_pix] - phi + 4*!PI) MOD (2*!PI) ; in 0,2pi
                endif else if (do_polvector) then begin
                    planmap[ystart*xsize+disc]         = data[id_pix] ; temperature
                    planvec[ystart*xsize+disc]         = pol_data[id_pix,0] ; amplitude
                    planvec[(ystart*xsize+n_uv)+disc]  = (pol_data[id_pix,1] - phi + 4*!PI) MOD (2*!PI) ; angle
                endif else begin ; temperature only or amplitude only
                                ;planmap[ystart*xsize+disc]         = data[id_pix] ; temperature
                    planmap[ystart*xsize+disc]      = sample_sparse_array(data,id_pix,in_pix=pixel_list,default=!healpix.bad_value) ; temperature
                endelse
            endelse
        endelse
        if (do_shade) then shademap[ystart*xsize+disc] = shade
    endif
    disc = 0 & id_pix = 0 & shade = 0
endfor

;-----------------------------------
; export in FITS and as an array the original mollweide map before alteration
;----------------------------------------------

; planmap -> IDL array (implies small_file)
if (do_map_out) then map_out = proj2map_out(planmap, offmap=off_disc, bad_data=bad_data)

; planmap -> FITS file (implies small_file)
if keyword_set(fits) then begin 
    reso_arcmin = 60.d0 * 180.d0/xsize * fudge * 2.d0 / !dpi
    ;;print,reso_arcmin,xsize
    proj2fits, planmap, fits, $
               projection = 'ORTH', flip=flip, $
               rot = rot_ang, coord=coord_out, reso = reso_arcmin, unit = sunits, min=mindata, max = maxdata, $
               half_sky = half_sky
endif

if (small_file) then begin
    data = 0 & pol_data = 0 & data_tc = 0
endif else begin
; file larger than final map, make
; costly coloring operation on the Orthographic map
    data_plot = temporary(data)
    pol_data = 0
    find_min_max_valid, planmap, mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
    if (do_poldirection) then begin
        min_set = 0.
        max_set = 2*!pi
    endif
    if (truetype eq 2) then begin
        ; truecolors=2 map each field to its color independently
        color = bytarr(xsize,ysize,zsize)
        for i=0,zsize-1 do begin
            find_min_max_valid, planmap[*,*,i], mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
            color[0,0,i] = COLOR_MAP(planmap[*,*,i], mindata, maxdata, Obs, $
                          color_bar = color_bar, mode=mode_col, silent=silent)
        endfor
        planmap = color
    endif else begin
        ; same for truecolors=1 and false colors:
        planmap = COLOR_MAP(planmap, mindata, maxdata, Obs, $
                            color_bar = color_bar, mode=mode_col, $
                            minset = min_set, maxset = max_set, silent=silent)
    endelse
    for i=0,zsize-1 do planmap[plan_off+i*n_uv] = !p.background
    if (do_polvector) then begin ; rescale polarisation vector in each valid pixel
        planvec[*,*,0] = vector_map(planvec[*,*,0], Obs, vector_scale = vector_scale)
        planvec[plan_off] = -1
    endif
    Obs = 0 & plan_off = 0
    Tmin = mindata & Tmax = maxdata
endelse

return
end

