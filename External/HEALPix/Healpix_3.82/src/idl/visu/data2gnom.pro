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
pro data2gnom, data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
               color, Tmax, Tmin, color_bar, dx, planvec, vector_scale, $
               PXSIZE=pxsize, PYSIZE=pysize, ROT=rot_ang, LOG=log, HIST_EQUAL=hist_equal, $
               MAX=max_set, MIN=min_set, $
               RESO_ARCMIN=reso_arcmin, FITS = fits, $
               FLIP=flip, DATA_plot = data_plot, $
               POLARIZATION=polarization, SILENT=silent, PIXEL_LIST=pixel_list, ASINH=asinh, $
               TRUECOLORS=truecolors, DATA_TC=data_tc, MAP_OUT=map_out

;+
;==============================================================================================
;     DATA2GNOM
;
;     turns a Healpix or Quad-cube map into in Gnomonic rectangular map
;
;     DATA2GNOM,  data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat,
;          color, Tmax, Tmin, color_bar, dx, planvec, vector_scale,
;          pxsize=, pysize=, rot=, log=, hist_equal=, max=, min=,
;          reso_arcmin=, fits=, flip=, data_plot=, polarization=, silent=,
;          pixel_list=, TRUECOLORS=, DATA_TC=, MAP_OUT=
;
; IN :
;      data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat
; OUT :
;      color, Tmax, Tmin, color_bar, dx, planvec, vector_scale
; KEYWORDS
;      Pxsize, Pysize, Rot, Log, Hist_equal, Max, Min, Reso_arcmin,
;      Fits, flip, data_plot, polarization, pixel_list, asinh, map_out
;
;  called by gnomview
;
;  HISTORY; Feb 2005: added small_file to avoid pol direction variation within pixels
; Sep 2007: added /silent
; April 2008: added pixel_list=
; July 2008: added asinh
; May 2009: can deal with maps without any valid pixel
;==============================================================================================
;-

do_true = keyword_set(truecolors) 
truetype = do_true ? truecolors : 0
proj_small = 'gnomic'
du_dv = 1.    ; aspect ratio
fudge = 1.00  ; 
if keyword_set(flip) then flipconv=1 else flipconv = -1  ; longitude increase leftward by default (astro convention)
if undefined(polarization) then polarization=0
do_polamplitude = (polarization[0] eq 1)
do_poldirection = (polarization[0] eq 2)
do_polvector    = (polarization[0] eq 3)

!P.BACKGROUND = 1               ; white background
!P.COLOR = 0                    ; black foreground

mode_col = keyword_set(hist_equal)
mode_col = mode_col + 2*keyword_set(log) + 4*keyword_value(asinh, default=0, min=0, max=2)

obs_npix = n_elements(data)
npix_full = (pix_type eq 'Q') ? 6*(4L)^(pix_param-1) : nside2npix(pix_param)
bad_data= !healpix.bad_value

if (do_poldirection or do_polvector) then begin
    ; compute new position of pixelisation North Pole in the plot coordinates
    north_pole = [0.,0.,1.]
    if (do_conv) then north_pole = SKYCONV(north_pole, inco= coord_in, outco=coord_out)
    if (do_rot) then north_pole = north_pole # transpose(eul_mat)
endif
; -------------------------------------------------------------
; create the rectangular window
; -------------------------------------------------------------
if defined(pxsize) then xsize = pxsize*1L else xsize = 500L
if defined(pysize) then ysize = pysize*1L else ysize = xsize
if defined(reso_arcmin) then resgrid = reso_arcmin/60. else resgrid = 1.5/60.
dx      = resgrid * !DtoR
zsize = (do_true) ? 3 : 1
N_uv = xsize*ysize
indlist = (n_elements(pixel_list) eq n_elements(data[*,0]))
small_file = ((!pi*4./dx^2 GT npix_full && do_poldirection))
dtype = size(data,/type) eq 5 ? 5 : 4 ; double (5) or float (4) by default


if (~keyword_set(silent)) then begin
    print,'Input map  :  ',3600.*6.d0/sqrt(!dpi*npix_full),' arcmin / pixel ',form='(a,f8.3,a)'
    print,'gnomonic map :',resgrid*60.,' arcmin / pixel ',xsize,'*',ysize,form='(a,f8.3,a,i4,a,i4)'
endif

if (small_file) then begin
    ; file smaller than final map, make costly operation on the file
    ; initial data is destroyed and replaced by color
    if (do_poldirection or do_polvector) then begin
        phi = 0.
        if (do_rot or do_conv) then begin
            ; position of each map pixel after rotation and coordinate changes
            if (indlist) then begin
                id_pix = pixel_list
            endif else begin
                id_pix = lindgen(npix_full)
            endelse
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
        find_min_max_valid, data, mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
        data    = COLOR_MAP(data, mindata, maxdata, Obs, $
                         color_bar = color_bar, mode=mode_col, $
                         minset = min_set, maxset = max_set, silent=silent )
    endelse
    if (do_polvector) then begin ; rescale polarisation vector in each valid pixel
        pol_data[0,0] = vector_map(pol_data[*,0], Obs, vector_scale = vector_scale)
    endif
    if defined(Obs) then Obs = 0
    Tmin = mindata & Tmax = maxdata
    color = MAKE_ARRAY(/BYTE, xsize, ysize, zsize, Value = !P.BACKGROUND) ; white
    grid  = MAKE_ARRAY(type=dtype, xsize, ysize)
endif else begin ; large
    grid  = MAKE_ARRAY(type=dtype, xsize, ysize, zsize)
endelse
if do_polvector then planvec = MAKE_ARRAY(type=dtype,xsize,ysize, 2, Value = bad_data) 
; -------------------------------------------------------------
; makes the projection around the chosen contact point
; -------------------------------------------------------------
; position on the planar grid  (1,u,v)
x0 = +1.
xll= 0 & xur =  xsize-1
yll= 0 & yur =  ysize-1
xc = 0.5*(xll+xur)  ; & deltax = (xur - xc)
yc = 0.5*(yll+yur)  ; & deltay = (yur - yc)

yband = LONG(5.e5 / FLOAT(xsize))
for ystart = 0, ysize - 1, yband do begin 
    yend   = (ystart + yband - 1) < (ysize - 1)
    nband = yend - ystart + 1
    npb = xsize * nband
    u = flipconv*(FINDGEN(xsize) - xc)# REPLICATE(dx,nband)   ; minus sign = astro convention
    v =           REPLICATE(dx,xsize) # (FINDGEN(nband) + ystart - yc)
    x = replicate(x0, npb)
    vector = [[x],[reform(u,npb,/over)],[reform(v,npb,/over)]] ; non normalised vector
    ; --------------------------------
    ; deal with polarisation direction
    ; --------------------------------
    if (do_poldirection or do_polvector) then begin
        phi = 0.
        if (do_rot or do_conv) then begin
            vector = vector / (sqrt(total(vector^2, 2))#replicate(1,3)) ; normalize vector
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
          ; we go from the final Gnomonic map (system coord_out) to
          ; the original one (system coord_in)
    ; -------------------------------------------------------------
    ; converts the position on the sphere into pixel number
    ; and project the corresponding data value on the map
    ; -------------------------------------------------------------
    case pix_type of
        'R' : VEC2PIX_RING, pix_param, vector, id_pix ; Healpix ring
        'N' : VEC2PIX_NEST, pix_param, vector, id_pix ; Healpix nest
        'Q' : id_pix = UV2PIX(vector, pix_param) ; QuadCube (COBE cgis software)
        else : print,'error on pix_type'
    endcase
    if (small_file) then begin ; (data and data_pol are already rescaled and color coded)
        if (do_true) then begin
            for i=0,zsize-1 do color[ystart*xsize+i*n_uv] = data[id_pix,i]
            grid[ystart*xsize]  = data_plot[id_pix] ; unaltered data
        endif else begin
            color[ystart*xsize] = data[id_pix]
            grid[ystart*xsize]  = data_plot[id_pix] ; unaltered data
            if (do_polvector) then begin
                planvec[ystart*xsize]       = pol_data[id_pix,0] ; amplitude
                planvec[ystart*xsize+n_uv]  = pol_data[id_pix,1] ; direction
            endif
        endelse
    endif else begin            ; (large file : do the projection first)
        if (do_true) then begin
            for i=0,zsize-1 do grid[ystart*xsize+i*n_uv] = data_tc[id_pix,i]
        endif else begin
            if (do_poldirection) then begin
                grid[ystart*xsize] = (data[id_pix] - phi + 4*!PI) MOD (2*!PI) ; in 0,2pi
            endif else if (do_polvector) then begin
                grid[ystart*xsize]         = data[id_pix]
                planvec[ystart*xsize]      = pol_data[id_pix,0]
                planvec[ystart*xsize+n_uv] = (pol_data[id_pix,1] - phi + 4*!PI) MOD (2*!PI) ; angle
            endif else begin
;;;            grid[ystart*xsize] = data[id_pix]
                grid[ystart*xsize] = sample_sparse_array(data, id_pix, in_pix=pixel_list, default= !healpix.bad_value)
            endelse
        endelse
    endelse
endfor
u = 0 & v = 0 & x = 0 & vector = 0

; -------------------------------------------------------------
; Test for unobserved pixels
; -------------------------------------------------------------
if (small_file) then begin
    data = 0 & pol_data = 0
endif else begin
    data_plot = temporary(data)
    pol_data = 0
    find_min_max_valid, grid, mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
endelse

;-----------------------------------
; export in FITS and as an array the original gnomic map before alteration
;-----------------------------------

; grid -> IDL array
if arg_present(map_out) then map_out = proj2map_out(grid, bad_data=bad_data)

; grid -> FITS file
if keyword_set(fits) then begin 
    proj2fits, grid, fits, $
               projection = 'GNOM', flip=flip, $
               rot = rot_ang, coord=coord_out, reso = resgrid*60., unit = sunits, min=mindata, max = maxdata
endif

; -------------------------------------------------------------
; set min and max and computes the color scaling
; -------------------------------------------------------------
if (small_file) then begin

endif else begin
    if (do_poldirection) then begin
        min_set = 0.
        max_set = 2*!pi
    endif
    if (truetype eq 2) then begin
        ; truecolors=2 map each field to its color independently
        color = bytarr(xsize,ysize,zsize)
        for i=0,zsize-1 do begin
            find_min_max_valid, grid[*,*,i], mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
            color[0,0,i] = COLOR_MAP(grid[*,*,i], mindata, maxdata, Obs, $
                          color_bar = color_bar, mode=mode_col, silent=silent)
        endfor
    endif else begin
        ; same for truecolors=1 and false colors:
        color = COLOR_MAP(grid, mindata, maxdata, Obs, $
                          color_bar = color_bar, mode=mode_col, $
                          minset = min_set, maxset = max_set, silent=silent)
    endelse
        
    if (do_polvector) then begin ; rescale polarisation vector in each valid pixel
        planvec[*,*,0] = vector_map(planvec[*,*,0], Obs, vector_scale = vector_scale)
    endif
    Obs = 0
    grid = 0
    Tmin = mindata & Tmax = maxdata
endelse

return
end

