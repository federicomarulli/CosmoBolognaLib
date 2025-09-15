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
pro oplot_graticule, graticule, eul_mat, $
                     projection=projection, $
                     mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, orthographic=orthographic, $
                     flip = flip, _extra = oplot_kw, $
                     half_sky=half_sky, coordsys=coordsys, charsize=charsize, reso_rad=reso_rad, $
                     latex=latex, ltxstc=ltxstc
;+
; NAME:
;       OPLOT_GRATICULE
;
; PURPOSE:
;       overplots graticule (ie, spherical coordinates grid) on top
;       of existing map
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;       oplot_graticule, graticule, eul_mat, $
;          [projection=,mollweide=,gnomic=,cartesian=, $
;           orthographic=,flip=,half_sky=,coordsys=, reso_rad=, + all oplot keywords]
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
;         Feb 2003, corrected bug occuring when combining astrophysical coordinates
;         changes with arbitrary rotation
;-

identify_projection, projtype, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, orthographic=orthographic

if keyword_set(flip) then flipconv=1 else flipconv = -1  ; longitude increase leftward by default (astro convention)

do_ig = (n_elements(coordsys) eq 2) 

; define default graticule spacing, in Degrees
if (projtype eq 2) then begin ; gnom
    dlong = 5. 
    dlat  = 5. 
    gratmin = 0.
endif
if (projtype eq 1) then begin ; moll
    dlong = 45.
    dlat  = 45.
    gratmin = 10.
endif
if (projtype eq 3) then begin ; cart
    dlong = 45. 
    dlat  = 45. 
    gratmin = 2.
endif
if (projtype eq 4) then begin ; ortho
    dlong = 45.
    dlat  = 45.
    gratmin = 10.
endif

; read in user defined grid spacings
if (n_elements(graticule) eq 2 and min(graticule) ge gratmin) then begin
    dlong = float(graticule[0]) & dlat = float(graticule[1])
endif else begin
    if (graticule[0] gt gratmin) then begin 
        dlong = float(graticule[0]) & dlat = dlong
    endif
endelse

fsgrat = (defined(reso_rad)) ? !dtor/(reso_rad+1.d-11) : 1. ; number of pixels / degree
fsgrat = long(fsgrat) > 1 < 5
; define variables
epsilon = 1.d-5
nmerid = fix(181./dlong)
nparal = fix(90./dlat)-1
nv = (projtype eq 2 || projtype eq 3) ? 721*fsgrat : 361 ; more points for partial projections
vector = DINDGEN(nv)/(nv-1.d0) * (1.d0-2.d0*epsilon) + epsilon ; from epsilon to 1-epsilon
bounds = [[-nmerid,nmerid-1],[-nparal,nparal]]

do_rot = (n_elements(eul_mat) eq 9)
form = '(i4)'
if (abs(round(dlong)-dlong) gt 1.e-2 || abs(round(dlat)-dlat) gt 1.e-2) then form = '(f6.1)'

case projtype of
2: begin  ; gnomic : straightforward
    for jg=0,1 do begin
        for i = bounds[0,jg],bounds[1,jg] do begin
;             if (jg eq 0) then ang2vec, vector*!pi,        replicate(i*dlong*!DtoR,nv), vv ; meridians
;             if (jg eq 1) then ang2vec, replicate((90.-i*dlat)*!DtoR,nv), vector*(2.*!pi), vv ; parallels
            if (jg eq 0) then begin
                mylong = i*dlong ; longitude in Deg
                linelabel = strtrim(string(mylong,form=form),2)
                ang2vec, vector*!pi,        replicate(mylong*!DtoR,nv), vv ; meridians
            endif
            if (jg eq 1) then begin
                mylat = i*dlat ; latitude in Deg
                linelabel = strtrim(string(mylat,form=form),2)
                ang2vec, replicate((90.-mylat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            endif
            if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
            if (do_rot) then vv = vv # transpose(eul_mat)

            k = where(vv(*,0) gt 0, nk)
            if (nk gt 0) then begin
                u = vv(k,1)/vv(k,0)
                v = vv(k,2)/vv(k,0)
                good = where(abs(u) lt !x.crange[1]*1.1 and abs(v) lt !y.crange[1]*1.1 ,ng)
                ; reorder points to have one continuous segment across the plot
                bad = where(good-shift(good,1) ne 1, nbad)
                if (nbad gt 1) then good = shift(good, bad[1])
;                oplot, flipconv * u, v, _extra = oplot_kw
               if (ng gt 1) then oplot_sphere, flipconv *u[good], v[good], _extra = oplot_kw, $
                 linelabel=linelabel,/flush, charsize=charsize[jg], latex=latex, ltxstc=ltxstc
            endif
        endfor
    endfor
end
1: begin  ; mollweide : deal with boundaries
    for jg=0,1 do begin
        for i = bounds[0,jg],bounds[1,jg] do begin
            if (jg eq 0) then begin
                mylong = i*dlong ; longitude in Deg
                linelabel = strtrim(string(mylong,form=form),2)
                ang2vec, vector*!pi,        replicate(mylong*!DtoR,nv), vv ; meridians
            endif
            if (jg eq 1) then begin
                mylat = i*dlat ; latitude in Deg
                linelabel = strtrim(string(mylat,form=form),2)
                ang2vec, replicate((90.-mylat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            endif
            if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
            if (do_rot) then vv = vv # transpose(eul_mat)

            vec2moll, vv, u, v
            oplot_sphere, -flipconv * u, v, _extra = oplot_kw, linelabel=linelabel, charsize=charsize[jg], $
                          latex=latex, ltxstc=ltxstc
;;            oplot_sphere, flipconv * u, v, _extra = oplot_kw
        endfor
    endfor
end
4: begin  ; orthographic : deal with boundaries
    if keyword_set(half_sky) then begin
        nd = 1 ; number of half-sky disc
        c0 = 0
    endif else begin
        nd = 2
        c0 = 1
    endelse
    for jg=0,1 do begin
        for i = bounds[0,jg],bounds[1,jg] do begin
            if (jg eq 0) then begin
                mylong = i*dlong ; longitude in Deg
                linelabel = strtrim(string(mylong,form=form),2)
                ang2vec, vector*!pi,        replicate(mylong*!DtoR,nv), vv ; meridians
            endif
            if (jg eq 1) then begin
                mylat = i*dlat ; latitude in Deg
                linelabel = strtrim(string(mylat,form=form),2)
                ang2vec, replicate((90.-mylat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            endif
            if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
            if (do_rot) then vv = vv # transpose(eul_mat)

            for sign = 1,1-nd,-2 do begin ; either (1,-1) or (1)
                k = where(vv[*,0]*sign ge 0, nk)
                if (nk gt 0) then begin
                    u = vv[k,1]
                    v = vv[k,2]
                    ;oplot_sphere, flipconv*(u+c0)*sign, v, _extra = oplot_kw
                    oplot_sphere,  flipconv*(u+c0)*sign, v, _extra = oplot_kw, linelabel=linelabel, charsize=charsize[jg], $
                                   latex=latex, ltxstc=ltxstc
                endif ; nk>0
            endfor ; loop on sign
        endfor
    endfor
end
3: begin  ; cartesian : straightforward
    for jg=0,1 do begin
        for i = bounds[0,jg],bounds[1,jg] do begin
;             if (jg eq 0) then ang2vec, vector*!pi,        replicate(i*dlong*!DtoR,nv), vv ; meridians
;             if (jg eq 1) then ang2vec, replicate((90.-i*dlat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            if (jg eq 0) then begin
                mylong = i*dlong ; longitude in Deg
                linelabel = strtrim(string(mylong,form=form),2)
                ang2vec, vector*!pi,        replicate(mylong*!DtoR,nv), vv ; meridians
            endif
            if (jg eq 1) then begin
                mylat = i*dlat ; latitude in Deg
                linelabel = strtrim(string(mylat,form=form),2)
                ang2vec, replicate((90.-mylat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            endif
            if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
            if (do_rot) then vv = vv # transpose(eul_mat)

            phi = atan(vv[*,1],vv[*,0]) ; in [0,2pi]
	    theta = asin(vv[*,2])       ; in [0,pi]
            ; OPLOT,-flipconv*phi,theta, _extra = oplot_kw
;            oplot_sphere, flipconv * phi, theta, _extra = oplot_kw
;            oplot_sphere, flipconv* phi, theta, _extra = oplot_kw,
;            linelabel=linelabel, charsize=charsize
            good = where(abs(phi) lt !x.crange[1]*1.1 and abs(theta) lt !y.crange[1]*1.1 ,ng)
                                ; reorder points to have one continuous segment across the plot
            bad = where(good-shift(good,1) ne 1, nbad)
            if (nbad gt 1) then good = shift(good, bad[1])
            if (ng gt 1) then oplot_sphere, flipconv *phi[good], theta[good], _extra = oplot_kw, $
              linelabel=linelabel,charsize=charsize[jg], latex=latex, ltxstc=ltxstc
        endfor
    endfor
end
endcase


;-----------------------
;     FOR i=-nmerid,nmerid-1 DO begin  ; meridians
;         ang2vec, vector*!pi, replicate(i*dlong*!DtoR,n_elements(vector)), vv
;         if (do_rot) then vv = vv # transpose(eul_mat)

;         k = where(vv(*,0) gt 0, nk)
;         if (nk gt 0) then begin
;             u = vv(k,1)/vv(k,0)
;             v = vv(k,2)/vv(k,0)
;             OPLOT, flipconv * u, v, COLOR = !P.COLOR
;         endif
;     endfor
;     FOR i=-nparal,nparal DO begin  ; parallels
;         ang2vec, replicate((90.-i*dlat)*!DtoR,n_elements(vector)), vector*2*!pi, vv
;         if (do_rot) then vv = vv # transpose(eul_mat)

;         k = where(vv(*,0) gt 0, nk)
;         if (nk gt 0) then begin
;             u = vv(k,1)/vv(k,0)
;             v = vv(k,2)/vv(k,0)
;             OPLOT, flipconv * u, v, COLOR = !P.COLOR
;         endif
;     endfor
; endif

; if (do_moll) then begin
;     ; mollweide : deal with boundaries
;     FOR i=-nmerid,nmerid-1 DO begin  ; meridians
;         ang2vec, vector*!pi, replicate(i*dlong*!DtoR,nv), vv
;         if (do_rot) then vv = vv # transpose(eul_mat)

;         vec2moll, vv, u, v
;         bad = where(abs(u-shift(u,1)) gt .1, nbad)
;         if (nbad eq 0) then begin
;             OPLOT, flipconv * u, v, COLOR = !P.COLOR
;         endif else begin
;             bad = [0,bad,n_elements(u)-1]
;             for j=0,nbad do begin
;                 if (bad[j+1] gt bad[j]) then $
;                   oplot, flipconv * u[bad[j]:bad[j+1]-1], v[bad[j]:bad[j+1]-1], color=!p.color
;             endfor
;         endelse
;     endfor
;     FOR i=-nparal,nparal DO begin  ; parallels
;         ang2vec, replicate((90.-i*dlat)*!DtoR,nv), vector*2*!pi, vv
;         if (do_rot) then vv = vv # transpose(eul_mat)

;         vec2moll, vv, u, v
;         bad = where(abs(u-shift(u,1)) gt .1, nbad)
;         if (nbad eq 0) then begin
;             OPLOT, flipconv * u, v, COLOR = !P.COLOR
;         endif else begin
;             bad = [0,bad,n_elements(u)-1]
;             for j=0,nbad do begin
;                 if (bad[j+1] gt bad[j]) then $
;                   oplot, flipconv * u[bad[j]:bad[j+1]-1], v[bad[j]:bad[j+1]-1], color=!p.color
;             endfor
;         endelse
;     endfor
; endif

return
end
