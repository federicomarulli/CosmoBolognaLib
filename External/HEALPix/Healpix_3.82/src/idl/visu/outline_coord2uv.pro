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
pro ocuv_assert_scalar, variable, name, error=error

if (undefined(error)) then error = 0
nv = n_elements(variable) 
if (nv gt 1) then begin
    print,'Error: Argument '+name+' in OUTLINE structure has '+string(nv)+' elements.'
    print,'It should be a scalar'
    error = 1
endif

return
end

pro outline_coord2uv, outline, coord_out, eul_mat, outline_uv, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, orthographic=orthographic, flip = flip, show=show, half_sky=half_sky, thick_default = thick_default, _extra = oplot_kw
;+
; NAME:
;     outline_coord2uv
;
;
; PURPOSE:
;     turn a outline in astro coordinate, into a outline in uv coordinate
;     both input and output are structures,
;     if Show is set, overplot the result on current plot
;     Note : the outline may be open, the procedure will NOT close it
;
; CATEGORY:
;     internal processing for data visualisation
;
;
; CALLING SEQUENCE:
;     OUTLINE_COORD2UV, Outline, Coord_Out, Eul_Mat, [Outline_Uv]
;           Projection=,Mollweide=, Gnomic=, Cartesian=, Orthographic=,Flip=,Show=,half_sky
;
; 
; INPUTS:
;     Outline : structure containing
;       - (COO*) coordinate
;       - (RA*) RA (or longitude) vector
;       - (DEC*) Dec (or latitude) vector
;       - (LINE*) line style 
;          +2 : black dashes
;          +1 : black dots
;           0 : black solid
;          -1 : black dots on white background
;          -2 : black dashes on white background
;          (default = 0 : solid)
;       - (PSY*) symbol used to represent vertices of outline
;            (same meaning as standard PSYM in IDL,
;             if 9<=abs(psym)<=46, D. Fanning's cgSYMCAT symbols
;             definition will be used, for example psym=9 is an open circle)
;          if <=0, the vertices are represented with the chosen symbols, and
;          connected by arcs of geodesics.
;          if >0, only the vertices are shown
;          (default = 0)
;       - (SYM*) symbol size (same meaning as SYMSIZE in IDL)
;       - (COL*) 1 or 2 color index in [0,255]
;
;     Coord_Out, Eul_Mat
;
;
; KEYWORD PARAMETERS:
;    Mollweide =
;    Gnomic =
;    Cartesian =
;    Orthographic =
;    Flip = 
;    Projection = either 'MOLL*', 'GNOM*', 'CART*' or 'ORTH*'
;
; OPTIONAL OUTPUTS:
;    Outline_Uv
;
;
; COMMON BLOCKS:
;    none
;
;
; SIDE EFFECTS:
;    oplot contour on existing plot if Show is set
;
;
; PROCEDURE:
;   makes a UV projection in either Mollweide or Gnomic scheme
;
;   calls: datatype, vec2moll
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;       2002-03-21, E. Hivon, IPAC, v1.0
;       2004-11-09, E. Hivon, IPAC : cleaning, bug correction (for orthview), 
;                                    addition of point plotting (ie, not
;                                    connected by curves)
;       2007-03-19: can now deal with single point, exits gracefully if
;       PSYM,LINES,SYMSIZE are not scalar
;       2007-05-14: test for ambiguous tag name in outline structure
;       2011-01-12: calls SYMCAT if 9 <= ABS(psym) <= 46
;       2013-02-08: replaced SYMCAT with CGSYMCAT
;       2014-06-23: allow different sub-structure to use different customized
;         symbols (via CGSYMCAT when PSYM in [9,46])
;       2018-??: introduced COLOR
;       2018-11-30: allow for large number of substructures
;       2018-12-04: avoid creation of intermediate structure, 
;          slightly changed number of interpolated points
;-

; outline : from astro coordinate to uv plan

identify_projection, projtype, projection=projection, mollweide=mollweide, gnomic=gnomic

flipconv = keyword_set(flip) ? 1 : -1  ; longitude increase leftward by default (astro convention)

if (datatype(outline) ne 'STC') then return
n_outlines = n_tags(outline)
no_sub = 0
if (datatype(outline.(0)) ne 'STC') then begin
    n_outlines = 1
    no_sub = 1
endif
do_show   = keyword_set(show)
do_outstc = arg_present(outline_uv)

ist = 0L
for i=0L,n_outlines-1 do begin
    if (no_sub eq 1) then c1 = outline else c1 = outline.(i)
    ;--------- parse structure fields ----------
    nc1 = n_tags(c1)
    names = tag_names(c1)
;     cont_coord = strmid(strupcase(c1.(0)),0,1)
;     cont_ra = c1.(1)
;     cont_dec = c1.(2)
    invalid_tags = replicate(1, nc1)
    expected_tags = ['Required:       COORD RA DEC','Optional:       LINESTYLE PSYM SYMSIZE']

    iw = index_word(names,'COO',err=errword)
    if (errword eq 0) then begin
        cont_coord = strupcase(strmid(c1.(iw),0,1))
        nc1 = nc1 - 1
        invalid_tags[iw] = 0
    endif
    ocuv_assert_scalar,cont_coord,'COORD', error=err_assert

    iw = index_word(names,'RA',err=errword)
    if (errword eq 0) then begin
        cont_ra = c1.(iw)
        nc1 = nc1 - 1
        invalid_tags[iw] = 0
    endif

    iw = index_word(names,'DEC',err=errword)
    if (errword eq 0) then begin
        cont_dec = c1.(iw)
        nc1 = nc1 - 1
        invalid_tags[iw] = 0
    endif

    if (undefined(cont_coord) or undefined(cont_ra) or undefined(cont_dec)) then begin
        print,'Invalid tag names in structure provided for outline :',names
        print,expected_tags,form='(a)'
        message,/info,'No outline is plotted'
        return
    endif

    cont_type = 0
    iw = index_word(names,'LINE',err=errword)
    if (errword eq 0) then begin
        cont_type = c1.(iw)
        nc1 = nc1 - 1
        invalid_tags[iw] = 0
    endif
    ocuv_assert_scalar,cont_type,'LINESTYLE', error=err_assert

    sym_type = 0
    iw = index_word(names,'PSY',err=errword)
    if (errword eq 0) then begin
        sym_type = c1.(iw)
        nc1 = nc1 - 1
        invalid_tags[iw] = 0
    endif
    ocuv_assert_scalar,sym_type,'PSYM', error=err_assert
;     if (abs(sym_type) ge 9 && abs(sym_type) le 46) then begin
;         mysign = (sym_type lt 0) ? -1 : 1
;         sym_type = mysign * cgsymcat(abs(sym_type)) ; call usersym and set sym_type to 8 or -8
;     endif
    
    ssize = 1
    iw = index_word(names,'SYM',err=errword)
    if (errword eq 0) then begin
        ssize = c1.(iw)
        nc1 = nc1 - 1
        invalid_tags[iw] = 0
    endif
    ocuv_assert_scalar,ssize,'SYMSIZE', error=err_assert


    colors = [!p.color, !p.background]
    iw = index_word(names,'COL',err=errword)
    if (errword eq 0) then begin
        colors[0] = c1.(iw)[0]
        if n_elements(c1.(iw)) ge 2 then colors[1] = c1.(iw)[1]
        nc1 = nc1 - 1
        invalid_tags[iw] = 0
    endif

    thick_factor = 1.0
    iw = index_word(names,'THI',err=errword)
    if (errword eq 0) then begin
        thick_factor = c1.(iw)
        nc1 = nc1 - 1
        invalid_tags[iw] = 0
    endif
    ocuv_assert_scalar,thick_factor,'THICK', error=err_assert
    ;;;;;;print,'THICK: ',thick_factor,thick_default,!P.thick
    
    if (err_assert ne 0) then begin
        message,/info,'No outline is plotted'
        return
    endif

    if (nc1 lt 0) then begin
        print,'Duplicated/ambiguous tag names in structure provided for outline :',names
        print,expected_tags,form='(a)'
        message,/info,'No outline is plotted'
        return
    endif
    bad_tags = where(invalid_tags, nc1)
    if (nc1 gt 0) then begin
        print,'WARNING: Unknown/Duplicated tag name(s) found in structure provided for outline :',names[bad_tags]
        print,'         will be ignored.'
        print,expected_tags,form='(a)'
    endif

    ; angles -> vector
    if (n_elements(cont_ra) ne n_elements(cont_dec)) then begin
        message,'Inconsistent coordinates in outline structure'
    endif
    ang2vec, /astro, cont_dec, cont_ra, cont_vec1

    ktype_min = 0
    ktype_max = 1
    if (sym_type eq 0) then ktype_max = 0 ; outline only
    if (sym_type gt 0) then ktype_min = 1 ; vertices only
    for ktype = ktype_min, ktype_max do begin

              ; interpolate intermediate points between vertices (along geodesics)
              ; this will avoid big step between vertices to be interpreted as boundary effect
              ; on the uv plan
        np_in = n_elements(cont_vec1[*,0])
        if (ktype eq 0) then begin
            pcont_type = cont_type
            psym_type = 0
            rf = 20
            if (np_in gt 500) then rf = 10
            ;np_out = rf*np_in + 1 ! 2018-12-04
            np_out = rf*(np_in-1L) + 1
        endif else begin
            pcont_type = 0
            psym_type = abs(sym_type)
            rf = 1             ; no interpolation when plotting vertices
            np_out = rf*np_in
        endelse
        if (np_in gt 1) then begin
            cont_x = interpol(cont_vec1[*,0], np_out) ; linear interpolation in X
            cont_y = interpol(cont_vec1[*,1], np_out) ; .. Y
            cont_z = interpol(cont_vec1[*,2], np_out) ; .. Z
        endif else begin
            one_tmp = replicate(1, np_out)
            cont_x = cont_vec1[0,0] * one_tmp
            cont_y = cont_vec1[0,1] * one_tmp
            cont_z = cont_vec1[0,2] * one_tmp
        endelse
        cont_norm = sqrt( cont_x*cont_x + cont_y*cont_y + cont_z*cont_z ) ; norm of vectors 
        cont_vec = [[cont_x/cont_norm],[cont_y/cont_norm],[cont_z/cont_norm]]
                                ; go from contour coordinates to plot coordinates
        cont_vec = rotate_coord(cont_vec, inco= cont_coord, outco=coord_out)
        cont_vec = rotate_coord(cont_vec, euler_matrix = transpose(eul_mat))

                                ; vector -> (U,V)
        case projtype of
            1: begin            ; mollweide
                vec2moll, cont_vec, u_cont, v_cont, flip=flip
                nkeep = n_elements(u_cont)
            end
            2: begin            ; gnomic
                keep = where(cont_vec[*,0] gt 0., nkeep)
                if (nkeep gt 0) then begin
                    u_cont =  flipconv*cont_vec[keep,1]/cont_vec[keep,0] ; U=Y/X
                    v_cont =           cont_vec[keep,2]/cont_vec[keep,0] ; U=Z/X
                endif
            end
            3: begin            ; cartesian
                u_cont = flipconv * atan(cont_vec[*,1],cont_vec[*,0])
                v_cont =            asin(cont_vec[*,2])
                nkeep = n_elements(u_cont)
            end
            4: begin            ; orthographic
                if (keyword_set(half_sky)) then begin
                    good = where(cont_vec[*,0] ge 0, nkeep)
                    if (nkeep) gt 0 then begin
                        u_cont = flipconv *  cont_vec[good,1]
                        v_cont =             cont_vec[good,2]
                    endif
                endif else begin
                    s_cont = (cont_vec[*,0] gt 0)*2 - 1 ; +1 for x>0, -1 for x<0
;                   u_cont = flipconv * (cont_vec[*,1] + s_cont)
                    u_cont = flipconv * (cont_vec[*,1] + 1)*s_cont
                    v_cont =             cont_vec[*,2]
                    nkeep = n_elements(u_cont)
                endelse
            end
            else: begin
                message,'Unsupported projection: '+projection
            end
        endcase
               
        ; make plot and/or create output structure
        if (nkeep gt 0) then begin
            if do_outstc then begin ;  store (U,V) in structure
                c1uv = {U:u_cont, V:v_cont, T:pcont_type, ST:psym_type, SS:ssize, COL:colors}
                tag = 's'+string(ist,form='(i6.6)')
                if (ist eq 0) then begin
                    outline_uv = create_struct(tag,c1uv) 
                endif else begin
                    outline_uv = create_struct(outline_uv,tag,c1uv)
                endelse
                ist += 1
            endif 
            if do_show then begin
                if (psym_type ge 9 && psym_type le 46) then  psym_type = cgsymcat(psym_type) ; call usersym and set psym_type to 8
                oplot_sphere, u_cont, v_cont, line_type=pcont_type, _extra = oplot_kw, psym=psym_type, symsize=ssize, color=colors, thick= thick_default * thick_factor
            endif
        endif
    endfor ; loop on type
endfor ; loop on outline


return
end

