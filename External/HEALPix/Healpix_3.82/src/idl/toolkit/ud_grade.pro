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
; function prograde_qu, map_in, nside_in, nside_out

; ; number of small pixels within the big ones
; rat2 = (nside_out/nside_in)^2
; ramp = lindgen(rat2)
; one  = replicate(1,rat2)

; npix = n_elements(map_in[*,0])
; ncol = n_elements(map_in[0,*])

; qu = [0,1]
; do_t = 0
; if (ncol eq 3) then begin
;     qu = [1,2]
;     do_t = 1
; endif

; map_out = make_array(npix*rat2, ncol, type=size(map_in,/type))
; ; build templates
; ntemplates = nside2ntemplates(nside_in)
; print,ntemplates

; bit_manipulation

; for it = 0, ntemplates-1 do begin
;     ; find all pixels having the same template
;     same_shape_pixels_nest, nside_in, it, list, refl, nr=nr
;     print,it,nr

;     ; build (Q,U) transform for template pixel
;     p0 = list[0]
;     r0 = refl[0]
;     small_pix = p0*rat2 + ramp

;     pix2vec_nest, nside_in,  p0,        vec_big
;     pix2vec_nest, nside_out, small_pix, vec_small

;     sum = vec_small + vec_big ## one
;     vz = sum[*,2]
;     s2 = total(sum*sum,2)
;     sindphi = vec_small[*,0]*vec_big[0,1]-vec_small[*,1]*vec_big[0,0]
;     sindphi = sindphi / sqrt( (1.d0-vec_small[*,2]^2) * (1.d0-vec_big[0,2]) )
;     sinpsi  = - 0 * sindphi * 2.d0 * vz / s2 ; <<<<<<<<<<<<<<<
;     sin2psi = 2.d0 * sinpsi * sqrt( 1.d0 - sinpsi^2 )
;     cos2psi = 1.d0 - 2.d0 * sinpsi * sinpsi
;     sum = 0 & vz = 0 & s2 = 0 & sindphi = 0 & sinpsi = 0

;     ; apply transform to template pixel
;     if (do_t) then map_out[small_pix,0] = map_in[p0,0]
;     map_out[small_pix,qu[0]] = map_in[p0,qu[0]]*cos2psi + map_in[p0,qu[1]]*sin2psi
;     map_out[small_pix,qu[1]] = map_in[p0,qu[1]]*cos2psi - map_in[p0,qu[0]]*sin2psi

;     ; apply transform to all other pixels with that template
;     for i=1, nr-1 do begin

;         p1 = list[i]
;         r1 = refl[i]
;         ; reshufle small pixels according to reflexions
;         ;          mapping template onto big pixel
;         rf = r1 XOR r0
;         case rf of
;             0: rp = ramp
;             1: rp = swapLSBMSB(ramp)
;             2: rp = invswapLSBMSB(ramp) + rat2
;             3: rp = invLSBMSB(ramp)
;         end
;         small_pix = p1*rat2 + rp
;         sgn = 1
;         if (rf eq 1 or rf eq 2) then sgn = -1

;         if (do_t) then map_out[small_pix,0] = map_in[p1,0]
;         map_out[small_pix,qu[0]] = map_in[p1,qu[0]]*cos2psi + (map_in[p1,qu[1]]*sgn)*sin2psi
;         map_out[small_pix,qu[1]] = map_in[p1,qu[1]]*cos2psi - (map_in[p1,qu[0]]*sgn)*sin2psi
;     endfor
; endfor

; return, map_out
; end




function sub_ud_grade, map_in, nside_in, nside_out, bad_data=bad_data, power = power, pessimistic=pessimistic
;
; sub_ud_grade : does the effective up/degradation of a nested scheme map
;
; map_out = map_in * (Nside_out/Nside_in)^power
;  intensitive quantity (eg. Temperature): power = 0
;  extensive quantity (number of hits) : power = -2
;                   rms                : power = 1
; power is 0 by default
;

do_pess = keyword_set(pessimistic)
npix_in  = nside2npix(nside_in)
npix_out = nside2npix(nside_out)
nmaps = n_elements(map_in[0,*])

if undefined(bad_data) then bad_data = !healpix.bad_value
;bad = where(map_in eq bad_data, nbad)
; test for bad pixels less sensitive to input map precision (2010-02-17)
threshold = abs(1.e-7 * bad_data)
bad = where(abs(map_in - bad_data) lt threshold  or finite(map_in,/nan), nbad)

ratio = 1
if defined(power) then ratio = (float(nside_out)/float(nside_in)) ^ power 

if nside_out gt nside_in then begin ; upgrade resolution (more, smaller pixels)
    rat2 = npix_out/npix_in
    if nbad gt 0 then map_in[bad] = !values.f_nan

    map_out = reform(replicate(ratio,rat2)#map_in[*],npix_out,nmaps)

;;    bad = where(map_out eq !values.f_nan, nbad)
    bad = where( finite(map_out, /nan), nbad)  ; corrected 2009-05-07
    if nbad gt 0 then map_out[bad] = bad_data
endif else begin ; degrade resolution (fewer, larger pixels)
    rat2 = npix_in/npix_out
    hit = replicate(1B,rat2,npix_out*nmaps)
    if (nbad gt 0) then begin
        map_in[bad] = 0.
        hit[bad] = 0B
    endif
    hit_out = total(hit,1) ; number of good small pixels in each large one
    if (do_pess) then begin
        bad = where(hit_out ne rat2, nbad)
    endif else begin
        bad = where(hit_out eq 0, nbad)
    endelse
    hit_out = (hit_out > 1B) ; to allow division
    map_in = reform(map_in, rat2, npix_out,nmaps, /over)
    map_out = total(map_in, 1) * ratio / hit_out
    if (nbad gt 0) then map_out[bad] = bad_data
endelse

return,map_out
end
;_____________________________________________________________________
;
pro ud_grade, map_in, map_out, nside_out=nside_out, order_in=order_in, order_out=order_out, bad_data = bad_data, pessimistic = pessimistic, help=help
;+
; NAME:
;  ud_grade
;
; PURPOSE:
;  upgrade or degrade a full sky Healpix map
;  the data in the map is supposed to be intensive 
;  (doesn't scale with pixel area) just like temperature
;
; CATEGORY:
;  Healpix toolkit
;
; CALLING SEQUENCE:
;  ud_grade, map_in, map_out, [bad_data=, help=, nside_out=, order_in=, order_out=, pessimistic=]
;
; INPUTS:
;  map_in : either 
;        - a fits file 
;            (with either the temperature only or temperature + polarisation(Q,U)
;            or the DMR like information : PIXEL, SIGNAL, N_OBS and SERROR, with
;            their specific scaling law),
;        - a memory vector containing a Healpix full sky map (T or U or Q)
;         
;  map_out : upgraded or degraded map
;        - if map_in is a FITS file, map_out should be a string containing the
;          name of the output file
;        - if map_in is a vector, map_out should be an IDL variable
;
; KEYWORD PARAMETERS:
;   bad_data : flag value of missing pixel
;        default : !healpix.bad_value = -1.6375e30
;
;   help: if set, this documentation header is printed out and the code exits
;
;   nside_out : output resolution parameter (can be larger or smaller than the
;        input one)
;        default : same as input (map unchanged)
;
;   order_in : input map ordering (either 'RING' or 'NESTED')
;        default : if map_in is a Healpix FITS file, use information in header
;                  if map_in is a vector, order_in has to be defined
;
;   order_out : output map ordering (either 'RING' or 'NESTED')
;        default : same as order_in
;
;   pessimistic : if set to 1, during degradation each big pixel containing one
;             bad or missing small pixel is also considered as bad,
;                 if set to 0, each big pixel containing at least one good pixel
;             is considered as good (optimistic)
;             default = 0 (:optimistic)
;
; OUTPUTS:
;   none
;
; OPTIONAL OUTPUTS:
;   none
;
; COMMON BLOCKS:
;   none
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;   reorder to NEST before up/degrading
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;    2000-01-08, Eric Hivon, Caltech
;    2000-11   : deal with cut sky format files, added pessimistic
;    2006-11-20: check for validity of order_in and order_out
;    2009-05-07: correctly flags bad output pixels with bad_data value when
;               upgrading maps
;    2010-02-17: bad pixels correctly identified in DP map
;    2012-03-16: added HELP keyword
;    2017-05-30: replaced obsolete DATATYPE() with IDL's SIZE(/TNAME)
;    2020-02-19: process partial sky FITS files containing TQU data
;-

routine = 'UD_GRADE'

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if n_params() ne 2 then begin
    print,'syntax : '+routine+', map_in, map_out, '
    print,'  [bad_data=, help=, nside_out=, order_in=, order_out=, pessimistic=]'
    print
    return
endif

dtype_in  = size(/tname,map_in)
dtype_out = size(/tname,map_out)

if (dtype_in eq 'STRING' && ~(dtype_out eq 'STRING' || dtype_out eq 'UNDEFINED')) then begin
    print,'inconsistent type for map_in and map_out in '+routine
    print,'******** ABORT *********'
    return
endif

if (dtype_in ne 'STRING' && dtype_out eq 'STRING') then begin
    print,'inconsistent type for map_in and map_out in '+routine
    return
endif

if (dtype_in ne 'STRING' && undefined(order_in)) then begin
    print,'the input ordering has to be defined for online input map'
    print,'******** ABORT *********'
    return
endif

; ----- test resolution parameters Nside ---------
if defined(nside_out) then begin
    npix_out = nside2npix(nside_out,error = error)
    if (error ne 0) then begin
        print,'invalid nside_out = ',nside_out,' in '+routine
        print,'******** ABORT *********'
        return
    endif
endif
;----------------
; file to file
;----------------
if dtype_in eq 'STRING' then begin
    junk = getsize_fits(map_in, ordering = order_in, nside = nside_in, type = ftype_in)
    fits_info, map_in, /silent, n_ext = n_ext

    if ~(keyword_set(order_out) || arg_present(order_out)) then order_out = order_in
    ord_in  = decode_ordering(order_in,  full=ford_in) ; strupcase(strmid(order_in,0,4))
    ord_out = decode_ordering(order_out, full=ford_out) ; strupcase(strmid(order_out,0,4))

    npix_in = nside2npix(nside_in, err = error)
    if (error ne 0) then begin
        print,'wrong Nside in '+map_in
        print,'******** Abort *********'
        return
    endif
    if undefined(nside_out) then nside_out = nside_in

    if (ftype_in eq 3) then begin
        ; ------- cut sky -------
        names = sxpar(headfits(map_in, ext=1),'TTYPE*')
        cut4  = array_equal(strupcase(strmid(names[0:3],0,1)),['P','S','N','S'])
        if (cut4) then begin
            for i_ext = 0, n_ext-1 do begin
                if (n_ext gt 1) then print,'Extension: ',i_ext
                read_fits_cut4, map_in, pixel, signal, n_obs, serror, xh = xhdr, ext= i_ext
                                ; update header
                sxaddpar,          xhdr,'NSIDE',nside_out
                add_ordering_fits, xhdr, ordering=ord_out
                sxaddpar,          xhdr,'HISTORY','    PROCESSED BY '+routine
                if (nside_out ne nside_in) then sxaddpar,          xhdr,'HISTORY', $
                  string('       NSIDE: ',nside_in,' -> ',nside_out,form='(a,i5,a,i5)')
                if (ord_in ne ord_out) then sxaddpar,          xhdr,'HISTORY', $
                  '    ORDERING: '+ford_in+' -> '+ford_out
                                ; process data
                ud_grade_cut4, pixel, signal, n_obs, serror, nside_in = nside_in, $
                               nside_out = nside_out, order_in = ord_in, order_out = ord_out                
                                ; write file
                write_fits_cut4, map_out, pixel, signal, n_obs, serror, xh = xhdr, ext= i_ext
            endfor
        endif else begin
            if undefined(bad_data) then bad_data = !healpix.bad_value
            for i_ext = 0, n_ext-1 do begin
                if (n_ext gt 1) then print,'Extension: ',i_ext
                read_fits_partial, map_in, pixel, tmp, xh = xhdr, ext= i_ext
                                ; update header
                sxaddpar,          xhdr,'NSIDE',nside_out
                add_ordering_fits, xhdr, ordering=ord_out
                sxaddpar,          xhdr,'HISTORY','    PROCESSED BY '+routine
                if (nside_out ne nside_in) then sxaddpar,          xhdr,'HISTORY', $
                  string('       NSIDE: ',nside_in,' -> ',nside_out,form='(a,i5,a,i5)')
                if (ord_in ne ord_out) then sxaddpar,          xhdr,'HISTORY', $
                  '    ORDERING: '+ford_in+' -> '+ford_out
                                ; process data
                ncol = (size(tmp,/dim))[1]  & mtype = size(tmp,/type)
                map1 = make_array(npix_in,  ncol, type=mtype)
                map1[*] = bad_data
                map1[pixel,0:ncol-1] = tmp
                map2 = make_array(npix_out, ncol, type=mtype)
                index_in  = reorder(lindgen(npix_in), in=ord_in,out='NEST')
                index_out = reorder(lindgen(npix_out),in='NEST',out=ord_out)
                for i=0,ncol-1 do begin
                    map2[0,i] = sub_ud_grade(map1[index_in,i], nside_in, nside_out, bad_data=bad_data, pessimistic=pessimistic)
                    map2[0,i] = map2[index_out,i]
                endfor
                pixel_out = where(map2[*,0] ne bad_data and finite(map2[*,0]), ngood)
                                ; write file
                write_fits_partial, map_out, pixel_out, map2[pixel_out,0:ncol-1], xh=xhdr, ext = i_ext
            endfor
        endelse
    endif else if (ftype_in eq 0 || ftype_in eq 2) then begin
        if (n_ext gt 1) then message,'Can not currently deal with multiple extension in this format'
        ; ----- full sky -----
        if (ftype_in eq 0) then begin ; image
            read_fits_s, map_in, xten
        endif else begin           ;  bin table
            read_fits_s, map_in, prim, xten
        endelse
        ncol = n_tags(xten) - 1
        xhdr = xten.(0)
        hdr = prim.(0)
        xnames = tag_names(xten)

        wmap = strupcase(strtrim(sxpar(hdr,'TELESCOP')))
        junk = where(xnames eq 'PIXEL', nk_pix)
        junk = where(xnames eq 'SIGNAL', nk_sig)
        junk = where(xnames eq 'N_OBS', nk_obs)
        junk = where(xnames eq 'TEMPERATURE', nk_temp)

        wmap_format = 0
        if (((nk_temp+nk_obs) eq 2) && (ncol eq 2)) then begin
            wmap_format = 1
        endif else begin
            if ((nk_pix+nk_sig+nk_obs) gt 0) then begin
               message,'wrong file type : '+map_in
            endif
        endelse
                                ; update header
        sxaddpar,          xhdr,'NSIDE',nside_out
        add_ordering_fits, xhdr, ordering=ord_out
        sxaddpar,          xhdr,'HISTORY','    PROCESSED BY '+routine
        if (nside_out ne nside_in) then sxaddpar,          xhdr,'HISTORY', $
          string('       NSIDE: ',nside_in,' -> ',nside_out,form='(a,i5,a,i5)')
        if (ord_in ne ord_out) then sxaddpar,          xhdr,'HISTORY', $
          '    ORDERING: '+ford_in+' -> '+ford_out
        xtenout=create_struct(xnames[0],xhdr)

        if ((ncol eq 1) || (wmap_format eq 1)) then begin ; Temperature
            map = reorder(xten.(1),in=ord_in,out='NEST')
            map = sub_ud_grade(map, nside_in, nside_out, bad_data=bad_data, pessimistic=pessimistic)
            map = reorder(map,in='NEST',out=ord_out)
                                ; create or upgrade ouput structure
            xtenout=create_struct(xtenout,xnames[1],map)

        endif else begin ; Temperature + Polarization
            index_in  = reorder(lindgen(npix_in), in=ord_in,out='NEST')
            index_out = reorder(lindgen(npix_out),in='NEST',out=ord_out)
            for i=1,ncol do begin
                map = sub_ud_grade((xten.(i))[index_in], nside_in, nside_out, bad_data=bad_data, pessimistic=pessimistic)
;                 map = reorder(map,in='NEST',out=ord_out)
                map = map[index_out]
                xtenout=create_struct(xtenout,xnames[i],map)
            endfor
        endelse

        if (wmap_format eq 1) then begin ; N_obs in WMAP maps
            map = reorder(xten.(2),in=ord_in,out='NEST')
            map = sub_ud_grade(map, nside_in, nside_out, bad_data=bad_data, pessimistic=pessimistic, power=-2)
            map = reorder(map,in='NEST',out=ord_out)
                                ; create or upgrade ouput structure
            xtenout=create_struct(xtenout,xnames[2],map)

            add_nside_fits,    hdr, nside=nside_out
            add_ordering_fits, hdr, ordering=ord_out
            sxaddpar,hdr,'RESOLUTN',nint(alog(nside_out)/alog(2.))
            prim1 = create_struct('HDR',hdr)
            if (n_tags(prim) ge 1) then begin
                for i=1,n_tags(prim)-1 do prim1=create_struct(prim1,(tag_names(prim))[i],prim.(i))
            endif
            prim = prim1
        endif
                                ; write file
        write_fits_sb, map_out, prim, xtenout
    endif else begin
        print,'wrong FITS format in '+map_in
        print,'******** Abort *********'
        return
    endelse


;------------------
; vector to vector
;------------------
endif else begin
    ; find input resolution
    npix = n_elements(map_in[*,0])
    ncol = n_elements(map_in[0,*])
    nside_in = npix2nside(npix,err=error)
    if (error eq 1) then begin
        print,routine+' : the input map is not a Healpix full sky map'
        print,npix,nside_in
        print,'******** ABORT *********'
        return
    endif

    if ((defined(nside_out) || defined(order_out)) && undefined(order_in)) then begin
        print,routine+': the input order (order_in) should be defined for an input vector map'
        print,'******** ABORT *********'
        return
    endif
    if ~(keyword_set(order_out) || arg_present(order_out)) then order_out = order_in
    if undefined(nside_out) then nside_out = nside_in

    ord_in  = decode_ordering(order_in)
    ord_out = decode_ordering(order_out)

    ; same resolution and order : do nothing
    if (nside_in eq nside_out && ord_in eq ord_out) then begin
        map_out = map_in
        return
    endif

    if (ncol eq 1) then begin ; temperature only
        map_out = reorder(map_in,in=ord_in,out='NEST')
        map_out = sub_ud_grade(map_out, nside_in, nside_out, bad_data=bad_data, pessimistic=pessimistic)
        map_out = reorder(map_out,in='NEST',out=ord_out)
    endif else begin
        map_out = reorder(map_in,in=ord_in,out='NEST')
        map_out = sub_ud_grade(map_out, nside_in, nside_out, bad_data=bad_data, pessimistic=pessimistic)
;         if (nside_out gt nside_in) then begin
;             map_out = prograde_qu( map_out, nside_in, nside_out)
;         endif else begin
;             message,'Not available'
;         endelse
        map_out = reorder(map_out,in='NEST',out=ord_out)
    endelse

endelse

return
end

