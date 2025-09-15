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
pro change_polcconv, file_in, file_out, i2c=i2c, c2i=c2i, c2c=c2c, i2i=i2i, force=force
;+
; NAME:
;          change_polcconv
;
; PURPOSE:
;          change the coordinate convention in FITS file containing a polarised sky map.
;          The main effect is to change the sign of the U Stokes parameter,
;          and add/update the POLCCONV FITS header with either
;          COSMO or IAU value
;
; CATEGORY:
;
; CALLING SEQUENCE:
;          CHANGE_POLCCONV, File_In , File_Out, [/I2C, /C2I, /C2C, /I2I, /FORCE]
;
; INPUTS:
;           File_In: Input FITS file
;
;           File_Out: Output FITS file, must be different from File_In
;
; KEYWORD PARAMETERS:
;          one and only one of the following four should be set
;
;          I2C: changes from IAU to COSMO coordinate convention
;              -if POLCCONV is not found or found with value 'IAU', 
;               or FORCE is set, 
;               it is added/replaced with value 'COSMO', 
;               and the sign of the U stokes parameter map is changed;
;              -if POLCCONV already has value 'COSMO', 
;               and FORCE is NOT set,
;               File_In is copied unchanged into File_Out,
;
;          C2I: changes from COSMO to IAU coordinate convention
;              -if POLCCONV is not found or found with value 'COSMO', 
;               or FORCE is set, 
;               it is added/replaced with value 'IAU', 
;               and the sign of the U stokes parameter map is changed;
;              -if POLCCONV already has value 'IAU',  
;               and FORCE is NOT set,
;               File_In is copied unchanged into File_Out.
;
;          C2C: does NOT change coordinate system
;              -if POLCCONV is found with value 'IAU', and FORCE is NOT set,
;               program will issue error message and no file is written;
;              -in all other cases POLCCONV is set/reset to 'COSMO' in File_Out, 
;               but data is NOT changed
;
;          I2I: does NOT change coordinate system
;              -if POLCCONV is found with value 'COSMO', and FORCE is NOT set,
;               program will issue error message and no file is written,
;              -in all other cases POLCCONV is set/reset to 'IAU', but
;               data is NOT changed
;
;          FORCE: if set, the value of POLCCONV read from the FITS header is
;          ignored.
;              The sign of U is swapped (if used with C2I or I2C), and/or the
;              keyword is updated to IAU (if used with I2I or C2I)
;              or to COSMO (if used with C2C or I2C).
;
;
;
; OUTPUTS:
;          none
;
; OPTIONAL OUTPUTS:
;          none
;
; COMMON BLOCKS:
;          none
;
; SIDE EFFECTS:
;          write a new FITS file
;
; RESTRICTIONS:
;          only works with recognised FITS formats
;
; PROCEDURE:
;          read file, modify header, write file
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;          v1.0, May 2005, Eric Hivon
;          v1.1, Oct 2006, EH: added /force
;          v2.0, Oct 2017, EH: supports Planck formats
;          v2.1, Feb 2020, EH: deal with cut sky file containing IQU in 1st extension
;-

t0 = systime(1)
syntax = 'CHANGE_POLCCONV, File_In , File_Out [, /I2C, /C2I, /C2C, /I2I, /FORCE]'

if (n_params() ne 2) then begin
    print,syntax
    return
endif

do_i2c = keyword_set(i2c) 
do_c2i = keyword_set(c2i)
do_c2c = keyword_set(c2c)
do_i2i = keyword_set(i2i)
do_force = keyword_set(force)
sum = do_i2c + do_c2i + do_c2c + do_i2i
if (sum ne 1) then begin
    print,syntax
    message,'One and only one keyword out of /I2C, /C2I, /C2C, /I2I should be set'
endif
if do_c2c then edit = 'C2C'
if do_c2i then edit = 'C2I'
if do_i2c then edit = 'I2C'
if do_i2i then edit = 'I2I'


if ~file_test(file_in) then begin
    message,'File '+file_in+' not found'
endif

if file_same(file_in, file_out) then begin
    print, syntax
    message,'For safety reasons, the input and output files must be different.'
endif

; full sky TQU format (with extra extension for correlation)
; WMAP (I,Q,U) and (I,Q,U,S) formats (extra column and/or extra extension with correlation)
; Planck polarized maps (with extra columns for correlation)
; cut sky format (multi extension ?)

kw_pol = 'POLCCONV'
kw_iau = 'IAU'
kw_cos = 'COSMO'
fits_info, File_In, /silent, n_ext = n_ext
junk = getsize_fits(File_In, type = type)

print,file_in+' -> ['+edit+'] -> '+file_out

if ~ (type eq 2 || type eq 3) then begin
    message,/info,'Input file is not a HEALPix map'
    message,/info,'File will be copied unchanged'
    file_copy, file_in, file_out, /overwrite
    return
endif

my_type = type
if (type eq 3) then begin
    xhdr  = headfits(File_In,ext=1)
    names = sxpar(xhdr,'TTYPE*')
    if ( array_equal(strupcase(strmid(names[0:3],0,1)),['P','S','N','S']) ) then my_type = 30
    if ( array_equal(strupcase(strmid(names[0:3],0,1)),['P','T','Q','U']) ) then my_type = 31
    if (my_type lt 30) then begin
        message,/info,'Input file could not be properly identified'
        message,/info,'File will be copied unchanged'
        file_copy, file_in, file_out, /overwrite
        return
    endif
endif

; if (n_ext ne 3 and type eq 3) then begin
;     message,'Unpolarised format cut sky format. No changes made',/info
;     return
; endif

modify_map = (do_c2i || do_i2c)

for i=0,n_ext-1 do begin
    ; read data
    case my_type of
        30: read_fits_cut4,    File_In, pixel, signal, n_obs, serror, hdr=hdr, xhdr=xhdr, ext = i
        31: read_fits_partial, File_In, pixel, iqu,                   hdr=hdr, xhdr=xhdr, ext = i
        else:begin
            read_fits_s, File_In, prim, xtn, ext=i
            hdr = prim.hdr  & xhdr = xtn.hdr
        end
    endcase
    names   = sxpar(xhdr,'TTYPE*')
    extname = sxpar(xhdr,'EXTNAME')
    nside   = long(sxpar(xhdr,'NSIDE'))
    xhdr_tmp = [hdr,xhdr]
    polcconv = strupcase(sxpar(xhdr_tmp,KW_POL,count=count))
    in_xxx = 1 & in_iau = 0 & in_cos = 0
    str_ext = 'extension #'+strtrim(i,2)+' ('+strtrim(extname,2)+')'
    ;nmaps = n_elements(tqu[0,*])
    nmaps = n_elements(names)
    xxQU= (nmaps ge 4) ? product(strcmp(names[2:3],['Q','U'],1,/fold_case),/preserve):0 ; first 4 column names are *,*,Q*,U*
    xQU = (nmaps ge 3) ? product(strcmp(names[1:2],['Q','U'],1,/fold_case),/preserve):0 ; first 3 column names are *,Q*,U*
     QU = (nmaps eq 2) ? product(strcmp(names[0:1],['Q','U'],1,/fold_case),/preserve):0 ; first 2 column names are Q*,U*
    
    case count of
        0: in_xxx = 1
        1: begin
            in_iau = (strtrim(polcconv) eq KW_IAU)
            in_cos = (strtrim(polcconv) eq KW_COS)
            in_xxx = ~ (in_iau || in_cos)
        end
        else: begin
            print,'FITS file: ',strtrim(File_In)
            message,'several '+KW_POL+' keywords in header of '+str_ext
        end
    endcase
    done = 0
    non_pol = 0 ; assume polarization
    do_edit = 0
    head_only = 0
    if (~ do_force) then begin
        if ((in_cos && do_i2i) || (in_iau && do_c2c)) then begin
            print,'in FITS file: ',strtrim(File_In)
            print,KW_POL+'= ',polcconv
            print,'can not perform requested '+edit+' edition'
            message,'Use the /FORCE keyword to only update '+KW_POL+' value.'
        endif
        done = ( (in_cos && (do_i2c || do_c2c)) || (in_iau && (do_c2i || do_i2i)) )
    endif

    if (modify_map && ~done) then begin
        ; change sign of U and related quantities
        case my_type of
            30: begin ; cut sky format (pixel, T, hT, sT; pixel, Q, hQ, sQ; pixel, U, hU, sU)
                if (i eq 2) then begin
                    clist   = [1] ;change sign of U (2nd column of last extension)
                    do_edit = 1
                endif else begin
                    head_only = (n_ext eq 3) ; update header of 1st and 2nd extensions
                endelse
            end
            31: begin
                if (xxQU) then begin
                    clist   = [3] ;change sign of U (4th column)
                    do_edit = 1
                endif
            end
            else: begin
                clist = [-1]

                if (nmaps le 1) then non_pol = 1
                
                if (nmaps eq 2) then begin
                    if (n_ext eq 1 && QU) then begin
                                ; Planck 2 columns (Q*, U*)
                        clist   = [1] ; change sign of U
                        do_edit = 1
                    endif else begin
                        non_pol = 1
                    endelse
                endif
            
                if (nmaps eq 3) then begin
                                ; standard Healpix full sky format (T, Q, U;   TT, QQ, UU;  QU, TU, TQ)
                                ; should NOT affect unpolarized Planck maps (I_Stokes ,HITS, II_COV) 
                    if (i eq 0) then begin
                        if (xQU) then begin
                            clist   = [2] ; change sign of U
                            do_edit = 1
                        endif else begin
                            non_pol = 1
                        endelse
                    endif
                    if (i eq 2) then begin
                        clist   = [0, 1] ; change sign of QU and of TU
                        do_edit = 1
                    endif
                endif
                
                if (n_ext ge 1 && nmaps eq 5) then begin
                                ; Planck 5 cols (I_STOKES, Q_STOKES, U_STOKES, TMASK, PMASK) + BEAM_TF extension
                    if (i eq 0 && xQU) then begin
                        clist   = [2] ; change sign of U
                        do_edit = 1
                    endif
                endif
                
                if (n_ext eq 1 && nmaps eq 10) then begin
                                ; Planck R3 (2017) 10 cols (I_STOKES, Q_STOKES, U_STOKES, HIT, II_COV, IQ_COV, IU_COV, QQ_COV, QU_COV, UU_COV)
                    if (i eq 0 && xQU) then begin
                        clist = [2, 6, 8] ; change sign of U, IU, QU
                        do_edit = 1
                    endif
                endif


                if (nside le 512 && n_ext eq 2) then begin
                                ; WMAP: https://lambda.gsfc.nasa.gov/product/map/dr5/skymap_file_format_info.cfm
                    if (nmaps ge 5) then begin
                                ; WMAP (I,Q,U,S) format (T, Q, U, S, N_OBS ; 
                                ;      N_OBS, M11=SS, M12=SQ, M13=SU, M23=QU, M22=QQ, M33=UU)
                        if (i eq 0 && nmaps eq 5  && xQU) then begin
                            clist   = [2] ; change sign of U
                            do_edit = 1
                        endif
                        if (i eq 1 && nmaps eq 7) then begin
                            clist   = [3,4] ; change sign of SU and QU
                            do_edit = 1
                        endif
                    endif
                
                    if (nmaps eq 4) then begin 
                                ; WMAP (I,Q,U) format (T, Q, U, N_OBS ; N_OBS, QQ, QU, UU)
                        clist   = [2] ; change sign of U and of N_QU
                        do_edit = 1
                    endif
                endif
                
                if (~do_edit) then begin
                    if (non_pol) then begin
                        print,' Unpolarised format, no change in '+str_ext
                    endif else begin
                        print,'WARNING: unrecognised format in '+str_ext
                        print,nmaps,n_ext
                        print,names
                        done = 1
                    endelse
                endif
            
            end
        endcase

        if (do_edit) then begin
            sxaddpar, xhdr, 'HISTORY',strtrim(file_in)
            sxaddpar, xhdr, 'HISTORY','edited on '+today_fits(/utc,/time)+' UTC'
            for kc=0,n_elements(clist)-1 do begin
                column = clist[kc]
                case my_type of
                    30:   signal         *= -1
                    31:   iqu[*,column-1] *= -1
                    else: xtn.(column+1) *= -1
                endcase
                str_col = 'column #'+strtrim(column+1,2)+' ('+strtrim(names[column],2)+')'
                print,' Flip sign of '+str_col+' in '+str_ext
                sxaddpar, xhdr, 'HISTORY','modified '+str_col+' to match coord. conv. ('+kw_pol+')'
            endfor
        endif
    endif
    if ( (do_edit || in_xxx || do_force || head_only) && ~non_pol) then begin
        new_cconv = (do_c2c || do_i2c) ? KW_COS : KW_IAU
        sxaddpar, xhdr, KW_POL, new_cconv,' Coord. convention for polarisation (COSMO/IAU)'
        print, (in_xxx ? ' Add ' : ' Update ') + 'keyword '+strtrim(KW_POL)+'='+new_cconv+' in header of '+str_ext
    endif else begin
        print, ' No change in header of '+str_ext
    endelse

    ; write new data
    case my_type of
        30: write_fits_cut4,    File_Out, pixel, signal, n_obs, serror, hdr=hdr, xhdr=xhdr, /pol, ext=i
        31: write_fits_partial, File_Out, pixel, iqu,                   hdr=hdr, xhdr=xhdr,       ext=i
        else: begin
            ; if (nmaps eq 1) then tqu = reform(tqu, n_elements(tqu), 1, /overwrite)
            ; write_tqu,      File_Out, tqu, hdr=hdr, xhdr=xhdr, ext = i
            write_fits_sb, File_Out, prim, xtn, ext = i, xthdr=xhdr, /nothealpix
        end
    endcase
endfor

spawn,/sh,'ls -l '+strtrim(file_in,2)+' '+strtrim(file_out,2)
t1 = systime(1)
print, "Time [s] = ",t1-t0
print,'     ----------------        '

return
end
