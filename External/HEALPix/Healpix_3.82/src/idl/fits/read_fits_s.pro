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
pro read_fits_s, filename, prim_stc, xten_stc, merge=merge, extension = extension_id, columns=columns, help=help
;+
; NAME:
;       READ_FITS_S
;
; PURPOSE:
; 	reads a multi column binary or ascii table extension (eg temperature
; 	and polarisation or polarisation power spectrum) of FITS file filename
; 	if doesn't find an extension, look in the image part
;
; CALLING SEQUENCE:
; 	READ_FITS_S, Filename, Prim_Stc [, Xten_Stc, Merge=, Extension=, Columns=]
; 
; INPUTS:
;	filename = String containing the name of the file to be read.
;
; OUTPUTS:
;	Prim_Stc = structure containing
;                  - the primary header, tag : 0, tag name : HDR
;                  - the primary image (if any) as an array,
;                    tag : 1, tag name : IMG, 
;
;       if MERGE is set Prim_Stc contains
;                  - the concatenated primary and extension header : HDR
;                  - the primary image if any, tag name : IMG
;                  - the data column, tag name : name given in TTYPEi (with all
;                    spaces removed and only letters, digits and underscore)
;       and xten_stc is 0
;                    
; OPTIONAL OUTPUTS:
;       if MERGE is not set
;	xten_stc = structure containing
;                  - the extension header, tag 0, tag name HDR
;                  - the data column
;                    tag : i, tag name : name given in TTYPEi (with all
;                    spaces removed and only letters, digits and underscore)
;
; OPTIONAL INPUT KEYWORDS:
;      Merge = if set, the content of the primary and secondary units
;      are merged
;
;      Extension : extension unit to be read from FITS file: 
;        either its 0-based ID number (ie, 0 for first extension after primary array)
;        or the case-insensitive value of its EXTNAME keyword.
;	 If absent, first extension (=0) will be read
;
;      Columns = list of columns to be read from a binary table 
;        can be a list of integer (1 based) indexing the columns positions
;        or a list of names matching the TTYPE* of the columns
;        by default, all columns are read
;
; PROCEDURES USED:
;	mrdfits, fits_read (astron library)
;
; TIPS:
;    to plot column 5 vs. column 1 of 'file' irrespective of their name
;        read_fits_s, 'file', prims, exts
;        x = exts.(1)
;        y = exts.(5)
;        plot,x,y
;
;    to plot the column SIGNAL versus the column NOISE without knowing
;    their position
;        read_fits_s, 'file', prims, exts
;        x = exts.signal
;        y = exts.noise
;        plot,x,y
;
; MODIFICATION HISTORY:
;  March 1999, EH, version 1.0
;  March 21, modified to deal with large files
;  Feb 17, 2003, upgraded to deal faster with WMAP map format
;  Oct 2004, cosmetic changes
;       May 2005, EH, replaces FINDFILE by FILE_TEST
;  Jan 2007: will exit gracefully if file contains image, but not requested extension
;  Jan 2008, EH: calls tbfree to remove heap pointer created by TBINFO
;  June 2008, EH: can deal with file with large TFORM
;  Jan 2009: calls init_astrolib
;  Jan 2013: allows EXTENSION to be a string
;  Mars 2013: added HELP
;
; requires the THE IDL ASTRONOMY USER'S LIBRARY 
; that can be found at http://idlastro.gsfc.nasa.gov/homepage.html
;
;-

code = 'READ_FITS_S'
syntax =  'Syntax : '+code+', Filename, Prim_Stc, Xten_Stc, [Merge=, Extension=, Columns=, Help=]'

if keyword_set(help) then begin
    doc_library,code
    return
endif

if N_params() eq 0 then begin
    print, syntax
    retall
endif

; run astrolib routine to set up non-standard system variables
init_astrolib

if N_params() LT 1 or N_params() GT 3 then begin
    prim_stc = 0 & xten_stc = 0
    message,syntax
endif

if (not file_test(filename)) then begin
    prim_stc = 0 & xten_stc = 0
    message, 'file '+filename+' not found'
endif

merge = keyword_set(merge)
expect_extension = n_params() ge 3 || merge

; ------ primary unit : header and extension ------

image = MRDFITS(filename,0,hdr,/silent, status=status)
image_found = (status eq 0)
prim_stc = CREATE_STRUCT('HDR',hdr)
if (image_found) then prim_stc = CREATE_STRUCT(prim_stc,'IMG',image)

xten_stc = 0
fits_info, filename, /silent, n_ext=n_ext, extname=extnames
if size(extension_id,/TNAME) eq 'STRING' then begin
    extension_idp1 = extension_id[0]
    junk = where(strupcase(strtrim(extnames,2)) eq strupcase(extension_idp1), count)
    not_found = (count eq 0)
endif else begin
    extension_idp1 = defined(extension_id) ? extension_id[0] + 1 : 1
    not_found = (extension_idp1 gt n_ext)
endelse
if (not_found) then begin ; no extension available
    if (image_found) then begin
        if (expect_extension) then print,'WARNING: Found 1 image but not the requested extension in '+filename
        return
    endif else begin
        ; no image, no extension: return in Error
        message,' Can not access requested extension in '+filename
    endelse
endif

; ----- if there is the required extension ------
table = MRDFITS(filename,extension_idp1,xthdr,range=1,/silent, columns=columns)  ; first row
tags = TAG_NAMES(table)
n_tag = N_TAGS(table)
bitpix  = ABS(ROUND(SXPAR(xthdr,'BITPIX'))) ; bits per 'word'
n_wpr   =     ROUND(SXPAR(xthdr,'NAXIS1')) ; 'word' per row
n_rows  =     ROUND(SXPAR(xthdr,'NAXIS2')) ; number of rows
byt_row = bitpix * n_wpr / 8 ; bytes per row
;;;;;stride = ((1024L^2 * 20L) / byt_row) > 1 ; strides of 20MB or 1 row
stride = ((1024L^2 * 20L) / byt_row) > 1 ; strides of 20MB or 1 row
ishift = 0
n_entry = intarr(n_tag)

; start building the structure with the header
case (merge) of
    0 : xten_stc = CREATE_STRUCT('HDR',xthdr)
    1 : begin
        prim_stc = CREATE_STRUCT('HDR',[hdr,xthdr])
        if (image NE 0) then begin
            prim_stc = CREATE_STRUCT(prim_stc,'IMG',image)
            ishift = 1
        endif
    end
endcase

; mrdfits is slower than fits_read + tbget
size = (bitpix/8.) * float(n_wpr) * float(n_rows) / 1024.^2 ; size in MB
print,'File size (MB) ',size
if (size le 600.) then begin ; smaller than 600MB, read in one sitting
    if defined(columns) then begin
        table = MRDFITS(filename,extension_idp1,/SILENT, columns=columns)
        for i=0L, n_tag-1 do begin
            nn = n_elements(table.(i))
            ni = n_entry(i)
            ;;help,prim_stc.(i+ishift+1)(r_start*ni:r_end*ni+ni-1),reform(table.(i),nn)
            case (merge) of
                0:xten_stc = create_struct(xten_stc, tags[i], (table.(i)[*]))
                1:prim_stc = create_struct(prim_stc, tags[i], (table.(i)[*]))
            endcase
        endfor
    endif else begin
        if (size(extension_idp1,/tname) eq 'STRING') then begin
            fits_read,filename, data, header_fr, extname=extension_idp1, /no_pdu
        endif else begin
            fits_read,filename, data, header_fr, exten_no=extension_idp1, /no_pdu
        endelse
        tbinfo, header_fr, tb_str
        for i=0L, n_tag-1 do begin
            case (merge) of
                0:xten_stc = create_struct(xten_stc, tags[i], (tbget(tb_str,data,i+1))[*])
                1:prim_stc = create_struct(prim_stc, tags[i], (tbget(tb_str,data,i+1))[*])
            endcase
        endfor
        tbfree, tb_str
    endelse
endif else begin
; build up the final structure according to data 
    for i=0L, n_tag-1 do begin      
        type = datatype(table.(i),2)
        n_entry(i) = n_elements(table.(i))
        n_el = n_entry(i)*n_rows
        case (merge) of
            0 : xten_stc = CREATE_STRUCT(xten_stc, tags(i), MAKE_ARRAY(n_el,TYPE=type,/nozero))
            1 : prim_stc = CREATE_STRUCT(prim_stc, tags(i), MAKE_ARRAY(n_el,TYPE=type,/nozero))
        endcase
    endfor

; read the data by piece to avoid overloading memory and fill in the structure
    r_start = 0L
    while (r_start LE (n_rows-1) ) do begin
        r_end = (r_start + stride - 1L) < (n_rows-1)
        table = MRDFITS(filename,extension_idp1,range=[r_start,r_end],/SILENT, columns=columns)
        for i=0L, n_tag-1 do begin
            nn = n_elements(table.(i))
            ni = n_entry(i)
            ;;help,prim_stc.(i+ishift+1)(r_start*ni:r_end*ni+ni-1),reform(table.(i),nn)
            case (merge) of
                0 : xten_stc.(i+ishift+1)[r_start*ni:r_end*ni+ni-1] = (table.(i))[*]
                1 : prim_stc.(i+ishift+1)[r_start*ni:r_end*ni+ni-1] = (table.(i))[*]
            endcase
        endfor
        r_start = r_end + 1L
    endwhile

endelse




return
end

