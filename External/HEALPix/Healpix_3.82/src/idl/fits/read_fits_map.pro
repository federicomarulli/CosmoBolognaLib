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


PRO READ_FITS_MAP, filename, T_sky, hdr, exthdr, $
                   HELP=help, PIXEL=pixel, SILENT=silent, $
                   NSIDE=nside, ORDERING=ordering, COORDSYS=coordsys, $
                   EXTENSION=extension_id

;+
; NAME:
;       READ_FITS_MAP
;
; PURPOSE:
; 	reads the temperature map and possibly the polarisation map in binary table
; 	extension of FITS file filename
; 	if doesn't find an extension, look in the image part
;
;       T_sky is returned as a 2-dim array 
;       1st dim = number of pixels
;       2nd dim = number of maps present in the file
;
; CALLING SEQUENCE:
; 	READ_FITS_MAP, filename, T_sky, [hdr, exthdr, SILENT=, PIXEL=, NSIDE=,
; 	ORDERING=, COORDSYS=, HELP=]
; 
; INPUTS:
;	filename = String containing the name of the file to be read.
;
; OPTIONAL OUTPUTS:
;	hdr = String array containing the header from the FITS file.
;
;	exthdr = String array containing the header from the extension
;
;       NSIDE = Healpix resolution parameter read from FITS file, set to -1 if
;       not found
;
;       ORDERING = pixel ordering, as read from FITS header, either 'RING' or
;       'NESTED' or ' ' (unkwnown)
;
;       COORDSYS = astrophysical coordinate system used, as read from FITS
;           header (value of keywords COORDSYS or SKYCOORD)
;
;
; OPTIONAL KEYWORDS:
;
;    EXTENSION : extension unit to be read from FITS file: 
;     either its 0-based ID number (ie, 0 for first extension after primary array)
;     or the case-insensitive value of its EXTNAME keyword.
;	If absent, first extension (=0) will be read
;
;    HELP=  : if set, this help header is printed out and the routine exits
;
;    PIXEL      = pixel number to read from or pixel range to read
;                 (in the order of appearance in the file), starting from 0.
;               if >= 0 scalar        : read from pixel to the end of the file
;               if two elements array : reads from pixel[0] to pixel[1] (included)
;               if absent             : read the whole file
;
;    SILENT : if set, no message is issued during normal execution
;

;
; PROCEDURES USED:
;	several (see below)
;
; MODIFICATION HISTORY:
;  October 1997:
;  *strip down edited by E. Hivon from the READ_FITS_MAP.pro 
;  written by A. Banday for COBE data analysis.
;  *added the image reading if no extension, EH		Nov-97
;  *made it faster by using READFITS,/exten		Dec-97
;  addition of /SLOW (=previous implementation using FXBREAD)
;  *made it compatible with multiple entry per row and per column
;  (eg TFORM = '1024E') wich makes the reading faster.  Jan-98
;  *addition of PIXEL                                   Aug-02
;        Jan 2005: added NISDE and ORDERING as outputs
;  May 2005, EH, replaces FINDFILE by FILE_TEST
;  Jan 2008, EH: calls tbfree to remove heap pointer created by TBINFO
;  Jan 2009: calls init_astrolib
;  Feb 2012: calls SELECTREAD which is much faster than FITS_READ
;    and can read Nside=8192 sky maps
;  Jan 2013: added HELP= and EXTENSION= keywords
;  Aug 2017: produces a full-sky map (with empty pixels set to 0) when reading
;           a FITS file in cut-sky format
;  Feb 2020: uses READ_FITS_PARTIAL instead of READ_FITS_CUT4 for partial maps
; requires the THE IDL ASTRONOMY USER'S LIBRARY 
; that can be found at http://idlastro.gsfc.nasa.gov/homepage.html
;
;-

routine = 'read_fits_map'
syntax = 'Syntax : READ_FITS_MAP, filename, T_sky [, hdr, exthdr, PIXEL=, SILENT=, NSIDE=, ORDERING=, COORDSYS=, HELP=, EXTENSION=]'
syntax2 = '   No File Read, Abort !'

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if N_params() LT 2 then begin
    print,syntax
    print,syntax2
    retall
endif

; run astrolib routine to set up non-standard system variables
init_astrolib

; get the primary header information
if (datatype(filename) ne 'STR') then begin
    print,syntax
    print,'Invalid file name'
    message,syntax2
endif

if (not file_test(filename)) then begin
    prim_stc = 0 & xten_stc = 0
    message, 'file '+filename+' not found'
endif

hdr = headfits(filename, errmsg=errmsg)
if (strtrim(errmsg) ne '') then message,errmsg

; checks for the existence of an extension
; if none, reads the image
fits_info,filename, /silent, n_ext=n_ext
if (n_ext eq 0) then begin
    T_sky = readfits(filename,hdr, silent=silent)
    exthdr = ' '
    return
endif

; identify the format in the row (number of entry per row)
;exthdr = headfits(filename, exten=1)
if size(extension_id,/tname) eq 'STRING' then begin
    fits_read, filename, junk, exthdr, /header_only, /no_pdu, extname=extension_id
endif else begin
    xx = defined(extension_id) ? extension_id+1 : 1
    fits_read, filename, junk, exthdr, /header_only, /no_pdu, exten_no=xx
endelse

nside  = LONG(       SXPAR(exthdr,'NSIDE',count=count  ))
if (count eq 0) then nside = -1
ordering =  strupcase(strtrim(SXPAR(exthdr,'ORDERING', count=count),2))
if (count eq 0) then ordering = ' '
coordsys = sxpar(exthdr,'COORDSYS', count=count)
if (count eq 0) then coordsys = sxpar(exthdr,'SKYCOORD', count=count)
if (count eq 0) then coordsys = sxpar(hdr,'COORDSYS', count=count)
if (count eq 0) then coordsys = sxpar(hdr,'SKYCOORD', count=count)
if (count eq 0) then coordsys = ' ' else coordsys = strtrim(coordsys,2)

; actual reading
obsnpix = getsize_fits(exthdr,/header,type=filetype)
if (filetype eq 3) then begin
;     message,/info,'WARNING: input file '+filename+' is in cut-sky format, which can be read properly with READ_FITS_CUT4.'
;     message,/info,'Here a full sky map will be generated, replacing missing data with 0.'
;     read_fits_cut4, filename, pix, c2, c3, c4, extension=extension_id
;     ncol = 1 + defined(c3) + defined(c4)
;     T_sky = make_array(nside2npix(nside), ncol, type=size(c2,/type))
;     T_sky[pix,0] = c2
;     if defined(c3) then T_sky[pix,1] = c3
;     if defined(c4) then T_sky[pix,2] = c4   
    message,/info,'WARNING: input file '+filename+' is in partial sky format, which can be read properly with READ_FITS_PARTIAL.'
    message,/info,'Here a full sky map will be generated, replacing missing data with 0.'
    read_fits_partial, filename, pix, tmp, extension=extension_id
    ; ncol = (size(tmp,/dim))[1]  ; corrected 2022-06-29
    ncol = (size(tmp,/n_dims) gt 1) ? (size(tmp,/dim))[1] : 1
    T_sky = make_array(nside2npix(nside), ncol, type=size(tmp,/type))
    for i=0,ncol-1 do T_sky[pix,i] = tmp[*,i]
endif else begin
    selectread, filename, T_sky, extension=extension_id, /no_pdu
endelse

; find pixels to keep
;max_pix  = n_elements(T_sky[*,0]) - 1L ; costly in memory
;max_pix = (size(T_sky))[1] - 1LL
max_pix = (size(T_sky,/dim))[0] - 1LL
firstpix = 0LL
lastpix = max_pix
cut = 0
if (defined(pixel)) then begin
    if (n_elements(pixel) ge 1 && pixel[0] ge 0) then firstpix  = LONG64(pixel[0])
    if (n_elements(pixel) eq 2) then begin
        if (pixel[1] lt pixel[0]) then begin
            print,pixel
            message,'Invalid pixel range'
        endif
        lastpix   = LONG64(pixel[1])
    endif
    cut = 1
endif
if (cut) then begin
    if (firstpix gt max_pix) then begin
        print,'range : ',firstpix,lastpix,',  available : ',max_pix+1
        message,'Invalid pixel range'
    endif
    if (lastpix gt max_pix) then begin
        print,'Only ',max_pix+1,' pixels available in '+filename+', will truncate the output array'
        lastpix  = lastpix  < max_pix
    endif
    T_sky = T_sky[firstpix:lastpix, *]
endif


; Exit, stage left ....
return
end

