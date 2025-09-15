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
pro selectread_wrapper, file, map, header=header, extension=extension_id, no_pdu=no_pdu, level=level,silent=silent
obsnpix = getsize_fits(file, type=filetype,nside=nside)
if (filetype eq 3) then begin
;     if ~keyword_set(silent) then begin
;         message,/info,level=level,'WARNING: input file '+file+' is in cut-sky format, which can be read properly with READ_FITS_CUT4.'
;         message,/info,level=level,'Here a full sky map will be generated, replacing missing data with 0.'
;     endif
;     read_fits_cut4, file, pix, c2, c3, c4, extension=extension_id, xhdr=header
;     map = make_array(nside2npix(nside), 3, type=size(c2,/type))
;     map[pix,0] = c2
;     map[pix,1] = c3
;     map[pix,2] = c4   
    if ~keyword_set(silent) then begin
        message,/info,level=level,'WARNING: input file '+file+' is in partial sky format, which can be read properly with READ_FITS_PARTIAL.'
        message,/info,level=level,'Here a full sky map will be generated, replacing missing data with 0.'
    endif
    read_fits_partial, file, pix, tmp, extension=extension_id, xhdr=header
    ; ncol = (size(tmp,/dim))[1] ; corrected 2022-06-29
    ncol = (size(tmp,/n_dims) gt 1) ? (size(tmp,/dim))[1] : 1
    map = make_array(nside2npix(nside), ncol, type=size(tmp,/type))
    for i=0,ncol-1 do map[pix,i] = tmp[*,i]
endif else begin
    selectread, file, map, extension=extension_id, no_pdu=no_pdu, header=header
endelse

return
end

PRO read_tqu, fitsfile, TQU, EXTENSION=extension_id, HDR = hdr, XHDR = xhdr, HELP = help, NSIDE=nside, ORDERING=ordering, COORDSYS=coordsys

;+
; NAME:
;    read_tqu
;
; PURPOSE:
;    reads a temperature+polarization Healpix map (T,Q,U)  from a
;    binary table FITS file
;    with optionally the error (dT,dQ,dU) and correlation (dQU, dTU, dTQ)
;    from separate extensions
;
; CALLING SEQUENCE:
;    READ_TQU, fitsfile, TQU, EXTENSION=, HDR = , XHDR = , HELP =, NSIDE=, ORDERING=, COORDSYS=
; 
; INPUTS:
;    fitsfile = String containing the name of the file to be read   
;
; OPTIONAL INPUTS:
;    Extension : extension unit to be read from FITS file: 
;     either its 0-based ID number (ie, 0 for first extension after primary array)
;     or the case-insensitive value of its EXTNAME keyword.
;	If absent, all available extensions are read.
;
; KEYWORDS:
;    Help : if set, provides help on routine
;
; OUTPUTS:
;
;   TQU : array of Healpix maps of size (npix,3,n_ext) where npix is the total
;   number of Healpix pixels on the sky, and (n_ext <=3) the number of
;   extensions read
;     3 maps are present in each extension of the FITS file : 
;      -the temperature+polarization Stokes parameters maps (T,Q,U) in extension 0
;      -the error maps (dT,dQ,dU) in extension 1 (if applicable)
;      -the correlation maps (dQU, dTU, dTQ) in extension 2 (if applicable)
;
;
; OPTIONAL OUTPUT KEYWORDS:
;       HDR      = String array containing the header of the FITS file.
;       XHDR     = String array containing the extension header(s). If 
;                  several extensions are read, then the extension 
;                  headers are returned appended into one string array.
;
;       NSIDE = Healpix resolution parameter read from FITS file, set to -1 if
;          not found
;
;       ORDERING = pixel ordering, as read from FITS header, either 'RING' or
;          'NESTED' or ' ' (unkwnown)
;
;       COORDSYS = astrophysical coordinate system used, as read from FITS
;          header (value of keywords COORDSYS or SKYCOORD)
;
; EXAMPLE:
;
;
; PROCEDURES CALLED:
;       DOC_LIBRARY, FILE_TEST, FITS_INFO, HEADFITS, SELECTREAD
;
; MODIFICATION HISTORY:
;       Dec 2002, E. Hivon (McMurdo)
;       Feb 2005, EH, replaces READ_FITS_S by faster SELECTREAD
;       May 2005, EH, replaces FINDFILE by FILE_TEST
;       Jan 2009: calls init_astrolib
;       Jan 2013: allows Extension to be a string
;       Aug 2017: produces a full-sky map (with empty pixels set to 0) when reading
;           a FITS file in cut-sky format
;       Feb 2020: uses READ_FITS_PARTIAL instead of READ_FITS_CUT4 for partial maps
;
; requires the THE IDL ASTRONOMY USER'S LIBRARY 
; that can be found at http://idlastro.gsfc.nasa.gov/homepage.html
;
;-

if (keyword_set(help)) then begin
    doc_library,'read_tqu'
    return
endif

syntax = 'SYNTAX: READ_TQU, fitsfile, TQU, [Extension=, Hdr=, Xhdr=, Help=, NSIDE=, ORDERING=, COORDSYS=]'
if n_params() eq 0 then begin
    print,syntax,form='(a)'
    return
end
if n_params() ne 2 then begin
    print,syntax,form='(a)'
    message,'Not enough argument, no file read'
endif

if (not file_test(fitsfile)) then message,'file '+fitsfile+' not found'

; run astrolib routine to set up non-standard system variables
init_astrolib

fits_info, fitsfile, /silent, n_ext = n_ext, extname=extnames

xtname = ''
if undefined(extension_id) then begin
    x0 = 0
    xf = n_ext - 1
endif else begin
    if size(extension_id,/tname) eq 'STRING' then begin
        xtname = extension_id[0]
        x0 = 0
        xf = x0
        n_ext = 1
    endif else begin
        if (extension_id + 1) gt n_ext then begin
            print,' Requested extension ',extension_id,' (0 based) from ',fitsfile
            print,' Only found ',n_ext,' extensions.'
            message,' Abort'
        endif else begin
            x0 = extension_id[0]
            xf = x0
            n_ext = 1
        endelse
    endelse
endelse
hdr  = HEADFITS(fitsfile)
;print,' Reading ',n_ext,' extension(s) from ',fitsfile

; -------- extension -----------

; simply read the extensions from the FITS file
savehdr = ''
for i = x0,xf do begin

    exten = (xtname ne '' && n_ext eq 1) ? xtname : i
    ; read data
    if (n_ext eq 1) then begin
        ; if only 1 extension: read directly into TQU
        selectread_wrapper, fitsfile, tqu, header = xhdr, exten=exten, /no_pdu,level=-1,silent=0
        nmaps = n_elements(tqu[0,*])
    endif else begin
        ; if several extensions, read into tmp array, TO BE IMPROVED
        selectread_wrapper, fitsfile, tmp, header = xhdr, exten=exten, /no_pdu,level=-1,silent=(i ne x0)
        npix  = n_elements(tmp[*,0])
        nmaps = n_elements(tmp[0,*])
        if (i eq 0) then begin
            tqu = fltarr(npix,nmaps,n_ext)
        endif else begin
            if (npix ne npix_old or nmaps ne nmaps_old) then begin
;                message,'Extensions are not consistent in '+fitsfile
                print,' ERROR: '
                message,/info,fitsfile+' contains '+strtrim(n_ext,2)+' extensions which can not be merged into a single array because of their different format'
                message,/info,'Specify the extension you want to read by providing its 0-based index or its name. Available choices are: '
                for k=1,(n_ext<5) do begin
                    print, k-1, "  '"+strtrim(extnames[k],2)+"'",form="(i5,a)"
                endfor
                if (n_ext gt 5) then begin
                    print, '  ...    ....'
                    k = n_ext
                    print, k-1,"  '"+strtrim(extnames[k],2)+"'",form="(i5,a)"
                endif
                print,syntax,form='(a)'
                message,'Aborting.'
            endif
        endelse
        tqu[0,0,i] = tmp
        npix_old = npix
        nmaps_old = nmaps
    endelse

    if (nmaps ne 3) then begin
        print,' WARNING: '+strtrim(nmaps,2)+' maps available in extension '+ $
              strtrim(exten,2)+" ('"+strtrim(extnames[exten+1],2)+"')"+' of '+fitsfile
        print,'      Expected 3'
    endif

    ; deal with header
    hdr_merge = [hdr, xhdr] ; to be compatible with WMAP
    if (i eq x0) then begin
        nside  = LONG(       SXPAR(hdr_merge,'NSIDE',count=count  ))
        if (count eq 0) then nside = -1
        ordering =  strupcase(strtrim(SXPAR(hdr_merge,'ORDERING', count=count),2))
        if (count eq 0) then ordering = ' '
        coordsys = sxpar(hdr_merge,'COORDSYS', count=count)
        if (count eq 0) then coordsys = sxpar(hdr_merge,'SKYCOORD', count=count)
        if (count eq 0) then coordsys = sxpar(hdr,'COORDSYS', count=count)
        if (count eq 0) then coordsys = sxpar(hdr,'SKYCOORD', count=count)
        if (count eq 0) then coordsys = ' ' else coordsys = strtrim(coordsys,2)
    endif

    savehdr = [savehdr,xhdr]

endfor

xhdr = savehdr

; Exit routine
Exit:
return
end



