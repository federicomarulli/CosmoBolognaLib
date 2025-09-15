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
pro read_fits_partial, filename, pixel, iqu, $
HDR = header, XHDR = xheader, NSIDE=nside, ORDERING=ordering, COORDSYS=coordsys,$
 EXTENSION=extension_id, UNITS=units, HELP=help
;+
; NAME:
;          read_fits_partial
;
; PURPOSE:
;         reads a FITS file for incomplete (polarized) sky coverage
;        PIXEL, TEMPERATURE, Q_POLARISATION, U_POLARISATION
;
; CATEGORY:
;        fits I/O for Healpix
;
; CALLING SEQUENCE:
;        read_fits_partial, Filename, Pixel, IQU, 
;        [COORDSYS=, EXTENSION=, HDR=, HELP=, NSIDE=, ORDERING=, UNITS=, XHDR=]
; 
; INPUTS:
;       Filename  : STRING scalar,      
;                   input filename
;
; OPTIONAL INPUTS:
;       EXTENSION= number (0 based) of extension to be read
;
; OUTPUTS:
;       Pixel   : INT, LONG or LONG64 vector of length np
;                   index of Healpix pixel
;
;       IQU   : FLOAT or DOUBLE array of size (np, nc)
;                  I,Q and U Stokes parameters for each pixel (if nc=3)
;                 or I of each pixel (if nc=1)
;
; OPTIONAL OUTPUTS:
;      
;       COORDSYS= astrophysical coordinate system used, as read from FITS
;           header (value of keywords COORDSYS or SKYCOORD)
;
;       HDR= : STRING vector 
;                contains the FITS header of the primary unit
;
;       NSIDE= Healpix resolution parameter read from FITS file, set to -1 if
;          not found
;
;       ORDERING= pixel ordering, as read from FITS header, either 'RING' or
;          'NESTED' or ' ' (unkwnown)
;
;
;       UNITS = physical units of each column of the table (except PIXEL)
;
;       XHDR= : STRING vector 
;                contains the FITS header of the extension (the most informative)
;
; KEYWORD PARAMETERS:
;
;       HELP=: if set, an extensive help (this IDL header) is printed
;
; RESTRICTIONS:
;
; PROCEDURE:
;      calls fits_read
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;     Jan 2020: EH, IAP, adapted from write_fits_cut4
;-

routine = 'READ_FITS_PARTIAL'
if (keyword_set(help)) then begin
    doc_library,routine
    return
endif
syntax = 'Syntax : '+routine+', filename, pixel, IQU, [HDR=, XHDR=, COORDSYS=, EXTENSION=, HELP=, NSIDE=, ORDERING=, UNITS=]'
if N_params() LT 3 then begin
      print,syntax
      return
endif
if (datatype(filename) ne 'STR') then begin
      print, syntax
      message,'abort'
endif

; run astrolib routine to set up non-standard system variables
init_astrolib

fits_info, filename, /silent, n_ext = n_ext
x0 = (defined(extension_id)) ? extension_id : 0
xf = x0 & n_ext = 1

; get the primary header information
header = headfits(filename)

fits_open, filename, fcb
for i=x0,xf do begin

; get basic keywords
    fits_npix =getsize_fits(filename, nmaps=nmaps, obs_npix = obs_npix, type = ftype, ext=i)

    if (ftype ne 3) then begin
        print,routine+': '+filename+' does not have the right format'
        return
    endif

; open extension i+1 and get header
    xtnum = i+1
    fits_read, fcb, void, xhdr, exten_no = xtnum, /header_only, no_pdu = 1

    tfields=    long(sxpar(xhdr,'TFIELDS'))
    nside  =    long(sxpar(xhdr,'NSIDE',count=count))
    if (count eq 0) then nside = -1
    ttypes = strtrim(sxpar(xhdr,'TTYPE*'),2)
    tunits = strtrim(sxpar(xhdr,'TUNIT*', count=n_units),2)
    tforms = strtrim(sxpar(xhdr,'TFORM*'),2)
    ordering =  strupcase(strtrim(SXPAR(xhdr,'ORDERING', count=count),2))
    if (count eq 0) then ordering = ' '
    coordsys = sxpar(xhdr,'COORDSYS', count=count)
    if (count eq 0) then coordsys = sxpar(xhdr,'SKYCOORD', count=count)
    if (count eq 0) then coordsys = sxpar(header,'COORDSYS', count=count)
    if (count eq 0) then coordsys = sxpar(header,'SKYCOORD', count=count)
    if (count eq 0) then coordsys = ' ' else coordsys = strtrim(coordsys,2)

    if (obs_npix le 0) then obs_npix = fits_npix

    if (ttypes[0] ne 'PIXEL') then begin
        message,/info,'First column of '+filename+' is '+ttypes[0]+'  instead of the expected PIXEL.'
        message,/info,'Proceed at your own risk'
    endif

; find number of words in extension
    n_wpr  = (fcb.axis)[0,xtnum] ; words per row in extension
    n_rows = long64((fcb.axis)[1,xtnum]) ; number of rows
    n_words = n_wpr * n_rows

    ; list of non-PIXEL columns
    nmaps = tfields - 1
    cols = indgen(nmaps) + 2

    tbinfo, xhdr, tab_xhdr
    nentry = max(tab_xhdr.numval[cols-1])
    npix = nentry * n_rows

    tpix  = (tab_xhdr.idltype)[0]
    types = (tab_xhdr.idltype)[cols-1]
    type = max(types) ; find out type of output array to fit all data except PIXEL
    junk = where(types eq 7, n_string)
    if (n_string ge 1) then begin
        message,'Table in input file contains strings'
    endif

; create array receiving final data
    pixel = make_array(type=tpix ,npix, /nozero)
    iqu = make_array(type=type, npix, nmaps, /nozero)

; read data piece by piece and process each piece individually
    stride = 5.e6 > n_wpr       ; 5 MB or 1 row per piece (June-2008)
;stride = 500.e6 > n_wpr ; 500 MB or 1 row per piece (Dec-2011)
    stride = FLOOR(stride / n_wpr) * n_wpr
    w_start = long64(0)
    pstart = long64(0)
    while (w_start LE (n_words-1) ) do begin
                                ; read one piece
        w_end = (w_start + stride - 1L) < (n_words-1)
        fits_read, fcb, data, exten_no = xtnum, first=w_start, last=w_end
        nr = (w_end - w_start + 1) / n_wpr ; number of rows read
        np = nr * nentry        ; number of pixels read
        data = reform(data, n_wpr, nr, /overwrite) ; required by tbget
        pixel[pstart] = (tbget(tab_xhdr, data, 1))[*]
        for i=0,nmaps-1 do begin
            iqu[pstart, i] = (tbget(tab_xhdr, data, cols[i]))[*]
        endfor
                                ; get ready for next piece
        w_start = w_end + 1
        pstart = pstart + np
    endwhile
    data = 0
    tbfree, tab_xhdr

    ; deal with units
    if (arg_present(units)) then begin
        units = strarr(nmaps) ; '' as default value
        if (n_units gt 0 && n_elements(tunits) gt 1) then begin
            units = tunits[1:*] ; drop out units of PIXEL column
        endif
    endif

;     fxbread, lun, pixel, 1
;     nc = tfields-1
;     if (nc gt 0) then begin
;         col = lindgen(nc) + 2
;         junk = where( tforms[1:nc-1] ne tforms[1], nbad)
;         if (nbad ne 0) then begin
;             message,/info,'Mismatched column types in '+filename
;             message,/info,'Proceed at your own risk'
;         endif
;         fxbread, lun, iqu, col
;     endif

    xheader = xhdr
endfor
fits_close, fcb

return
end

