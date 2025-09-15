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
pro read_fits_cut4, filename, pixel, signal, n_obs, serror, $
HDR = header, XHDR = xheader, NSIDE=nside, ORDERING=ordering, COORDSYS=coordsys,$
 EXTENSION=extension_id, HELP=help
;+
; NAME:
;          read_fits_cut4
;
; PURPOSE:
;         reads a FITS file for incomplete sky coverage with a 4 column format
;        (similar to the one used by DMR) :
;        PIXEL, SIGNAL, N_OBS, SERROR
;
; CATEGORY:
;        fits I/O for Healpix
;
; CALLING SEQUENCE:
;        read_fits_cut4, Filename, Pixel, Signal, [N_obs, Serror, 
;        HDR=, XHDR=, NSIDE=, ORDERING=, COORDSYS=, EXTENSION=, HELP=]
; 
; INPUTS:
;       Filename  : STRING scalar,      
;                   input filename
;
; OPTIONAL INPUTS:
;       EXTENSION= number (0 based) of extension to be read
;
; OUTPUTS:
;       Pixel   : LONG vector
;                   index of Healpix pixel
;
;       Signal   : FLOAT vector (same length as Pixel)
;                  signal in each pixel (eg, average of many measurements)
;
; OPTIONAL OUTPUTS:
;       N_obs    : INT or LONG vector (same length as Pixel)
;                  number of measurements in each pixel
;      
;       Serror  : FLOAT vector (same length as Pixel)
;                 statistical error on Signal in each Pixel
;                 = rms of individual measurement / sqrt(N_obs)
;      
;       HDR= : STRING vector 
;                contains the FITS header of the primary unit
;
;       XHDR= : STRING vector 
;                contains the FITS header of the extension (the most informative)
;
;       NSIDE= Healpix resolution parameter read from FITS file, set to -1 if
;          not found
;
;       ORDERING= pixel ordering, as read from FITS header, either 'RING' or
;          'NESTED' or ' ' (unkwnown)
;
;       COORDSYS= astrophysical coordinate system used, as read from FITS
;           header (value of keywords COORDSYS or SKYCOORD)
;
; KEYWORD PARAMETERS:
;
;       HELP=: if set, an extensive help (this IDL header) is printed
;
; RESTRICTIONS:
;
; PROCEDURE:
;      calls fxbreadm
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;        Sept 2000, EH, Caltech
;        Aug  2002: use HDR and XHDR keywords to match other routines
;        Jan 2005: added NSIDE and ORDERING as outputs
;        Jan 2009: calls init_astrolib
;        Aug  2016: added HELP keyword
;        Oct 2018: more robust column parsing
;-

routine = 'READ_FITS_CUT4'
if (keyword_set(help)) then begin
    doc_library,routine
    return
endif
syntax = 'Syntax : '+routine+', filename, pixel, signal, [n_obs, serror, HDR=, XHDR=, NSIDE=, ORDERING=, COORDSYS=, EXTENSION=, HELP=]'
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

for i=x0,xf do begin

; get basic keywords
    fits_npix =getsize_fits(filename, nmaps=nmaps, obs_npix = obs_npix, type = ftype, ext=i)

    if (ftype ne 3) then begin
        print,routine+': '+filename+' does not have the right format'
        return
    endif

; open extension i+1 and get header
    fxbopen, lun, filename, i+1, xhdr

    tfields=    long(sxpar(xhdr,'TFIELDS'))
    nside  =    long(sxpar(xhdr,'NSIDE',count=count))
    ttypes = strtrim(sxpar(xhdr,'TTYPE*'),2)
    if (count eq 0) then nside = -1
    ordering =  strupcase(strtrim(SXPAR(xhdr,'ORDERING', count=count),2))
    if (count eq 0) then ordering = ' '
    coordsys = sxpar(xhdr,'COORDSYS', count=count)
    if (count eq 0) then coordsys = sxpar(xhdr,'SKYCOORD', count=count)
    if (count eq 0) then coordsys = sxpar(header,'COORDSYS', count=count)
    if (count eq 0) then coordsys = sxpar(header,'SKYCOORD', count=count)
    if (count eq 0) then coordsys = ' ' else coordsys = strtrim(coordsys,2)


; read extension data
    if (obs_npix le 0) then obs_npix = fits_npix
                                ; simultaneous read of 4 columns
    goodnames = ['PIXEL','SIGNAL','N_OBS','SERROR']
    fxbreadm, lun, goodnames, pixel, signal, n_obs, serror,  status = status
    fxbclose, lun
    ncols = round(total(status))
    if (ncols le 1) then begin
        if (~status[0]) then begin ; exit if PIXEL not found
            message,'No valid column found in '+filename
        endif
        message,/info,'WARNING: the columns do not have the expect names in '+filename
        message,/info,'Got:      '+string(ttypes,    form='(10(a,:,", "))')
        message,/info,'Expected: '+string(goodnames, form='(10(a,:,", "))')
        fxbopen, lun, filename, i+1 ; read ignoring column names
        fxbreadm,lun, lindgen(tfields)+1, pixel, signal, n_obs, serror, status = status2
        fxbclose,lun
        ncols2 = round(total(status2))
        if (ncols2 le 1) then begin
            message,/info,'Wrong format for fits file '+filename
            message,'only could read '+strtrim(ncols2,2)+' columns'
        endif
        status = status2
    endif
    if (n_elements(pixel) eq n_elements(signal)) then begin ; standard explicit format
        if (obs_npix lt fits_npix) then begin
            if (status[0]) then pixel  = temporary( pixel[0:obs_npix-1])
            if (status[1]) then signal = temporary(signal[0:obs_npix-1])
            if (status[2]) then n_obs  = temporary( n_obs[0:obs_npix-1])
            if (status[3]) then serror = temporary(serror[0:obs_npix-1])
        endif
    endif else begin            ; compressed explicit format
        message,'inconstitent pixel index and data value in '+filename
    endelse

    xheader = xhdr
endfor

return
end

