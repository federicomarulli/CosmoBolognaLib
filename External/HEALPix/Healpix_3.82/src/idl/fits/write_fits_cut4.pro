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
pro write_fits_cut4, filename, pixel, signal, n_obs, serror, hdr=prim_header, xhdr = info_header, $
Coordsys=coordsys, Nested=nested, Ring=ring, Ordering=ordering, Nside = nside_usr, $
Units = units, Extension = extension_id, polarisation=polar_in, Help=help
;+
; NAME:
;         write_fits_cut4
;
; PURPOSE:
;        writes a FITS file for incomplete sky coverage with a 4 column format
;        (similar to the one used by DMR) :
;        PIXEL, SIGNAL, N_OBS, SERROR
;
; CATEGORY:
;        fits I/O for Healpix
;
; CALLING SEQUENCE:
;      write_fits_cut4, Filename, Pixel, Signal, [N_obs, Serror, HDR=, XHDR=, 
;                       COORDSYS=, NESTED=, RING=, ORDERING=, NSIDE =,
;                       UNITS=, EXTENSION=, POLAR=,  HELP]]
;
; INPUTS:
;       Filename : STRING scalar,      
;                  output filename
;
;       Pixel    : INT, LONG or LONG64 vector, 
;                  index of Healpix pixel
;
;       Signal   : FLOAT or DOUBLE vector (same length as Pixel)
;                  signal in each pixel (eg, average of many measurements)
;
; OPTIONAL INPUTS:
;       N_obs    : INT, LONG or LONG64 vector (same length as Pixel)
;                  number of measurements in each pixel
;         if absent, N_OBS will be set uniformly to 1 in the output file
;         if set to a scalar constant, N_OBS will be set uniformly to this constant
;      
;       Serror  : FLOAT or DOUBLE vector (same length as Pixel)
;                 statistical error on Signal in each Pixel
;                 = rms of individual measurement / sqrt(N_obs)
;         if absent, SERROR will be set uniformly to 0.0 in the output file
;         if set to a scalar constant, SERROR will be set uniformly to this constant
;      
; KEYWORD PARAMETERS:
;       HDR=: primary header
;
;       XHDR=: STRING vector 
;                FITS like header to be inserted in the FITS file
;                and containing extra information on e.g, data origin, data
;                processing ...
;
;       COORDSYS=: STRING scalar, either 'C', 'G' or 'E' 
;                  coordinate system
;       NESTED=: if set, the Healpix NESTED scheme is used
;       RING=: if set, the Healpix RING scheme is used
;       ORDERING=: STRING scalar, either 'RING' or 'NESTED'
;            the information on the ORDERING scheme has to be given to the routine, either using
;            the keywords NESTED, RING or ORDERING, or filling in the Header
;            with the relevant keywords
;       
;       NSIDE=: Healpix resolution parameter
;            the information on the resolution NSIDE has to be given to the routine, either using
;            the keyword NSIDE, or filling in the Header with the relevant keywords
;       
;       UNITS=: STRING scalar, units of both Signal and Serror
;       
;       EXTENSION=: (0 based) extension number in which to write data
;           default = 0
;       
;       POLARISATION=: specify that file will contain the I, Q and U polarisation
;          Stokes parameter in extensions 0, 1 and 2 respectively
;       
;       HELP=: if set, an extensive help (this IDL header) is printed
;
; SIDE EFFECTS:
;        write a FITS file
;
; RESTRICTIONS:
;
; PROCEDURE:
;        calls write_fits_sb
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;        Sept 2000, EH, Caltech
;        Aug  2002 : uses Xhdr keyword to match other routines
;        Nov  2006: set default for Extension
;        Mar  2008: can deal with large nside, and double precision signal
;        Aug  2016: added HELP keyword
;
;-


routine = 'WRITE_FITS_CUT4'
if (keyword_set(help)) then begin
    doc_library,routine
    return
endif

syntax1 = 'SYNTAX: '+routine+', filename, pixel, signal, [n_obs, serror,  '
syntax2 = '                   HDR=, XHDR=, COORDSYS=, NESTED=, RING=, ORDERING=, NSIDE =, '
syntax3 = '                   UNITS=, EXTENSION=, POLARISATION=, HELP= ]'
if n_params() lt 3 then begin
    print,syntax1,syntax2,syntax3,format='(a)'
    if n_params() ne 0 then print,' ********* file not written ! ***********'
    return
endif

; check variable type
if undefined(info_header) then info_header = [' ']
if (size(/tname,filename) ne 'STRING' or $
    size(/tname,pixel) eq 'STRING' or $
    size(/tname,n_obs) eq 'STRING' or $
    size(/tname,signal) eq 'STRING' or $
    size(/tname,serror) eq 'STRING' or $
    size(/tname,info_header) ne 'STRING') then begin
    print,' wrong argument ordering:'
    print,syntax1,syntax2,syntax3,format='(a)'
    print,' ********* file not written ! ***********'
    return
endif

sigtype = size(/type,signal)
if (sigtype ge 6 && sigtype le 9) then begin
    print,' signal can not be complex or a string or a structure'
    print,syntax1,syntax2,syntax3,format='(a)'
    print,' ********* file not written ! ***********'
    return
endif

; check that nside is present
c_nside = 0
if defined(info_header) then nside_fits = sxpar(info_header,'NSIDE',count=c_nside)
if (undefined(nside_usr) and c_nside eq 0) then begin
    message,'NSIDE information is required (either as a keyword or in the FITS header)'
endif
nside = defined(nside_usr) ? nside_usr : nside_fits

; default value for N_obs and Serror
npix = n_elements(pixel)
nn_obs = n_elements(n_obs)
nserror = n_elements(serror)
if nn_obs  eq 0 then n_obs = replicate(1L,npix)
if nn_obs  eq 1 then n_obs = replicate(long(n_obs[0]),npix)
if nserror eq 0 then serror = replicate(0.,npix)
if nserror eq 1 then serror = replicate(float(serror[0]),npix)


; check consistency
nn_obs  = n_elements(n_obs)
nserror = n_elements(serror)
nsignal = n_elements(signal)
if (nsignal ne npix or nn_obs ne npix or nserror ne npix) then begin
    print,'inconsistent argument length in '+routine
endif

local_header = info_header
; insert units
if defined(units) then begin
    add_units_fits, local_header, units='      ', col=1, err=err
    add_units_fits, local_header, units=units[0], col=2, err=err
    add_units_fits, local_header, units='      ', col=3, err=err
    add_units_fits, local_header, units=units[0], col=4, err=err
    if (err ne 0) then message,'Error while writing header'
endif

add_ordering_fits, local_header, nested=nested, ring=ring, ordering=ordering,error=error
if (error ne 0) then message,'Error while writing header'

; insert polarisation related information if more than one extension
extname = ['TEMPERATURE', 'Q_POLARISATION', 'U_POLARISATION']
if undefined(extension_id) then extension_id = 0
if ((extension_id ge 0 and extension_id lt 3) or keyword_set(polar_in)) then begin
    sxaddpar,local_header, 'POLAR','T'
    sxaddpar,local_header, 'EXTNAME',extname[extension_id],after='TFIELDS'
    if (extension_id gt 0) then begin
        ; check consistency of NSIDE and ORDERING among extensions
        hdr1 = headfits(filename, ext=1, /silent)
        nside1 = long(sxpar(hdr1,'NSIDE'))
        if (nside1 ne nside) then begin
            print, 'NSIDE = ',nside1, nside
            message,'expects same NSIDE for all extensions of FITS file: '+strtrim(filename)
        endif
        ord1 = strupcase(strtrim(sxpar(hdr1,'ORDERING'),2))
        ord  = strupcase(strtrim(sxpar(local_header,'ORDERING'),2))
        if (ord1 ne ord) then begin
            print, 'ORDERING = ',ord1,ord, form="(3a10)"
            message,'expects same ORDERING for all extensions of FITS file: '+strtrim(filename)
        endif
    endif
endif

; create structures
prim_st = defined(prim_header) ? {HDR:prim_header} : 0
; exten_st = create_struct('HDR',local_header,'PIXEL',round(pixel),'SIGNAL',float(signal),'N_OBS',round(n_obs),'SERROR',float(serror))
exten_st = create_struct('HDR',local_header,'PIXEL',round(pixel),'SIGNAL',signal*1.0,'N_OBS',round(n_obs),'SERROR',serror*1.0)


; write file
write_fits_sb, filename, prim_st, exten_st,  $
  Coordsys=coordsys, Nside = nside_usr, /partial, extension=extension_id 


return
end


