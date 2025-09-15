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
pro write_fits_partial, filename, pixel, iqu,  $
                     hdr=prim_header, xhdr = info_header, $
                     Coordsys=coordsys, $
                     Nested=nested, Ring=ring, Ordering=ordering, $
                     Nside = nside_usr, $
                     Units = units, Help=help, colnames=colnames_usr, verbose=verbose, $
                     Extension = extension_id
;+
; NAME:
;         write_fits_partial
;
; PURPOSE:
;        writes a FITS file for incomplete (polarized) sky coverage with a 4 column format
;        PIXEL, TEMPERATURE, Q_POLARISATION, U_POLARISATION
;
; CATEGORY:
;        fits I/O for Healpix
;
; CALLING SEQUENCE:
;      write_fits_partial, Filename, Pixel, IQU, [COLNAMES=, COORDSYS=, 
;                       EXTENSION=,  HDR=, HELP, 
;                       NESTED=, NSIDE =, ORDERING=, RING=, 
;                       UNITS=, VERBOSE=,  XHDR=]
;
; INPUTS:
;       Filename : STRING scalar,      
;                  output filename
;
;       Pixel    : INT, LONG or LONG64 vector of length np
;                  index of Healpix pixel
;
;       IQU   : FLOAT or DOUBLE array of size (np, nc)
;                  I,Q and U Stokes parameters for each pixel (if nc=3)
;                 or I of each pixel (if nc=1)
;
;      
; KEYWORD PARAMETERS:
;       COLNAMES=: STRING vector with column names (beside PIXEL) (not case sensitive: [A-Z,0-9,_])
;         Default:  TEMPERATURE (for 1 column), 
;                   TEMPERATURE, Q_POLARISATION, U_POLARISATION  (for 3 columns)
;                   C01, C02, C03, C04, ... (otherwise)
;         If provided the number of COLNAMES must be >= the number of columns
;
;       COORDSYS=: STRING scalar, either 'C', 'G' or 'E' 
;                  coordinate system
;
;       EXTENSION=: (0 based) extension number in which to write data
;           default = 0
;
;       HDR=: primary header
;
;       HELP=: if set, an extensive help (this IDL header) is printed
;
;       NESTED=: if set, the Healpix NESTED scheme is used
;
;       ORDERING=: STRING scalar, either 'RING' or 'NESTED'
;            the information on the ORDERING scheme has to be given to the routine, either using
;            the keywords NESTED, RING or ORDERING, or filling in the Header
;            with the relevant keywords
;
;       RING=: if set, the Healpix RING scheme is used
;       
;       NSIDE=: Healpix resolution parameter
;            the information on the resolution NSIDE has to be given to the routine, either using
;            the keyword NSIDE, or filling in the Header with the relevant keywords
;       
;       UNITS=: STRING scalar or vector defining physical units of columns (except PIXEL)
;          if scalar, same for all
;          if vector, each column can have its own units
;          if needed, the last UNITS provided will be replicated for the remaining columns
;
;       VERBOSE=: if set, be verbose while the file is written
;      
;       XHDR=: STRING vector 
;                FITS like header to be inserted in the FITS file
;                and containing extra information on e.g, data origin, data
;                processing ...
;       
;       VERBOSE=: if set, be verbose
;
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
;       nside = 512
;       pixel = lindgen(nside2npix(nside)/10)
;       write_fits_partial,'part_T.fits',  pixel, pixel*100.,            nside=nside, units='K'
;         ; creates a FITS binary table with the columns PIXEL and TEMPERATURE

;       write_fits_partial,'part_TQU.fits',pixel, pixel#[100.,1.,1.],    nside=nside, units='K'
;         ; creates a FITS binary table with the columns PIXEL, TEMPERATURE, Q_POLARISATION and U_POLARISATION

;       write_fits_partial,'part_xxx.fits',pixel, pixel#[1.,-2.,3.,-4.], nside=nside, colnames=['c1',B2','xx','POWER'],units=['K','m','s','W']
;         ; creates a FITS binary table with the columns PIXEL, C1, B2, XX and POWER
;
; MODIFICATION HISTORY:
;     Jan 2020: EH, IAP, adapted from write_fits_cut4
;
;-


routine = 'WRITE_FITS_PARTIAL'
if (keyword_set(help)) then begin
    doc_library,routine
    return
endif

syntax1 = 'SYNTAX: '+routine+', filename, pixel, IQU, '
syntax2 = '                   [HDR=, XHDR=, COLNAMES=, COORDSYS=, EXTENSION=, HELP=, '
syntax3 = '                   NESTED=, NSIDE =, ORDERING=, RING=, UNITS=, VERBOSE= ]'
if n_params() lt 3 then begin
    print,syntax1,syntax2,syntax3,format='(a)'
    if n_params() ne 0 then print,' ********* file not written ! ***********'
    return
endif

; check variable type
if undefined(info_header) then info_header = [' ']
if (size(/tname,filename) ne 'STRING' or $
    size(/tname,pixel) eq 'STRING' or $
    size(/tname,iqu) eq 'STRING' or $
    size(/tname,info_header) ne 'STRING') then begin
    print,' wrong argument ordering:'
    print,syntax1,syntax2,syntax3,format='(a)'
    print,' ********* file not written ! ***********'
    return
endif

sigtype = size(/type,iqu)
if (sigtype ge 6 && sigtype le 9) then begin
    print,' iqu can not be complex or a string or a structure'
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


; check pixel
pix_type = size(/tname, pixel)
pix_min  = min(pixel, max=pix_max)
if ( (pix_type ne 'INT' && pix_type ne 'LONG' && pix_type ne 'LONG64') || (pix_min lt 0)) then begin
    message,/info,'PIXEL array must be a 2, 4 or 8-byte signed integer >=0'
    message,/info,'Currently: type = '+pix_type
    message,/info,'           min  = '+strtrim(pix_min)
    message,/info,'           max  = '+strtrim(pix_max)
    message,'Abort'
endif

; check consistency
npix = n_elements(pixel)
siqu = size(iqu,/dim)
nc = n_elements(iqu)/siqu[0]
if (npix ne siqu[0]) then begin
    print,'size of IQU does not match number of pixels in '+routine
    print,syntax1,syntax2,syntax3,format='(a)'
    print,' ********* file not written ! ***********'
    return
endif


local_header = info_header

; remove TFORM*, TTYPE*, TUNIT*
sxdelpar, local_header, ['TFORM1','TTYPE1','TUNIT1',$
                         'TFORM2','TTYPE2','TUNIT2',$
                         'TFORM3','TTYPE3','TUNIT3',$
                         'TFORM4','TTYPE4','TUNIT4']

; insert units
add_units_fits, local_header, units='      ', col=1, err=err ; no units for PIXEL
if defined(units) then begin
    final_units = units[ lindgen(nc) < (n_elements(units)-1) ]
    for icol=0,nc-1 do begin
        add_units_fits, local_header, units=final_units[icol], col=2+icol, err=err
    endfor
    if (err ne 0) then message,'Error while writing header'
endif

add_ordering_fits, local_header, nested=nested, ring=ring, ordering=ordering,error=error
if (error ne 0) then message,'Error while writing header'

; determine column names
if defined(colnames_usr) then begin
    colnames = colnames_usr
    if (n_elements(colnames) lt nc) then begin
        message,'Number of colnames must be as large as number of columns to write'
        print,syntax1,syntax2,syntax3,format='(a)'
        print,' ********* file not written ! ***********'
        return
    endif
endif else begin
    colnames = ['TEMPERATURE', 'Q_POLARISATION', 'U_POLARISATION']
    if (nc eq 3) then sxaddpar, local_header, 'POLAR', 'T'
    if (nc eq 1) then sxaddpar, local_header, 'POLAR', 'F'
    if (nc ne 1 && nc ne 3) then begin
        colnames = 'C'+string(lindgen(nc)+1,form='(i2.2)') ; [C01, C02, C03, ...]
    endif
endelse

if keyword_set(verbose) then begin
    print,'Writing '+filename+' of '+strtrim(siqu[0],2)+ ' rows with columns (and units)'
    print,['PIXEL',strupcase(colnames)],form='(30(a15,:))'
    if defined(final_units) then print,[' ',final_units],form='(30(a15,:))'
endif

; create structures
prim_st  = defined(prim_header) ? {HDR:prim_header} : 0
exten_st = {HDR: local_header, PIXEL: pixel}
for ic=0,nc-1 do begin
    exten_st = create_struct(exten_st, colnames[ic], iqu[*,ic])
endfor

; write file
write_fits_sb, filename, prim_st, exten_st,  $
  Coordsys=coordsys, Nside = nside_usr, /partial, extension=extension_id 


return
end


