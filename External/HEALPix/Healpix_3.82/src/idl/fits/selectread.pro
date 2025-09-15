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
pro selectread, file, array, polvec, header=exthdr, columns=columns, extension=extension_id, poltype=poltype, tonan=tonan, offset=offset, factor=factor, flip = flip, no_pdu= no_pdu, help=help
;+
; NAME:
;    selectread
;
; PURPOSE:
;    routine to read a FITS file containing temperature or temperature+
;    polarisation data. It has been optimized for speed and low memory usage
;
; CATEGORY:
;    healpix I/O
;
; CALLING SEQUENCE:
;    selectread, File, Array, [Polvec, HEADER=, COLUMNS=, EXTENSION=, POLTYPE=,
;    TONAN=, OFFSET=, FACTOR=, FLIP=, NO_PDU=, HELP=]
;
; INPUTS:
;     File : string containing FITS file name
;
; OPTIONAL INPUTS:
;
;      
; KEYWORD PARAMETERS:
;     COLUMNS= list of numbers specifying columns to be read (1 based).
;       By default all columns are read, unless POLTYPE is set, in which case
;       only the relevant columns are read and processed
;
;     EXTENSION : extension unit to be read from FITS file: 
;       either its 0-based ID number (ie, 0 for first extension after primary array)
;       or the case-insensitive value of its EXTNAME keyword.
;	If absent, first extension (=0) will be read
;
;     POLTYPE= processing to apply to polarised data before output
;        0: do nothing: output selected column(s) as read (after application of
;               OFFSET and FACTOR to each columns)
;        1: read the Q and U columns of the FITS file and outputs in Array the
;               norm = sqrt(U^2+Q^2) * factor + offset * factor
;        2: read the Q and U columns of the FITS file and output in Array the
;              angle psi = 1/2 atan(U/Q)
;        3: read the T,Q and U columns and outputs 
;              in Array (T + offset) * factor
;              polvec will contain 
;                 norm=sqrt(U^2+Q^2) * factor    and 
;                 psi=1/2 atan(U/Q)
;
;     TONAN= set flagged pixel to NaN. Is set when FACTOR and/or OFFSET are set
;
;     FACTOR= multiplicative factor applied to data, default = 1.
;
;     OFFSET= additive factor applied to data 
;      (ie, data_out = (data_fits + offset) * factor
;
;     FLIP= flip orientation of polarisation vector 
;      (ie, change the sign of psi, only applies when POLTYPE=2 or 3)
;
;     NO_PDU= if set, the primary data unit header is not included in the
;     returned header
;
;     HELP= if set, print out this help header and exit
;
; OUTPUTS:
;     Array: array of size (np, ncols) 
;        where np is the number of data points, ncols is the number of columns
;        read, and the array type is determined by the data read.
;
; OPTIONAL OUTPUTS:
;     Polvec: if POLTYPE=3, will contain the polarisation norm and direction in
;       an array of size (np, 2)
;
;     HEADER: contains the FITS header (combination of the primary unit header
;     and the one of the extension being read, unless NO_PDU is set)
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;         v1.0, EH, 2005-02-07
;                   2005-08-29: added no pdu
;  Jan 2008, EH: calls tbfree to remove heap pointer created by TBINFO
;  June 2008: EH, modified to deal with files with very large TFORM
;                 corrected documentation header wrt offset/factor
;  Oct 2008, EH: allows offsetting of polarization norm when poltype=1
;  Dec 2008, EH: use OR instead of '||' for non-scalar tests
;  Nov 2009, EH: accepts non-scalar Offset and Factor
;  Jan 2013, EH: accepts string type EXTENSION keyword
;  Mar 2013, EH: added HELP keyword
;  Jul 2014, EH: FACTOR and OFFSET no more required (and take default value) when TONAN is set
;-

code = 'selectread'
syntax = [code+', File, Array, [Polvec, HEADER=, COLUMNS=, EXTENSION=, EXTNAME= POLTYPE=,','                                  TONAN=, OFFSET=, FACTOR=, FLIP=, NO_PDU=, HELP=]']

if keyword_set(help) then begin
    doc_library,code
    return
endif

if (n_params() lt 2) then begin
    print,syntax,form='(a)'
    return
endif

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

do_rescale = (keyword_set(tonan) || keyword_set(offset) || keyword_set(factor))
polar = (keyword_set(poltype))      ? poltype          : 0
flipconv = (keyword_set(flip))      ? 1                : -1

my_factor = defined(factor) ? factor : [1]
my_offset = defined(offset) ? offset : [0]

; open file
fits_open, file, fcb

; identify right extension
if size(extension_id,/tname) eq 'STRING' then begin ; by name, if provided
    sxn = strupcase(extension_id)
    xtnum = (where(strupcase(fcb.extname) eq sxn, count))[0] ; 1 based extension number
    if count gt 1 then begin
        print,'WARNING: '+strtrim(count,2)+' extensions matching '+sxn+' found in '+file+'. The first one will be read.'
    endif
    if count eq 0 then begin
        message,/info,'Valid (case un-sensitive) EXTNAME choices in '+file+' are'
        print,' 0-based #          EXTNAME                    (XTENSION):'
        for i=1,fcb.nextend do print,string(i-1,'(i6)')+' '+string(fcb.extname[i],'(a30)')$
          +'          ('+fcb.xtension[i]+')'
        message,/info,'Extension matching '+sxn+' not found in '+file
        message,'Aborting.'
    endif
endif else begin ; otherwise by number, taking the first one by default
    xtnum   = (keyword_set(extension_id)) ? (extension_id+1) : 1 ; 1 based extension number
    if (xtnum gt fcb.nextend) then message,'Only '+strtrim(fcb.nextend,2)+' extensions in '+file
endelse

; find number of words in extension
n_wpr  = (fcb.axis)[0,xtnum] ; words per row in extension
n_rows = long64((fcb.axis)[1,xtnum]) ; number of rows
n_words = n_wpr * n_rows

; read header for extension + optionally PDU
fits_read, fcb, void, exthdr, exten_no = xtnum, /header_only, no_pdu = no_pdu
tfields  = round(float(sxpar(exthdr,'TFIELDS')))
bad_data =       float(sxpar(exthdr,'BAD_DATA', count=nbd))
if (nbd eq 0) then bad_data = !healpix.bad_value
polcconv = strtrim(strupcase(sxpar(exthdr,'POLCCONV',count=n_pccv)),2)
if (n_pccv eq 0) then polcconv=''

; read image if file is in deprecated format
if (xtnum eq 0 && fcb.nextend eq 0) then begin
    if polar gt 0 then message,'no polarisation information found in '+file
    fits_read, fcb, array, exten_no = xtnum
    if (do_rescale) then begin
        bad_pixels = where(array le (bad_data*0.9) or finite(array,/nan), nbad)
        if (nbad gt 0)    then array[bad_pixels] = !values.f_nan
        if (my_factor[0] ne 1.) then array = temporary(array) * my_factor[0] 
        if (my_offset[0] ne 0.) then array = temporary(array) + (my_factor[0]*my_offset[0])
    endif
    goto, all_done
endif

; check out file for polarisation
if (tfields lt 3 && polar gt 0) then begin
    message,'no polarisation information found in '+file
endif

tbinfo, exthdr, tab_xhdr
nentry = max(tab_xhdr.numval)
npix = nentry * n_rows

; columns to be output in array
cols = indgen(tfields) + 1
if defined(columns) then cols = columns
if (max(cols) gt tfields || min(cols) le 0) then begin
    print,'columns = ',cols
    message,'invalid choice of columns for '+file
endif

nmaps = n_elements(cols)
if (polar eq 1 || polar eq 2) then begin
    cols = [2,3] ; read Q and U
    nmaps = 1
endif
if (polar eq 3) then begin
    cols = [1,2,3]
    nmaps = 1
endif
if (polar gt 3 || polar lt 0) then begin
    print,'Poltype = ',polar
    message,'invalid poltype for '+file
endif
types = (tab_xhdr.idltype)[cols-1]
type = max(types) ; find out type of output array to fit all data
junk = where(types eq 7, n_string)
if (n_string ge 1) then begin
    message,'Table in input file contains strings'
endif


; create array receiving final data
array = make_array(type=type, npix, nmaps, /nozero)
if (polar eq 3) then polvec = make_array(type=type, npix, 2, /nozero)

; read data piece by piece and process each piece individually
stride = 5.e6 > n_wpr ; 5 MB or 1 row per piece (June-2008)
;stride = 500.e6 > n_wpr ; 500 MB or 1 row per piece (Dec-2011)
stride = FLOOR(stride / n_wpr) * n_wpr
w_start = long64(0)
pstart = long64(0)
n_factors = n_elements(my_factor)
n_offsets = n_elements(my_offset)
while (w_start LE (n_words-1) ) do begin
    ; read one piece
    w_end = (w_start + stride - 1L) < (n_words-1)
    fits_read, fcb, data, exten_no = xtnum, first=w_start, last=w_end
    nr = (w_end - w_start + 1) / n_wpr ; number of rows read
    np = nr * nentry ; number of pixels read
    data = reform(data, n_wpr, nr, /overwrite) ; required by tbget
    ; select useful columns and/or process data read
    if (polar eq 0) then begin ; standard case
        for i=0,nmaps-1 do begin
            if (do_rescale) then begin
                x = (tbget(tab_xhdr, data, cols[i]))[*]
                bad_pixels = where(x le (bad_data*0.9) or finite(x,/nan), nbad)
                if (nbad gt 0)    then x[bad_pixels] = !values.f_nan
                id_fact = (n_factors eq nmaps) ? i : 0
                id_offs = (n_offsets eq nmaps) ? i : 0
                if (my_factor[id_fact] ne 1.) then x *=  my_factor[id_fact] 
                if (my_offset[id_offs] ne 0.) then x += (my_factor[id_fact]*my_offset[id_offs])
                array[pstart,i] = x
            endif else begin
                array[pstart,i] = (tbget(tab_xhdr, data, cols[i]))[*]
            endelse
        endfor
    endif else begin
        if (w_start eq 0) then begin
            tmp  = make_array(type=type,np,3,/nozero) 
            norm = make_array(type=type,np,  /nozero) & psi = norm
        endif
        for i=0,2 do begin
            tmp[0,i] = (tbget(tab_xhdr, data, i+1))[*]
            if (do_rescale) then begin
                bad_pixels = where(tmp[*,i] le (bad_data*0.9) or finite(tmp[*,i],/nan), nbad)
                if (nbad gt 0)    then tmp[bad_pixels,i] = !values.f_nan
                if (my_factor[0] ne 1.) then tmp[*,i] *= my_factor[0]  ; rescale all fields
                if (my_offset[0] ne 0. && i eq 0) then tmp[*,i] +=  (my_factor[0]*my_offset[0]) ; only offset temperature
            endif
        endfor
        if (polar eq 1 || polar eq 3) then begin
            norm[0:np-1] = sqrt(tmp[0:np-1,1]^2+tmp[0:np-1,2]^2)  ; (Q,U) --> P
            if (polar eq 1 && my_offset[0] ne 0.) then norm[0:np-1] += (my_factor[0]*my_offset[0]) ; offset P, only when plotting P alone
        endif
        if (polar eq 2 || polar eq 3) then begin
            psi[0:np-1] = 0.5 * atan(tmp[0:np-1,2]*flipconv, tmp[0:np-1,1]) ; (Q,U) --> psi
            primer = '(http://healpix.sourceforge.net/pdf/intro.pdf)'
            case polcconv of
                'COSMO':  ; COSMO -> carry on silently
                '': begin ; absent -> assume COSMO and make a comment
                    print,'The POLCCONV keyword was not found in '+file
                    print,'COSMO (HEALPix default) will be assumed, and map is unchanged.'
                    print,'See HEALPix primer '+primer+' for details.'
                end
                'IAU':begin ; IAU -> issue a warning
                    message,/info,'POLCCONV=IAU found in '+file,level=-2
                    message,/info,'The sign of U (and direction of polarisation) is changed',level=-2
                    message,/info,'See HEALPix primer '+primer+' for details.',level=-2
                    psi[0:np-1] *= -1
                end
                else: begin ; none of the above -> stop
                    message,/info,'POLCCONV='+polcconv+' found in '+file,level=-2
                    message,/info,'It is neither COSMO nor IAU. Aborting',level=-2
                    message,      'See HEALPix primer '+primer+' for details.',level=-2
                end
            endcase
        endif
        if (polar eq 1) then array[pstart] = norm[0:np-1]
        if (polar eq 2) then array[pstart] = psi[0:np-1]
        if (polar eq 3) then begin
            array[pstart,0] = tmp[0:np-1,0]
            polvec[pstart,0] = norm[0:np-1]
            polvec[pstart,1] = psi[0:np-1]
        endif
    endelse
    ; get ready for next piece
    w_start = w_end + 1
    pstart = pstart + np
endwhile
data = 0
tmp = 0 & norm = 0 & psi = 0
tbfree, tab_xhdr

; close file
all_done:
fits_close, fcb



return
end
