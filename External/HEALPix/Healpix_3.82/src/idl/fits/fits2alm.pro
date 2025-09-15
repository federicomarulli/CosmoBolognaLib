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

;+
; NAME:
;       FITS2ALM
;
; PURPOSE:
;       read from a FITS file a binary table extension containing spherical
;       harmonic (scalar or tensor) coefficients together with their index 
;       and return them. Reads header information if required.
;
; CALLING SEQUENCE:
;       FITS2ALM, index, alm_array, fitsfile, 
;          [signal, HDR=, LMAX=, LMIN=, XHDR=, HELP= ]
; 
; INPUTS:
;       fitsfile = String containing the name of the file to be written      
;
; OPTIONAL INPUT:
;       signal   = String defining the signal coefficients to read
;                  Valid options: 'T' 'E' 'B' 'ALL' 
;                  Default: 'T'
;
; OUTPUTS:
;       index    = Integer array containing the index for the corresponding 
;                  array of alm coefficients (and erralm if required). The
;                  index {i} is related to {l,m} by the relation
;                  i = l^2 + l + m + 1
;      alm_array = real/double array of alm coefficients read from the file.
;                  This has dimension (nl,nalm,nsig) -- corresponding to
;                  nl   = number of {l,m} indices
;                  nalm = 2 for real and imaginary parts of alm coefficients;
;                         4 for above plus corresponding error values
;                  nsig = number of signals extracted (1 for any of T E B
;                         or 3 if ALL extracted). Each signal is stored
;                         in a separate extension.
;
; OPTIONAL INPUT KEYWORDS:
;       HELP=      If set, shows this documentation header
;       LMAX=      Largest multipole  to be output, among those available
;       LMIN=      Smallest multipole to be output, among those available
;              If LMIN (LMAX) is below (above) the range of l's present in the file,
;              it will be silently ignored
;
; OPTIONAL OUTPUT KEYWORDS:
;       HDR=      String array containing the header for the FITS file.
;       XHDR=     String array containing the extension header(s). If 
;                  ALL signals are required, then the three extension 
;                  headers are returned appended into one string array.
;
; EXAMPLE:
;       Read the B tensor harmonic coefficients into the real array, alm,
;       from the FITS file, coeffs.fits, with extension header information
;       contained in the string array, ext_text
;
;       IDL> fits2alm,alm,'coeffs.fits','B',XHDR=ext_txt
;
;
; PROCEDURES CALLED:
;       HEADFITS,  FITS_INFO, FITS_READ, TBINFO, TBGET
;
; MODIFICATION HISTORY:
;       May 1999: written by A.J. Banday (MPA)     
;       Dec 2004, EH: edited to avoid faulty /use_colnum keyword in MRDFITS  
;       Feb 2005, EH: replaced MRDFITS by faster FITS_READ+TBINFO+TBGET
;       May 2005, EH, replaces FINDFILE by FILE_TEST
;       Aug 2005, EH: make output alm_array of same precision as FITS file data
;  Jan 2008, EH: calls tbfree to remove heap pointer created by TBINFO
;       Jan 2009: calls init_astrolib
;       Nov 2009: added LMIN, LMAX and HELP
;       Feb 2012: issues warning if an extension is not found, 
;         and crashes if none of the requested data is available.
;
; requires the THE IDL ASTRONOMY USER'S LIBRARY 
; that can be found at http://idlastro.gsfc.nasa.gov/homepage.html
;
;-
pro select_good_index_alm, index, good, ngood=ngood, lmax=lmax, lmin=lmin, do_select=do_select

; useful range of index
llmin = keyword_set(lmin) ? long64(lmin[0]) : 0
llmax = keyword_set(lmax) ? long64(lmax[0]) : 10LL^7
my_range = [ llmin*llmin + 1, (llmax+1)^2]

; find range intersection
file_range = hpx_minmax(index)
final_range = [my_range[0]>file_range[0], my_range[1]<file_range[1]]

; select alm in specified range, only if necessary
do_select = ~array_equal(file_range, final_range)
if (do_select) then begin
    good = where( index ge final_range[0] and index le final_range[1], ngood)
endif else begin
    ngood = n_elements(index)
endelse

if (ngood eq 0) then begin
    index2lm, final_range, ll, mm
    message_patch,'No alm''s found in specified ell range '+string(ll[0])+string(ll[1]),level=-1
endif

return
end

;--------------------------

PRO FITS2ALM, index, alm_array, fitsfile, signal, $
              HDR = hdr, XHDR = xhdr, LMIN = lmin, LMAX = lmax, HELP=help

if (keyword_set(help)) then begin
    doc_library,'fits2alm'
    return
endif

if N_params() LT 3 or N_params() gt 4 then begin
      print,'Syntax : FITS2ALM, index, alm_array, fitsfile, [signal, HDR=, HELP=, LMAX=, LMIN=, XHDR=] '
      goto, Exit
endif

if (~file_test(fitsfile)) then message,'file '+fitsfile+' not found'

; run astrolib routine to set up non-standard system variables
init_astrolib

hdr  = HEADFITS(fitsfile)

; -------- extension -----------

usignal = defined(signal) ? strupcase(signal) : 'T'

CASE usignal OF
    'T'  : BEGIN
        extension = 1 
        nsig = 1
    END
    'E'  : BEGIN
        extension = 2
        nsig = 1
    END
    'B'  : BEGIN
        extension = 3 
        nsig = 1
    END
    'ALL': BEGIN
        extension = 1 
        nsig = 3                ; extension value initialises read-loop
    END
    
    else: BEGIN
        message,/info,' Incorrect signal selected: '+signal
        message,' File '+fitsfile+' not read.'
        ;goto, Exit
    END
ENDCASE


; count extensions
fits_info,fitsfile, /silent, n_ext=n_ext

; simply read the extensions from the FITS file
savehdr = ''
nrows_old = -1
nsig_eff = 0
for i = 0,nsig-1 do begin
    exten = extension+i
    if (exten gt n_ext) then begin
        message,/info,' WARNING: Required extension does not exist in file: '+fitsfile
        message,/info,'          Output data will be shorter than expected.'
        goto, skip
    endif else begin
        nsig_eff +=1
    endelse
    

    ; read data
    fits_read, fitsfile, tmpout, xhdr, /no_pdu, exten_no = exten
    nrows = sxpar(xhdr,'NAXIS2')
    savehdr = [savehdr,xhdr]
  
    if (i eq 0) then begin
        ; first extension, create arrays
        tbinfo, xhdr, tab_xhdr
        ncols = n_elements(tab_xhdr.(0))

        index = tbget(tab_xhdr, tmpout, 1)
        select_good_index_alm, index, good, ngood=ngood, lmax=lmax, lmin=lmin, do_select=do_select
        if (do_select) then index = index[good]
        alm_array = make_array(ngood, ncols-1, nsig, type=(tab_xhdr.idltype)[1])
    endif else begin
        ; other extensions, check consistency
;         if (nrows ne nrows_old) then begin
;             print,exten-1,nrows_old,exten,nrows,   $
;               form='("#",i1,":", i8,",  #",i1,":",i8)'
;             message,' ERROR: Extensions have different sizes'
;         endif
        indtmp = tbget(tab_xhdr, tmpout, 1)
        select_good_index_alm, indtmp, good, ngood=ngood, lmax=lmax, lmin=lmin, do_select=do_select
        if (do_select) then indtmp = indtmp[good]
        if (~array_equal(index, indtmp)) then begin
            message,' ERROR: Alm''s defined for different (l,m)'
        endif
    endelse
;     nrows_old = nrows

    if (do_select) then begin
        for jc=2,ncols do alm_array[*,jc-2,i] = (tbget(tab_xhdr, tmpout, jc))[good]
    endif else begin
        for jc=2,ncols do alm_array[*,jc-2,i] =  tbget(tab_xhdr, tmpout, jc)
    endelse

    skip:

endfor

if (nsig_eff eq 0) then begin
    message,'ERROR: None of the requested data found. Aborting'
endif

if (nsig_eff lt nsig) then alm_array = alm_array[*,*,0:nsig_eff-1]
alm_array = reform(alm_array, /Overwrite)
xhdr = savehdr
tbfree, tab_xhdr

; Exit routine
Exit:
return
end



