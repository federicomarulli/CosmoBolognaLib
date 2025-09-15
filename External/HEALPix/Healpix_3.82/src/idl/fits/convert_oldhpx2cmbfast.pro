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
; -*- idl -*-
pro convert_oldhpx2cmbfast, file_in, file_out, no_renorm=no_renorm
;+
; NAME:
;   convert_oldhpx2cmbfast
;
; PURPOSE:
;   change normalization of polarization spectra in FITS file
;   from Healpix 1.1 convention to Healpix 1.2 (which is the same as cmbfast).
;   A keyword POLNORM = CMBFAST is added to the header to keep track of which
;   files have been processed.
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;   convert_oldhpx2cmbfast, file_in [, file_out, no_renorm=no_renorm]
; 
; INPUTS:
;   file_in : FITS file to be converted
;
; OPTIONAL INPUTS:
;   file_out : FITS file where new power spectra will be written.
;     if absent, the new file in written in file_in
;      
; KEYWORD PARAMETERS:
;   no_renorm : if set, the renormalization is not done.
;      but the keyword POLNORM = CMBFAST is added to the FITS header
;   (useful if the FITS file is already in CMBFAST format)
;
; OUTPUTS:
;   none
;
; OPTIONAL OUTPUTS:
;   none
;
; COMMON BLOCKS:
;   none
;
; SIDE EFFECTS:
;   a new file is written, or an existing file is overwritten
;
; RESTRICTIONS:
;   none
;
; PROCEDURE:
;   read, (renorm), add keyword, write
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   v1.0, EH, Feb 5 2003
;
;-

routine = 'convert_oldhpx2cmbfast'
syntax = routine+', file_in [, file_out, no_renorm=no_renorm]'

if (n_params() eq 0) then begin
    print,syntax
    return
endif

if (n_params() gt 2) then begin
    print,syntax
    message,'wrong number of arguments'
endif

if (n_params() eq 1) then file_out = file_in


fits2cl, clin, file_in, hdr=hdr, xhdr=xhdr, /silent
polnorm = strtrim(strupcase(sxpar(xhdr,'POLNORM',count=count)),2)
if (count eq 0) then polnorm = ''

clout = clin
addkey = 0
actual = 0
if (polnorm ne 'CMBFAST') then begin
     if keyword_set(no_renorm) then begin
         print,'No renormalization done for '+file_in
         addkey = 1
     endif else begin
        ncols = n_elements(clin[0,*])
        case (ncols) of
             1 : begin
                 print,'No polarization in file '+file_in
                 print,'no renormalization necessary'
             end
             4 : begin
                 print,'Renormalizing power spectra in '+file_in
                 print,' T -> T, E -> 2*E, B -> 2*B, TE -> sqrt(2)*TE '
                 addkey = 1
                 ; 0 = T : left as is
                 clout[*,1] = clout[*,1] * 2.  ; E or G
                 clout[*,2] = clout[*,2] * 2.  ; B or C
                 clout[*,3] = clout[*,3] * sqrt(2.) ; TE or CG
                 actual = 1
             end
             else : begin
                 print,'Number of columns : ',ncols,' in '+file_in
                 print,'I do not know how to renormalize spectra'
                 print,'File left unchanged'
                 return
             end
        endcase
    endelse
endif else begin ; already renormalized
    print,'File '+file_in+' is tagged as already renormalized.'
    print,'Spectra left unchanged'
endelse


if (actual eq 1) then begin
    sxaddpar, xhdr, 'HISTORY', 'Renormalization of pol. spectra done by '+routine
    sxaddpar, xhdr, 'HISTORY', 'on '+today_fits()
endif

print,'Writting output in '+file_out
cl2fits, clout, file_out, hdr=hdr, xhdr=xhdr, cmbfast=addkey

return
end
