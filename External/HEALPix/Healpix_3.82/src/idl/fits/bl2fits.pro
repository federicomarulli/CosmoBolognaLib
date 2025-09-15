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
PRO BL2FITS, bl_array, fitsfile, HDR = hdr, HELP = help, XHDR = xhdr

;+
; NAME:
;       BL2FITS
;
; PURPOSE:
;       write into a FITS file as an ascii table extension the beam window function 
;       passed to the routine . Adds additional headers if required.
;
; CALLING SEQUENCE:
;       BL2FITS, bl_array, fitsfile, [HDR = , HELP =, XHDR =]
; 
; INPUTS:
;       bl_array = real array of B(l) coefficients to be written to file.
;                  This has dimension:
;                    (lmax+1,1) for Temperature only   
;                      or
;                    (lmax+1,2) given in the sequence T E
;                    (lmax+1,3) given in the sequence T E B
;       Note: B(l) is a window function (most likely of an azymuthally
;       symmetrical beam), NOT a power spectrum, and its sign is arbitrary.
;
;       fitsfile = String containing the name of the file to be written      
;
; OPTIONAL INPUTS:
;
;
; OPTIONAL INPUT KEYWORDS:
;       HELP = if set, this header is printed out, no file is written
;
;       HDR      = String array containing the (non-trivial) header for
;                  the FITS file.
;       XHDR     = String array containing the (non-trivial) extension header 
;                  for the extension.
;                  NOTE: trivial keywords (such as 'BITPIX', 'NAXIS',
;                  'TFORM*', ...) 
;                  present in user-provided headers will be IGNORED
;                  and replaced by routine-generated values.
; PROCEDURES CALLED:
;      wrapper to CL2FITS, /BEAMWINDOW
;
; MODIFICATION HISTORY:
;       May 2007. EH: created
;-

code = 'BL2FITS'
syntax = 'Syntax : '+code+', bl_array, fitsfile, [HDR = , HELP = , XHDR = ] '

if keyword_set(help) then begin
      doc_library,code
      return
endif

if N_params() NE 2 then begin
      print,syntax
      print,'   file NOT written '
      goto, Exit
endif

if datatype(bl_array) eq 'STR' or datatype(fitsfile) ne 'STR' then begin
      print,syntax
      print,'   the array comes first, then the file name '
      print,'   file NOT written '
      goto, Exit
endif

cl2fits, bl_array, fitsfile, HDR = hdr, HELP = help, XHDR = xhdr, /BEAMWINDOW

Exit:
return
end



