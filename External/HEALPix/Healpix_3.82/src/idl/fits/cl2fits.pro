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
PRO CL2FITS, cl_array, fitsfile, HDR = hdr, HELP = help, XHDR = xhdr, CMBFAST = cmbfast, UNITS = units, BEAMWINDOW = beamwindow

;+
; NAME:
;       CL2FITS
;
; PURPOSE:
;       write into a FITS file as an ascii table extension the power 
;       spectrum coefficients passed to the routine . Adds additional
;       headers if required.
;
; CALLING SEQUENCE:
;       CL2FITS, cl_array, fitsfile, [HDR = , HELP =, XHDR =, CMBFAST=, UNITS=, BEAMWINDOW=]
; 
; INPUTS:
;       cl_array = real array of Cl coefficients to be written to file.
;                  This has dimension:
;                    (lmax+1,1) for Temperature only   
;                      or
;                    (lmax+1,4) given in the sequence T E B TxE,
;                     for temperature and polarisation
;                    (lmax+1,6) given in the sequence T E B TxE TxB ExB
;                     for temperature and polarisation
;                  The convention for the power spectrum is 
;                  cl = Sum_m {a_lm^2}/(2l+1)
;                  ie. NOT normalised by the Harrison-Zeldovich spectrum.
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
;       CMBFAST  = if set, the keyword POLNORM = CMBFAST is added to the header
;
;       UNITS =    String scalar containing units of power spectrum (eg, uK^2,
;         Kelvin**2, ...), to be put in keyword 'TUNIT*' of extension Header. 
;         If provided, will override the values present in XHDR (if any).
;
;       BEAMWINDOW = if set, the data to be written is a beam window function
;         (which can be of arbitrary sign) instead of a power spectrum
;
; EXAMPLE
;       Write the power spectrum coefficients in the real array, pwr,
;       into the FITS file, spectrum.fits, with extension header information
;       contained in the string array, ext_text
;
;       IDL> cl2fits,pwr,'spectrum.fits',XHDR=ext_txt
;
;
; PROCEDURES CALLED:
;       WRITEFITS, HEADFITS, SXADDPAR, FTCREATE, FTPUT, TODAY_FITS(), SXDELPAR
;
; MODIFICATION HISTORY:
;       May 1999: written by A.J. Banday (MPA)        
;       Feb 2000: accept one column array for temperature only,
;                 give explicit name to each column
;                 put in today_fits()   EH (Caltech) 
;                 added EXTNAME and fit file size to actual number of fields
;       Sep 2001, EH: ignore trivial/ancillary keywords in user-provided
;                headers.
;       Jan 2003, EH: furiously faster by using vectorial call to FTPUT
;       Feb 2003, EH: added CMBFAST keyword
;       Aug 2004, EH: added HELP keyword
;       May 2005, EH: added keywords in FITS header: MAX-LPOL, POLAR, BCROSS, UNITS
;       May 2007. EH: added beamwindow
;       Oct 2008. EH: prints out 'BL2FITS' in /beamwindow mode
;       Jan 2009: calls init_astrolib
;       Aug 2015: can write files with 9 columns
;-

code = 'CL2FITS'
syntax = 'Syntax : '+code+', cl_array, fitsfile, [HDR = , HELP = , XHDR = , CMBFAST = , UNITS = ] '

if keyword_set(help) then begin
      doc_library,code
      return
endif

if N_params() NE 2 then begin
      print,syntax
      print,'   file NOT written '
      goto, Exit
endif

if datatype(cl_array) eq 'STR' or datatype(fitsfile) ne 'STR' then begin
      print,syntax
      print,'   the array comes first, then the file name '
      print,'   file NOT written '
      goto, Exit
endif

; run astrolib routine to set up non-standard system variables
init_astrolib

if n_elements(hdr)  eq 0 then hdr  = ' '
if n_elements(xhdr) eq 0 then xhdr = ' '

do_beam = keyword_set(beamwindow)

; ------- primary unit ----------
; opens the file, write a minimal header and close it
WRITEFITS,fitsfile,0

; update date to Y2K
h0 = HEADFITS(fitsfile)
fdate = today_fits()
SXADDPAR,h0,'DATE',fdate,' Creation date (CCYY-MM-DD) of FITS header'

; add header information given by user
if (STRMID(hdr(0),0,1) NE ' ') then begin
    ; remove "trivial" keywords from user defined header
    trivial = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND','DATE']
    hdr1 = hdr
    sxdelpar, hdr1, trivial
    iend = WHERE( STRUPCASE(STRMID(h0,0,8)) EQ 'END     ', nend)
    if (nend eq 1) then begin
        iend = iend(0)
        h1 = [h0(0:iend-1), hdr1] 
    endif 
endif else begin
    h1 = h0
endelse

; opens the file, write the header and the image if any, and close it
WRITEFITS, fitsfile, 0, h1

; -------- extension -----------

; get the dimensions of the input array
info  = size(cl_array)
ndims = info(0) ; # of array dimensions
nrows = info(1) ; # of entries for l-range; nrows = lmax+1
if ndims eq 1 then nsigs = 1 else $
  nsigs = info(2) ; # of signals; should be 1 corresponding to T
;                                        or 4 corresponding to T E B TxE

if (do_beam) then begin
    code='BL2FITS'
    if (nsigs gt 3) then begin
        print,code+': Input array must have dimensions (lmax+1,1), (lmax+1,2) or (lmax+1,3)'
        print,'  No file written'
        goto, Exit
    endif
endif else begin
    if (nsigs ne 9 and nsigs ne 6 and nsigs ne 4 and nsigs ne 1) then begin
        print,code+': Input array must have dimensions (lmax+1,d) with d=1,4,6 or 9'
        print,'  No file written'
        goto, Exit
    endif
endelse

if (datatype(cl_array) eq 'DOU') then begin
    width = 25                  ;(assuming TFORM= E24.x)
    tform = 'E24.15'
endif else begin
    width = 16                  ;(assuming TFORM= E15.x)
    tform = 'E15.7'
endelse
; create the minimal extension header
ncols = width*nsigs-1 ; # of character columns 
FTCREATE,ncols,nrows,xthdr,tab

; add/update header information
if (do_beam) then begin
    SXADDPAR,xthdr,'EXTNAME','WINDOW FUNCTION',' (beam) window function : B(l) or W(l)'
endif else begin
    SXADDPAR,xthdr,'EXTNAME','POWER SPECTRUM',' Power spectrum : C(l)  '
endelse
SXADDPAR,xthdr,'COMMENT','-----------------------------------------------',before='EXTNAME'
SXADDPAR,xthdr,'COMMENT','-----------------------------------------------',after='EXTNAME'
SXADDPAR,xthdr,'TFIELDS',nsigs,' # of fields in file'
SXADDPAR,xthdr,'NAXIS1',width*nsigs-1,    ' # of characters in a row'
SXADDPAR,xthdr,'NAXIS2',nrows, ' # of rows in file'

if (do_beam) then begin
    ttype_kw = ['TEMPERATURE','GRADIENT','CURL']
    ttype_cm = [' Temperature B(l)',' Gradient (=ELECTRIC) polarisation B(l)',$
                ' Curl (=MAGNETIC) polarisation B(l)']
endif else begin
;    ttype_kw = ['TEMPERATURE','GRADIENT','CURL','G-T','C-T','C-G']
    ttype_kw = ['TEMPERATURE','GRADIENT','CURL', 'TG','TC','GC', 'GT','CT','CG']
    ttype_cm = [' Temperature C(l)',' Gradient (=ELECTRIC) polarisation C(l)',$
                ' Curl (=MAGNETIC) polarisation C(l)', $
;                 ' Gradient-Temperature cross terms', $
;                 ' Curl-Temperature cross terms', $
;                 ' Curl-Gradient terms']
                ' Temperature X Gradient cross terms', $
                ' Temperature X Curl cross terms', $
                ' Gradient X Curl cross terms', $
                ' Gradient1 X Temperature2 cross terms', $
                ' Curl1 X Tempeture2 cross terms', $
                ' Curl1 X Gradient2 cross terms' ]
endelse

for i = 0,nsigs-1 do begin
  keyword = 'TTYPE' + STRTRIM(STRING(i+1),2)
  SXADDPAR,xthdr,keyword,ttype_kw[i],ttype_cm[i]
  keyword = 'TBCOL' + STRTRIM(STRING(i+1),2)
  SXADDPAR,xthdr,keyword,1+width*i,'beginning column of field'
  keyword = 'TFORM' + STRTRIM(STRING(i+1),2)
  SXADDPAR,xthdr,keyword,tform
  if (defined(units)) then begin
      keyword = 'TUNIT' + STRTRIM(STRING(i+1),2)
      SXADDPAR,xthdr,keyword,strtrim(units)
  endif
endfor

lmax = nrows - 1
ft = ['F','T']
polar_flag = (nsigs ge 4) || (do_beam && (nsigs gt 1)) ; 0 or 1
sxaddpar,xthdr,'MAX-LPOL', lmax,          ' Maximum L multipole'
sxaddpar,xthdr,'POLAR',    ft[polar_flag],' Polarisation included (True/False)'
sxaddpar,xthdr,'BCROSS',   ft[nsigs ge 6],' Magnetic cross terms included (True/False)'
sxaddpar,xthdr,'ASYMCL',   ft[nsigs ge 9],' Asymmetric pol cross terms (XY vs YX) included'
nine = ['1','2','3','4','5','6','7','8','9']
; add the user defined header at the end of routine generated one
if (STRMID(xhdr(0),0,1) NE ' ') then begin
    ; remove "trivial" keywords from user defined header
    trivial = ['XTENSION','EXTNAME','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT', $
               'TFIELDS','TTYPE'+nine,'TBCOL'+nine,'TFORM'+nine]
    xhdr1 = xhdr
    sxdelpar,xhdr1,trivial
    iend = WHERE( STRUPCASE(STRMID(xthdr,0,8)) EQ 'END     ', nend)
    if (nend eq 1) then begin
        iend = iend(0)
        xthdr = [xthdr(0:iend-1), xhdr1]
    endif
endif

if (keyword_set(cmbfast) && ~ do_beam) then begin
    sxaddpar, xthdr, 'POLNORM', 'CMBFAST',' CMBFAST (EB) convention for polar. spectra'
endif

; add the data values
; for irow = 0L,LONG(nrows-1) do begin
;   for icol = 0L,LONG(nsigs-1) do begin
;     FTPUT,xthdr,tab,icol+1,irow,cl_array(irow,icol)
;   endfor
; endfor
; -------------------------
; vectorial loop, much faster
irow = lindgen(nrows)
for icol = 0L,LONG(nsigs-1) do begin
    FTPUT,xthdr,tab,icol+1,irow,cl_array[irow,icol]
endfor

; write the header 
WRITEFITS, fitsfile, tab, xthdr, /append
; Exit routine
Exit:
return
end



