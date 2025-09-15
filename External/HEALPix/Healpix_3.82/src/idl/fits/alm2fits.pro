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
PRO ALM2FITS, index, alm_array, fitsfile, HDR = hdr, XHDR = xhdr, HELP=help

;+
; NAME:
;       ALM2FITS
;
; PURPOSE:
;       write to a FITS file a binary table extension containing spherical
;       harmonic coefficients (and optioanl errors) and their index. Writes 
;       header information if required.
;
; CALLING SEQUENCE:
;       ALM2FITS, index, alm_array, fitsfile, [HDR=, XHDR=, HELP= ]
; 
; INPUTS:
;       fitsfile = String containing the name of the file to be written      
;       index    = Long array containing the index for the corresponding 
;                  array of alm coefficients (and erralm if required). The
;                  index {i} is related to {l,m} by the relation
;                  i = l^2 + l + m + 1
;      alm_array = real array of alm coefficients written to the file.
;                  This has dimension (nl,nalm,nsig) -- corresponding to
;                  nl   = number of {l,m} indices
;                  nalm = 2 for real and imaginary parts of alm coefficients;
;                         4 for above plus corresponding error values
;                  nsig = number of signals to be written (1 for any of T E B
;                         or 3 if ALL to be written). Each signal is stored
;                         in a separate extension.
;
; OPTIONAL INPUTS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUT KEYWORDS:
;       HDR      = String array containing the header for the FITS file.
;       XHDR     = String array containing the extension header(s). If 
;                  ALL signals are required, then the three extension 
;                  headers are returned appended into one string array.
;                  NOTE: optional header strings should NOT include the
;                  header keywords explicitly written by this routine.
;       HELP     = If set, this documentation header is shown
;
; EXAMPLE:
;       Write the T spherical harmonic coefficients with corresponding errors
;       stored in the array, alm, into the FITS file, coeffs.fits, as indexed 
;       by the array index, and extension header information contained in the 
;       string array, ext_text
;
;       IDL> alm2fits,index,alm,'coeffs.fits',XHDR=ext_txt
;
;
; PROCEDURES CALLED:
;       FXHMAKE, SXADDPAR, FXWRITE, FXBHMAKE, FXADDPAR, FXBCREATE, 
;       FXBWRITE, FXBFINISH
;       TODAY_FITS()
;
; MODIFICATION HISTORY:
;       May 1999: written by A.J. Banday (MPA)   
;       Feb 2000: replace today() by today_fits(), EH (Caltech)      
;       Aug 2002: added EXTNAME and TTYPE* keywords, EH (IPAC)
;       Aug 2005: replaced FXBWRITE by much faster FXBWRITM
;       Jan 2009: calls init_astrolib
;       Jan 2010: added HELP keyword
;       Aug 2018: user provided HDR is not ignored anymore
;       Mar 2019: fixed bug introduced on Aug 2018
;-

if (keyword_set(help)) then begin
    doc_library,'alm2fits'
    return
endif

if N_params() LT 3 or N_params() gt 5 then begin
  print,'Syntax : ALM2FITS, index, alm_array, fitsfile, [HDR=, XHDR=, HELP=]'
  goto, Exit
endif

; run astrolib routine to set up non-standard system variables
init_astrolib

if n_elements(hdr)  eq 0 then hdr  = ' '
if n_elements(xhdr) eq 0 then xhdr = ' '

; ------- primary unit ----------
; create the primary header information
fxhmake,h0,0,/extend,/date            ; default primary header

; update date to Y2K
fdate = today_fits()
SXADDPAR,h0,'DATE',fdate,' Creation date (CCYY-MM-DD) of FITS header'

; add header information given by user
h1 = h0  ; March 2019
if (STRMID(hdr[0],0,1) NE ' ') then begin
    ; clean user provided header
    myhdr = hdr
    sxdelpar,myhdr,'NAXIS'
    sxdelpar,myhdr,'EXTEND'
    sxdelpar,myhdr,'BITPIX'
    ; remove END from default user
    iend = WHERE( STRUPCASE(STRMID(h1,0,8)) EQ 'END     ', nend)
    if (nend eq 1) then h1 = h1[0:iend[0]-1]
    ; add user provided user, after default one
    h1 = [h1, myhdr]  ; March 2019
endif

; initial write of minimal primary header and null image
;fxwrite,fitsfile,h0 ; fixed on Aug 2018
fxwrite,fitsfile,h1

; -------- extension -----------

; array information
info  = size(alm_array)
type  = info[info[0]+1]
nrows = info[1]
ncols = info[2] + 1 ; # of columns in array, ADD one for output # of columns
if(info(0) gt 2)then nsigs = info(3) else nsigs = 1

nbc = 4 & form = '1E' & fcomm = ' 4-byte float'
if type eq 5 then begin ; double precision
    nbc = 8 & form = '1D' & fcomm = ' 8-byte float'
endif

extname = ['TEMPERATURE','GRAD/ELECTRIC','CURL/MAGNETIC']+' alm'
ttype = ['INDEX','REAL','IMAGINARY','REAL ERROR','IMAGINARY ERROR']
ttype_com = [' index = l*l+l+m+1',' real part of alm',' imaginary par of alm',' error on real part',' error on imaginary part']
; run a loop to create successive binary extension tables: 1 = T; 2 = E; 3 = B
for ix = 0,nsigs-1 do begin

; create the extension header information
    fxbhmake,xthdr,nrows        ; default extension header
    fxaddpar,xthdr,'EXTNAME',extname[ix],' complex alm coefficients'

; add/update header information
    nbytes = 4 + nbc*(ncols-1)
    fxaddpar,xthdr,'NAXIS1  ',nbytes,'Number of 8 bit bytes in each row'
    fxaddpar,xthdr,'NAXIS2  ',nrows,'Number of rows in file'
    fxaddpar,xthdr,'TFIELDS ',ncols,'Number of fields (columns) in the table'
    fxaddpar,xthdr,'TTYPE1  ',ttype[0],ttype_com[0]
    fxaddpar,xthdr,'TFORM1  ','1J        ',' 4-byte integer'
    for i = 1,ncols-1 do begin
        si = STRTRIM(STRING(i+1),2) + '  '
        fxaddpar,xthdr,'TTYPE'+si,ttype[i],ttype_com[i]
        keyword = 'TFORM' + si
        fxaddpar,xthdr,keyword,form, fcomm
    endfor

    index2lm, index, l, m
    lmax = max(l) & mmax = max(m)
    l = 0 & m = 0
    fxaddpar,xthdr,'MAX-LPOL',lmax,' Maximum L multipole order'
    fxaddpar,xthdr,'MAX-MPOL',mmax,' Maximum M multipole degree'

; add the user defined header 
    if (STRMID(xhdr(0),0,1) NE ' ') then begin
        iend = WHERE( STRUPCASE(STRMID(xthdr,0,8)) EQ 'END     ', nend)
        if (nend eq 1) then begin
            iend = iend(0)
            xthdr = [xthdr(0:iend-1), xhdr]
        endif
    endif

; write the binary table extension header
    fxbcreate,lun,fitsfile,xthdr ; add extension header to disk FITS file
    ; write several rows at a time
    nloop = (long(nrows) * ncols) / 1.e4 ; write 10^4 values at a time
    nloop = long(nloop+1) < nrows
    nw = nrows/nloop + 1L
    for i=0L, nloop-1 do begin
        low = long64(i * nw)
        if (low ge nrows) then goto, done
        hi  = long64(low + nw - 1) < (nrows-1)
        case ncols of
            3: fxbwritm, lun, [1,2,3],     long(index[low:hi]),alm_array[low:hi,0,ix],alm_array[low:hi,1,ix],row=[low+1L,hi+1L]
            5: fxbwritm, lun, [1,2,3,4,5], long(index[low:hi]),alm_array[low:hi,0,ix],alm_array[low:hi,1,ix],alm_array[low:hi,2,ix],alm_array[low:hi,3,ix],row=[low+1L,hi+1L]
            else: message,'wrong number of columns'
        endcase
        done:
    endfor
    if (hi ne (nrows-1)) then message,'error in FITS writing'

; free the file
    fxbfinish,lun

endfor                          ; loop over extension tables

; Exit routine
Exit:
return
end



