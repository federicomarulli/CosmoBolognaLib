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
pro write_fits_sb, filename, prim_st, exten_st, Coordsys=coordsys, Nested=nested, Ring=ring, Ordering=ordering, Partial=Partial_usr, Nside=nside_usr, Extension=extension_id, Nothealpix=nothealpix, help=help, hdr=hdr_u, xthdr=xthdr_u
;+
; writes a FITS file with data contained in a BINTABLE extension
;
; CALLING SEQUENCE:
;    WRITE_FITS_SB, filename, prim_st, exten_st, [Coordsys=, Ring=, Nested=,
;    Ordering=, Partial=, Nside=, Extension=, Nothealpix=, Help=, Xthdr=]
;
; INPUTS:
;    filename = name of the output file 
;
;    prim_st  = structure containing the primary unit
;    with prim_st.hdr = prim_st.(0) : primary header to be added to the 
;       minimal header (can be over-rules with HDR)
;    with prim_st.(1) = primary image (if any)
;
;    ** note ** the healpix standard is to write all the information in the extension units,
;  leaving the primary unit empty. To do so, simply put prim_st = 0, and put the
;  information in exten_st (see example below)
;
;   the data is assumed to represent a full sky data set with 
;   the number of data points npix = 12*nside^2
;   unless
;      A)    1) Partial is set OR the input fits header contains OBJECT =
;                  'PARTIAL'
;          AND 
;            2) the Nside is given a valid value OR the fits header contains
;                    a NSIDE
;      OR
;      B)   the keyword Nothealpix is set
;
;   if the keyword Nothealpix is set, the data set can be arbitrary
;
; OPTIONAL INPUT
;    exten_st = structure containing the extension unit if any
;    with exten_st.hdr = exten_st.(0) = extension header to be added to the
;        minimal header (will be over-ruled by XTHDR)
;    with exten_st.(1)   = 1st column of the BINTABLE
;    with exten_st.(2)   = 2nd column of the BINTABLE if any
;    ...
;
;    by default the name given to each colum (TTYPE*) is the tag name of the structure
;     
; KEYWORD
;    Ring- if set, add 'ORDERING= RING' to the extension fits header
;    Nested- if set, add 'ORDERING= NESTED' to the extension fits header
;    Ordering- if set to the string 'nested' or 'ring', set the
;              keyword 'ORDERING' to the respective value
;
;    one of them has to be present unless the ordering information is already
;    present in the fits header, or unless Nothealpix is set.
;
;    Nside- only useful if the data is not full sky and Partial is set,
;      by default Nside = sqrt(npix/12) where npix is the number of data point
;
;    Help: prints out this header and leaves
;
; OPTIONAL KEYWORD
;   Coord = if set to either 'C', 'E' or 'G' specifies that the
;   Healpix coordinate system is respectively Celestial=equatorial, Ecliptic or
;   Galactic
;   
;   Extension = extension unit in which to put the data (0 based), default=0
;
;   Nothealpix = if set, the data set to be output is not a Healpix map, and the
;   restriction on the number of pixels, Nside, Ordering, Ring, Nested do not apply
;
;   Hdr = FITS Header for primary HDU, will over-rule  prim_st.hdr
;   Xthdr = FITS Header for extension, will over-rule exten_st.hdr
;
; EXAMPLES
;  to write out in 'map.fits' the RING-ordered Healpix maps 
;      signal (real) and n_obs (int) in Galactic coordinates,
;      NB: signal and n_obs should have the same size
;
;  e = create_struct('HDR',' ',  'SIGNAL',signal, 'N_OBS',n_obs) 
;                ; empty header, 1st column = signal, 2nd column = n_obs
;  write_fits_sb, 'map.fits', 0, e, order='ring', coord='g'
;                ; fits file with empty primary unit and binary extension
;
; Note
; ----
; for consistency with Healpix format the user should 
; give missing pixels
; the value !healpix.bad_value = -1.6375e30 in the array data
;
; MODIFICATION HISTORY:
;    written by Eric Hivon, IAP, Nov 11 1997
;    * writes 1024 entries per row and column (Healpix format)
;      this is much faster anyway                              EH Jan-98
;    March 2 1999, Eric Hivon, Caltech
;       solved problems with file unit numbers
;       added the structure interface
;  Feb 2000, EH, added test on pixel number, 
;    and healpix information (nside) in file header,
;    replace today() by today_fits()
;  Sept 2000, EH, corrected various bugs,
;    uses tag names as default column names
;    uses FXBWRITM instead of FXBWRITE
;  Feb 2005, EH, makes sure xhdr is not defined if it is empty
;           added the Nothealpix keyword
;           improved user defined header handling
;  March 2005, EH, extreme speed improvement (x10) by writing only a few rows
;       at a time
;  May 2005, EH, replaces FINDFILE by FILE_TEST
;  Jan 2009: calls init_astrolib
;  Nov 2009: increased buffersize in fxbwritm for slightly faster writing
;  2010-05-11: adds BAD_DATA = -1.6375e30  in FITS header of *Healpix* data sets
;  2012-09-27: added HELP keyword. Work around for STRING(format='(:)') bug in GDL
;  2017-04-25: increased buffersize from 2MB to 1GB, 
;     which makes large difference on some systems (eg, NERSC's edison)
;     but none on others (eg, my mac)
;  2017-10-27: added hdr and xthdr,  use UTC date
;-
;
; NB : do NOT use 'T' or 'F' as the tag names
;-

routine = 'write_fits_sb'
if (keyword_set(help)) then begin
    doc_library,routine
    return
endif

if n_params() ne 3 then begin
    message,' Invalid number of arguments ',/nopref,/info 
    message,'Syntax : WRITE_FITS_SB, filename, prim_st, exten_st, $ ',/nopref,/info,/noname
    message,'    [Coordsys=, Ring=, Nested=, Ordering=, Partial=, Nside=, Extension=, Nothealpix=, Hdr=, Xthdr=, Help=]',/nopref,/noname,/info
    return
endif

init_astrolib

;-------------------------------------------------------------------------------
create = ~keyword_set(extension_id)
do_hpx = ~keyword_set(nothealpix)

if (create) then begin
    hdr = strarr(1)
    image = 0
    IF (N_TAGS(prim_st) GE 1) THEN hdr   = defined(hdr_u) ? hdr_u : prim_st.(0)
    IF (N_TAGS(prim_st) GE 2) THEN begin
        image = prim_st.(1)
        if (n_elements(image) gt 1) then begin
            print,'WARNING : writting data in the primary unit is NOT Healpix standard'
            print,'write_fits_sb, filename, prim_st, exten_st, [order=]'
        endif
    endif
    IF (N_TAGS(prim_st) GE 3) THEN begin
        print,'ERROR : multi column data should be written in the Extension, not the primary units'
        print,'write_fits_sb, filename, prim_st, exten_st, [order=]'
        print,'*** file not written ***'
        message,'Abort'
    endif


; minimal header for this image
    WRITEFITS, filename, 0
    h0 = HEADFITS(filename)
    CHECK_FITS, image, h0, /UPDATE, /FITS

; update date to Y2K
    fdate = today_fits(/utc)
    SXADDPAR,h0,'DATE',fdate,' Creation date (CCYY-MM-DD UTC) of FITS header'

; add header information given by user
    sxdelpar,hdr,['BITPIX','NAXIS','NAXIS1','NAXIS2','DATE','SIMPLE','EXTEND'] ; remove redundant information
    iend = WHERE( STRUPCASE(STRMID(h0,0,8)) EQ 'END     ', nend)
    iend = iend[0]
    if (strtrim(hdr[0],2) NE '') then begin
        h1 = [h0[0:iend-1], hdr,h0[iend]] 
;    print,h1,form='(a80)'
    endif else begin
        h1 = h0
    endelse

; opens the file, write the header and the image if any, and close it
    WRITEFITS, filename, image, h1
    image = 0 & hdr = '' & h0=0 & h1=0

endif
;-------------------------------------------------------------------------------
;------------------------------
; write the extension if any
;------------------------------
IF (N_ELEMENTS(exten_st) EQ 0) THEN message,'empty data set' ; undefined
IF (N_TAGS(exten_st)     EQ 0) THEN message,'invalid data structure' ; not a structure
;xhdr = exten_st.(0)
xhdr = defined(xthdr_u) ? xthdr_u : exten_st.(0)
if (n_elements(xhdr) eq 1) then xhdr = [xhdr] ; make sure it is an array
number = n_tags(exten_st)
tags = tag_names(exten_st)

if (do_hpx) then begin
    nside_fits = sxpar(xhdr,'NSIDE',count=c_nside)
    object_fits = STRTRIM(STRING(sxpar(xhdr,'OBJECT',count=c_object)),2)

    fullsky = 1
    if undefined(nside_usr) then nside_usr = 0
    junk1 = nside2npix(nside_fits,err=err1)
    junk2 = nside2npix(nside_usr,err=err2)
    if ((err1 eq 0 or err2 eq 0) and (object_fits eq 'PARTIAL' or keyword_set(partial_usr))) then begin
        if err2 eq 0 then nside = nside_usr else nside = nside_fits
        fullsky = 0
    endif
endif else begin
    fullsky = 0
endelse

npix = n_elements(exten_st.(1))
if (fullsky) then begin    
    nside = npix2nside(npix,err=errpix)
    if errpix ne 0 then begin
        print,'write_fits_sb: Non-Healpix data set'
        print,' npix = ',npix,' nside = ',nside
        print,' *** file NOT written !! ***'
        message,'Abort'
    endif
endif

; if xhdr is empty make it undefined
if (n_elements(xhdr) eq 1) then begin
    if (strtrim(xhdr[0],2) eq '') then junk=temporary(xhdr)
endif

if (do_hpx) then begin
; add ordering information according to user supplied value
    add_ordering_fits,xhdr, nested=nested, ring=ring, ordering=ordering,error=error
    if error ne 0 then begin
        print,' *** file NOT written '
        message,'Abort'
    endif

; add NSIDE information
    add_nside_fits, xhdr, nside=nside, partial = 1-fullsky, obs_npix = npix

; add BAD_DATA information
    add_bad_data_fits, xhdr
endif

; add coordsys information to user supplied extension header
if defined(coordsys) then add_coordsys_fits, xhdr, coordsys=coordsys

; remove reserved keywords that will be added automatically later on
if (undefined(xhdr)) then sxaddpar,xhdr,'COMMENT','  '
sxdelpar,xhdr,['XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','TFIELDS']
sxdelpar,xhdr,'TFORM'+strtrim(string(indgen(20)+1,form='(i2)'),2)

; Putting many elements per entry can speed up the IO.
; However FXB* does NOT support variable length arrays: the total number of elements for each
; field shoud be an exact multiple of the number of elements per table entry
nentry_healpix = 1024L ; number of elements per entry
if ((npix MOD nentry_healpix) EQ 0 and do_hpx) then begin
    nrows = npix / nentry_healpix 
    nentry = nentry_healpix
endif else begin
    nrows = npix
    nentry = 1L
endelse

; create the minimal extension header
FXBHMAKE,xthdr,nrows

for i=1, number-1 do begin
    data0 = (exten_st.(i))[0]

    ; update the header for 1 column with TFORMi = ' 1024x'
    FXBADDCOL,col,xthdr,REPLICATE(data0,nentry)

    ; put default value for TTYPE*
    si = strtrim(string(i,form='(i2)'),2)
    ttype_in = sxpar(xhdr,'TTYPE'+si,count=cttype_in) ; are TTYPE* given by user ?
    if (cttype_in eq 0) then begin ; if not, put structure tag name instead
        sxaddpar,xthdr,'TTYPE'+si,tags[i]
    endif
endfor

; add the user defined header (xhdr) to xthdr
iend = WHERE( STRUPCASE(STRMID(xthdr,0,8)) EQ 'END     ', nend)
xh_tmp = (nend eq 1) ? xthdr[0:iend[0]-1] : xthdr
for i=0, n_elements(xhdr)-1 do begin
    if (strtrim(xhdr[i]) ne '') then xh_tmp = [xh_tmp, xhdr[i]]
endfor
xthdr = xh_tmp

; if reopening, check that file exists
xtn_id = 0
if defined(extension_id) then begin
    xtn_id = extension_id
    if (file_test(filename) eq 0 and extension_id gt 0) then begin
        print,'********************************************************'
        print,filename+' not found, can not access requested extension'
        message,' Abort'
    endif

    fits_info, filename, /silent, n_ext = n_ext
    if (xtn_id gt n_ext) then begin
        print,'********************************************************'
        print,'Can not access requested extension in '+filename
        message,' Abort'
    endif

    if (xtn_id lt n_ext) then begin
        print,'********************************************************'
        print,'Can not modify existing extension in '+filename
        message,' Abort'
    endif
endif
; reopens the file, goes to the extension and puts the  header there
FXBCREATE, unit, filename, xthdr, xtn_id
naxis1 = sxpar(xthdr,'NAXIS1')

;nloop = (long(nrows) * nentry * (number-1L)) / 1.e4 ; write 10^4 values at a time
;buffersize = 32768L*64 ; number of bytes per chunk  (2MB)
buffersize = 32768LL^2 ; number of bytes per chunk  (1GB)
;;bytes_per_row = nentry * (number-1L) * 4 ; nuber of bytes per row in FITS file = NAXIS1
bytes_per_row = naxis1
nloop = (long64(nrows) * bytes_per_row) / buffersize  ; number of chunks of size buffersize
nloop = long(nloop+1) < nrows
nw = nrows/nloop + 1L ; number of rows in each chunk
buffersize = nw * bytes_per_row ; final size of chunk
ip1   = indgen(number-1)+1 ; from 1 to number-1
;print,nloop,nw,nw*nloop-nrows, buffersize

stcol = '['+string(ip1,form='(50(i2, :, ","))')+']'
stbuf = ',buffersize='+strtrim(buffersize,2)
for i=0LL, nloop-1LL do begin
    low = long64(i * nw)
    if (low ge nrows) then goto, done
    hi  = long64(low + nw - 1) < (nrows-1)
; uses FXBWRITM rather than FXBWRITE, and write a few rows at a time
    srange = strcompress(string(low*nentry,hi*nentry+nentry-1,form='("[",i15,":",i15,"]")'),/remove_all)
    if (is_gdl()) then begin
        starg = ''
        for ip=1,number-1 do starg += string(ip, form='(",exten_st.(",i0,")'+srange+'")')
        starg = strcompress(starg, /remove_all)
    endif else begin
        starg =  strcompress(string(ip1,form='(50(:, ",exten_st.(",i2,")'+srange+'"))'),/remove_all)
    endelse
    rows =   strcompress(string(low+1L,hi+1L,form='(",row=[",i12,",",i12,"]")'),/remove_all)
    command = 'FXBWRITM, unit, '+stcol+starg+rows+stbuf
    junk = execute(command)
    done:
endfor



if (hi ne (nrows-1)) then message,'error in FITS writing'

; closes the file
FXBFINISH,unit

return
end

