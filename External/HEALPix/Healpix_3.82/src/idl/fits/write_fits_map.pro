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
pro write_fits_map, filename, data, info_header, coordsys=coordsys, nested=nested, ring=ring, ordering=ordering, units=units, error=error, help=help
;+
; writes a FITS file with data contained in a BINTABLE extension
;
; CALLING SEQUENCE:
;    WRITE_FITS_MAP, filename, data [, info_header, Coordsys=, Nested=, Ring=,
;    Ordering=, Error=, Help=]
;
; INPUTS:
;    filename = name of the output file 
;
;    data = array containing the data
;
; OPTIONAL INPUT
;    info_header = extra header information to be added to the
;     minimal extension header
;
; OPTIONAL KEYWORDS
;    Ring   = if set, specifies in the fits header that the map is in the RING scheme
;    Nested = if set, specifies in the fits header that the map is in the NESTED scheme
;    Ordering = string can be set either to 'NESTED' or 'RING'
;
;     at most one of Ring, Nested or Ordering should be set
;     and at least one, unless the keyword ORDERING is already set in
;     the fits header info_header
;     the value given by Ring, Nested or Ordering overrides the value
;     present in the FITS header 'info_header'
;   in either cases, the actual ordering of the data is NOT altered
;   (use reorder if you want to do so)
;     the card PIXTYPE = 'HEALPIX ' is also added to the fits header
;
;    Coord = if set to either 'C', 'E' or 'G' specifies that the
;      Healpix coordinate system is respectively Celestial=equatorial, Ecliptic
;      or Galactic
;    Units = units of the data set
;
;   Help= if set, an extensive help (this IDL header) is printed
;
; OPTIONAL OUTPUT
;    Error = takes value 1 on output if error occurs
;
; EXAMPLES
;   writting a fits file file.fits containing a RING ordered healpix map in 
;     is : write_fits_map, 'file.fits', map, /ring
;     or : write_fits_map, 'file.fits', map, ordering = 'ring'
;
;   reading a RING ordered map, reshuffling it to NESTED and writting out the
;   new map is :
;     read_fits_map, 'file_ring.fits', map_ring, h0, h1
;     map_nest = reorder(map_ring,in='ring',out='nest')
;     write_fits_map, 'file_nest.fits', map_nest, h1, /nest
;
; RESTRICTION
;   this routine only writes 1-column binary table, for multi column
;   use WRITE_FITS_SB
;
; Note
; ----
; for consistency with Healpix format the user should 
; give missing pixels
; the value -1.6375e30 in the array data
;
; MODIFICATION HISTORY
;  written by Eric Hivon, IAP, Nov 11 1997
;  writes 1024 entries per row and column (Healpix format)
;   this is much faster anyway                              EH Jan-98
;  March 2 1999, Eric Hivon, Caltech
;   solved problems with file unit numbers
;   added info_header
;  April 1999, EH, added ring, nested and ordering keyword
;  Feb 2000, EH, added test on pixel number, 
;    and healpix information (nside) in file header,
;    replace today() by today_fits()
;  Oct 2001, EH, slightly improved removal of trivial FITS keyword from
;    user supplied header
;  Mar 2008, EH, use fxbwritm for better efficiency for large files
;  Jan 2009: calls init_astrolib
;  Nov 2009: added Error and Help keywords. Slightly faster writting by adapting
;  buffer size
;  2010-01-28: make sure that error is defined and =0 when everything is OK
;  2010-05-11: adds BAD_DATA = -1.6375e30  in FITS header
;
;-

code = 'WRITE_FITS_MAP'
if (keyword_set(help)) then begin
    doc_library,code
    return
endif
init_astrolib


syntax_string = ['Syntax : '+code+', filename, data $',$
'   [, info_header, Coordsys=, /Nested, /Ring, Ordering=, Units=, Error=, Help=]']
if N_params() lt 2 or N_params() gt 3 then begin
    print, syntax_string,form='(a)'
    return
endif

if datatype(filename) ne 'STR' or datatype(data) eq 'STR' then begin
    print, syntax_string,form='(a)'
    print,'   the file name comes first, then the data '
    print,'   file NOT written '
    return
endif

sz = size(data)
dim = sz[0]
npix = n_elements(data)
nside = npix2nside(npix,err=errpix)
if errpix ne 0 then begin
    if (dim gt 1) then begin
        npix    = sz[1]
        nside   = npix2nside(npix,err=errp1)
        if (errp1 ne 0) then begin
            goto, abort
        endif else begin
            print,'multi-dimensional Healpix map'
            if (sz[2] ge 3) then begin
                print,'try WRITE_TQU instead'
            endif else begin
                print,'try WRITE_FITS_SB instead'
            endelse
            goto, crash
        endelse
    endif
    abort:
    print,code+': Non-Healpix data set'
    print,' npix = ',npix,' nside = ',nside
    crash:
    print,' *** file NOT written !! ***'
    return
endif

count = 0
if defined(info_header) then begin
    info_hdr = info_header
    ttype = sxpar(info_hdr,'TTYPE1',count=count)
endif
if (count eq 0) then sxaddpar,info_hdr,'TTYPE1','UNKNOWN1','unknown content'

error = 0
; add ordering information to user supplied extension header
add_ordering_fits,info_hdr, nested=nested, ring=ring, ordering=ordering,error=error
if error ne 0 then return

; add NSIDE information
add_nside_fits,info_hdr,nside=nside,error=error
if error ne 0 then return

; add UNITS information
if defined(units) then begin
    add_units_fits,info_hdr,units=units,error=error
    if error ne 0 then return
endif

; add coordsys information to user supplied extension header
if defined(coordsys) then add_coordsys_fits, info_hdr, coordsys=coordsys

; add BAD_DATA keyword (with !healpix.bad_value value)
add_bad_data_fits, info_hdr

; remove reserved keywords that will be added automatically later on
nine = ['1','2','3','4','5','6','7','8','9']
sxdelpar,info_hdr,['XTENSION','BITPIX','NAXIS','NAXIS1','NAXIS2','PCOUNT','GCOUNT','TFIELDS','TBCOL'+nine,'TFORM'+nine]

; ------- primary unit ----------
; opens the file, write a minimal header and close it
WRITEFITS,filename,0

; update date to Y2K
h0 = HEADFITS(filename)
fdate = today_fits()
SXADDPAR,h0,'DATE',fdate,' Creation date (CCYY-MM-DD) of FITS header'

; opens the file, write the header and the image if any, and close it
WRITEFITS, filename, 0, h0

; -------- extension -----------

data = REFORM(data,N_ELEMENTS(data),/OVERWRITE)
nentry_healpix = 1024L
if (npix MOD nentry_healpix) EQ 0 then begin
    nrows = npix / nentry_healpix 
    nentry = nentry_healpix
endif else begin
    nrows = npix
    nentry = 1L
endelse

; create the minimal extension header
FXBHMAKE,xthdr,nrows

; update the header for 1 column of real*4 with TFORMi = ' 1024E'
FXBADDCOL,col,xthdr,REPLICATE(data(0),nentry)

; merge headers
iend = WHERE( STRUPCASE(STRMID(xthdr,0,8)) EQ 'END     ', nend)
if (nend eq 1) then begin
    iend = iend(0)
    xthdr = [xthdr(0:iend-1), info_hdr]
endif
; reopens the file, goes to the extension and puts the  header there
FXBCREATE, unit, filename, xthdr
naxis1 = sxpar(xthdr,'NAXIS1') ; bytes per row

; writes data in the table
; for row = 1LL, nrows do begin
;     ir = row - 1LL
;     ifirst = ir * nentry
;     ilast  = row * nentry - 1LL
;     FXBWRITE,unit, data(ifirst:ilast), col, row
; endfor

step = 12L ; 12 rows at once
if (nrows gt 24) then step = 24L
if (nrows gt 96) then step = 96L
buffersize = step * naxis1
nstep = (nrows + step - 1L)/step
for is=0LL, nstep-1 do begin
    mystep = step < (nrows - is * step)
    ifirst = (is * step) * nentry
    ilast =  ifirst + mystep * nentry - 1LL
    row = 1L + is*step + [0, mystep-1]
    ;;;print,is, ifirst,ilast ,row
    FXBWRITM,unit, col, data(ifirst:ilast), row=row, buffersize=buffersize
endfor

; closes the file
FXBFINISH,unit

return
end


