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
pro gnom2fits, map, file_out, rot = rot_ang, coord = coord, reso = reso, units = units, min=mindata, max=maxdata
;+
; 
;  writeout in a fits file a gnomic projection of a patch of sky
;
;  March 1999, EH caltech, 1.0
;  Feb 2000, EH, replaced today() by today_fits()
;  Jan 2007, EH: header variable (head_out) must be undefined for better results
;-

sz = size(map)
nn = sz(0)

if (nn NE 2) then begin 
    print,'map should be 2-dim in gnom2fits'
    return
endif

if undefined(rot_ang) then rot_ang = [0., 0.]
if undefined(units) then units = ' '
bad_data = -1.6375E30

n1 = sz(1)
n2 = sz(2)

bad = where( abs(map/bad_data - 1.) lt 0.01 or (1 - finite(map)), nbad)
map1 = float(map)
if nbad gt 0 then map1[bad] = bad_data

;mindata = MIN(map,MAX=maxdata)
scoord = decode_coord(coord, flon=flon, flat=flat)
ctype1 = flon+'TAN'
crval1 = rot_ang(0)
crpix1 = (n1+1.)/2.
reso_deg = reso / 60.

ctype2 = flat+'TAN'
crval2 = rot_ang(1)
crpix2 = (n2+1.)/2.

fdate = today_fits()

;head_out =strarr(60)
SXADDPAR,head_out,'SIMPLE','T'
SXADDPAR,head_out,'BITPIX',-32,' IEEE single precision floating point'
SXADDPAR,head_out,'NAXIS',nn
SXADDPAR,head_out,'NAXIS1',n1
SXADDPAR,head_out,'NAXIS2',n2
SXADDPAR,head_out,'DATE',fdate,' Creation date (CCYY-MM-DD) of FITS header'
SXADDPAR,head_out,'BAD_DATA',bad_data,' value for missing data',form='(e15.8)'
SXADDPAR,head_out,'BUNIT',units
;
SXADDPAR,head_out,'CTYPE1',ctype1, ' X-axis '
SXADDPAR,head_out,'CRVAL1',crval1, '    Origin coordinate'
SXADDPAR,head_out,'CRPIX1',crpix1, '    Origin pixel index (1 based)'
SXADDPAR,head_out,'CDELT1',-reso_deg,     '    Degrees/pixel'
;
SXADDPAR,head_out,'CTYPE2',ctype2, ' Y-axis '
SXADDPAR,head_out,'CRVAL2',crval2, '    Origin coordinate'
SXADDPAR,head_out,'CRPIX2',crpix2, '    Origin pixel index (1 based)'
SXADDPAR,head_out,'CDELT2',reso_deg,     '    Degrees/pixel'
;
if defined(mindata) then SXADDPAR,head_out,'DATAMIN',mindata
if defined(maxdata) then SXADDPAR,head_out,'DATAMAX',maxdata

SXADDPAR,head_out,'HISTORY',' '
SXADDPAR,head_out,'HISTORY','Gnomonic projection : '+scoord+' coordinates'
SXADDPAR,head_out,'HISTORY','created by Gnomview '
WRITEFITS, file_out, float(map1), head_out
map1 = 0.

return
end

