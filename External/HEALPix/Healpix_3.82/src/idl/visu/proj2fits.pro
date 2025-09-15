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
pro proj2fits, map, FITSFile, $
               coord = coord, $
               flip  = flip, $
               half_sky  = half_sky, $
               max   = maxdata, $
               min   = mindata, $
               projection  = projection, $
               reso_arcmin = reso_arcmin, $
               rot   = rot_ang, $
               units = units
;+
; 
;  writeout in a FITS file a projection of a (patch of) sky
;
;  requires: SXADDPAR, MAKE_ASTR, PUTAST, WRITEFITS
;
;  March 1999, EH caltech, 1.0: gnom2fits
;  Feb 2000, EH, replaced today() by today_fits()
;  Jan 2007, EH: header variable (head_out) must be undefined for better results
;  May 2010, EH: turned into proj2fits. Uses WCS routines
;  Jul 2015, EH: forced EQUINOX value to 2000.0
;  Aug 2015, EH: make sure putast.pro writes LONPOLE+LATPOLE or PV1_3+PV1_4 in output FITS header
;-

routine = 'proj2fits'
syntax = [routine+', map, fitsfile, $', '  [COORD=, /FLIP, /HALF_SKY, MAX=, MIN=, PROJECTION=, RESO_ARCMIN=, ROT=, UNITS=]']

if (n_params() ne 2) then begin
    print,syntax,form='(a)'
    return
endif


; check map size
sz = size(map)
nn = sz[0]
if (nn NE 2) then begin 
    message,/info,'map should be 2-dim for FITS output.'
    return
endif
n1 = sz[1]
n2 = sz[2]

; flag bad pixels
init_healpix
bad_data = !healpix.bad_value
bad = where( abs(map/bad_data - 1.) lt 0.01 or (~finite(map)), nbad)
case size(/tname, map) of
    'DOUBLE': begin
        map1 = map
        mynan = !values.d_nan
        bitpix = -64
        bit_comment = ' IEEE double precision floating point'
    end
    else: begin
        map1 = float(map)
        mynan = !values.f_nan
        bitpix = -32
        bit_comment = ' IEEE single precision floating point'
    end
endcase
if (nbad gt 0) then map1[bad] = mynan

; units
bunit = defined(units) ? units : ' '

; coordinates
scoord = decode_coord(coord, flon=flon, flat=flat)

; identify projection
identify_projection, projtype, projection=projection

; check rotation vector
my_rot_ang = defined(rot_ang) ? rot_ang : [0.d0, 0.d0, 0.d0]

if (projtype eq 1 || projtype eq 3) then begin
    ; third Euler angle meaningless at pole for cylindrical projections
    if (my_rot_ang[1] eq  90.) then my_rot_ang=[my_rot_ang[0]+my_rot_ang[2], my_rot_ang[1], 0.]
    if (my_rot_ang[1] eq -90.) then my_rot_ang=[my_rot_ang[0]-my_rot_ang[2], my_rot_ang[1], 0.]
endif
DDtoR = !DPI / 180.d0
ang_rad = my_rot_ang * DDtoR

; default
longpole = 180.d0
latpole  = 0.0d0
; polar projections (TAN/gnom, SIN/orth, ARC/azeq);   Native(0,0) = Pole
longpole1 = 180.d0 - my_rot_ang[2]
; cylindrical projections (MOL, CAR);  Native(0,0) = Equator
longpole2 = atan(cos(ang_rad[1])*sin(ang_rad[2]), sin(ang_rad[1])) / DDtoR
latpole2  = asin(cos(ang_rad[1])*cos(ang_rad[2])) / DDtoR

; print,projtype
case projtype of
    1: begin
        ctype = 'MOL'
        longpole = longpole2
        latpole  = latpole2
        proj_small = 'mollweide'
    end
    2: begin
        ctype = 'TAN'
        longpole = longpole1
        proj_small = 'gnomic'
    end
    3: begin
        ctype = 'CAR'
        longpole = longpole2
        latpole  = latpole2
        proj_small = 'cartesian'
    end
    4: begin
        ctype = 'SIN'
        longpole = longpole1
        proj_small = 'orthographic'
        if (~keyword_set(half_sky)) then begin
            message_patch,'FITS output of orthographic map only possible with HALF_SKY',level=-2
        endif
    end
    6: begin
        ctype = 'ARC'
        longpole = longpole1
        proj_small = 'azimequid'
    end
    else: begin
        message,'Unknown projection type'
    end
endcase
;;help, longpole1, longpole
ctype1 = [flon, flat]+ctype
OutFile = (DATATYPE(FITSFile) ne 'STR') ? 'plot_'+proj_small+'.fits' : FITSFile

; date
fdate = today_fits()


; standard FITS keywords
SXADDPAR,head_out,'SIMPLE', 'T'
SXADDPAR,head_out,'BITPIX', bitpix, bit_comment
SXADDPAR,head_out,'NAXIS',  nn
SXADDPAR,head_out,'NAXIS1', n1
SXADDPAR,head_out,'NAXIS2', n2
SXADDPAR,head_out,'DATE',   fdate,' Creation date (CCYY-MM-DD) of FITS header'
;;SXADDPAR,head_out,'BAD_DATA',bad_data,' value for missing data',form='(e15.8)'
SXADDPAR,head_out,'BUNIT',  bunit
sxaddpar,head_out,'COMMENT','----------------------------------------'
sxaddpar,head_out,'COMMENT','    WCS keywords                        '
sxaddpar,head_out,'COMMENT','----------------------------------------'
if (ctype eq 'TAN') then begin ; added 2015-08-05
; PUTAST.pro features (as of Oct 6, 2014):
;  - gnomic projection: add LONPOLE and LATPOLE slots that will be updated by the routine
;  - other projections: PV1_3 and PV1_4 will be written instead by the routine
    sxaddpar,head_Out,'LONPOLE',longpole,' [Deg] Native longitude of '+scoord+' pole'
    sxaddpar,head_Out,'LATPOLE',latpole ,' [Deg] Native latitude  of '+scoord+' pole'
endif

delt_deg = keyword_set(flip) ? [1,1] : [-1,1]
; WCS keywords
; print, flon, flat
; wcs_check_ctype, ctype1, _projtype, _coord_type

equinox = 2000.0d ; added 2015-07-02
;;radecsys = 'FK5'
make_astr, astr, $
           ctype = ctype1 $     ; projection type
           , delt  = delt_deg*(reso_arcmin/60.d0) $ ; resolution (Degrees/pixel)
           , crpix = [n1+1.,n2+1.]*0.5d0 $ ; index of center pixel (1-based)
           , crval = my_rot_ang[0:1] $ ; long-lat location of center (Degrees)
           , longpole = longpole $
           , latpole  = latpole $
           , equinox  = equinox
           ;;;;;; , radecsys = radecsys $

putast, head_out, astr, cd_type=0

sxaddpar,head_out,'COMMENT','----------------------------------------',after='HISTORY'
; data related keywords
if defined(mindata) then SXADDPAR,head_out,'DATAMIN',mindata, ' Data minimum value'
if defined(maxdata) then SXADDPAR,head_out,'DATAMAX',maxdata, ' Data maximum value'

SXADDPAR,head_out,'HISTORY',' '
SXADDPAR,head_out,'HISTORY', proj_small+' projection, '+scoord+' coordinates'
;;;;;SXADDPAR,head_out,'HISTORY','created by Gnomview '
WRITEFITS, OutFile, map1, head_out
map1 = 0.

return
end

