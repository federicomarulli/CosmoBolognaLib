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
pro write_kml_file, img_file, kml_file, box=box, rotation=rotation, title=title, silent=silent, subtitle=subtitle

; this default box works for an image with RA increasing towards the left
if undefined(box) then box ={north:90.d0,south:-90.d0,east:180.d0,west:-180.d0}

if undefined(rotation) then rotation = 0.

;;;help,img_file,kml_file
openw,lun,kml_file,/get_lun

printf,lun,'<?xml version="1.0" encoding="UTF-8"?>'
printf,lun,'<kml xmlns="http://earth.google.com/kml/2.1" hint="target=sky">'
printf,lun,'<!--'
printf,lun,'      To open this KML Folder under Google Earth: do'
printf,lun,'           File > Open and select this file.'
;;;;printf,lun,'       It will automatically switch to Google Sky.'
printf,lun,'     '+systime()
printf,lun,' -->'
printf,lun,'  <Folder>'
printf,lun,'    <name>Healpix sky data</name>'
printf,lun,'    <description>Healpix data in Google format</description>'
printf,lun,'    <GroundOverlay>'
;;printf,lun,'      <altitude>0</altitude>'
if defined(title)    then printf,lun,'      <name>'+title+'</name>'
if defined(subtitle) then printf,lun,'      <description>'+subtitle+'</description>'
printf,lun,'      <Icon>'
printf,lun,'        <href>'+img_file+'</href>'
printf,lun,'      </Icon>'
printf,lun,'      <LatLonBox>'
printf,lun,'      <!-- Swapping North and South (or East and West) numbers will flip the image accordingly-->'
printf,lun,'        <north>'+strtrim(box.north,2)+'</north>'
printf,lun,'        <south>'+strtrim(box.south,2)+'</south>'
printf,lun,'        <east>' +strtrim(box.east,2) +'</east>'
printf,lun,'        <west>' +strtrim(box.west,2) +'</west>'
printf,lun,'        <rotation>'+strtrim(rotation,2)+'</rotation>'
printf,lun,'      </LatLonBox>'
printf,lun,'    </GroundOverlay>'
printf,lun,'  </Folder>'
printf,lun,'</kml>'
free_lun,lun

if (~ keyword_set(silent)) then begin
    print,'KML file created in '+kml_file
    print,'Image in '+IMG_file
endif

return
end
;--------------------------------------------------
pro hpx2gs, file_in, select_in, $
                     ;gif = gif, $
                     coord_in = coord_in, $
                     help = help, $
                     kml = kml, $
                     png = png, $
                     _extra = other_kwds, $
                     reso_arcmin = reso_arcmin, $
                     silent = silent, $
                     subtitle = subtitle, $
                     titleplot = titleplot
;+
; NAME:
;    HPX2GS
;
; PURPOSE:
;    Turns a Healpix data set into a cartesian PNG file, 
;    and a KML folder file, usuable by Google Sky, Google Sky and Oculus VR headsets
;
; CATEGORY:
;
; CALLING SEQUENCE:
;   Hpx2gs, file [, select, COORD_IN=, /HELP, KML=, PNG=, RESO_ARCMIN=, SUBTITLE=,
;   TITLEPLOT=, 
;   + all remaining CARTVIEW's keywords 
;    (ASINH=, COLT=, FACTOR=, FLIP=, GRATICULE=, HBOUND=, HIST_EQUAL=, 
;    IGRATICULE=, LOG=, MAX=, MIN=, NESTED=, NO_DIPOLE = , NO_MONOPLE=, 
;    OFFSET=, ONLINE=, OUTLINE=, POLARIZATION=, PREVIEW=, QUADCUBE= ,SAVE=, SILENT=)
;
;    (*except*: CHARSIZE, CROP, FITS, GIF, HXSIZE, NOBAR, NOLABELS,
;    NOPOSITION, PS, PXSIZE, PYSIZE, ROT, TRANSPARENT, UNITS, XPOS, YPOS, WINDOW) ]
;
;
; INPUTS:
; 	File = 
;          by default,           name of a FITS file containing 
;               the healpix map in an extension or in the image field
;          or    name of a variable containing the healpix map
;          if Save is set   :    name of an IDL saveset file containing
;               the healpix map stored under the variable  data
;
;
; OPTIONAL INPUTS:
;       Select =  if the file read is a multi column BIN table, Select indicates
;                 which column is to be plotted (the default is to plot the
;                 first column containing some signal, as opposed to pixel index)
;               can be either a name : value given in TTYPEi of the FITS file
;                        NOT case sensitive and
;                        can be truncated, 
;                        (only letters, digits and underscore are valid)
;               or an integer        : number i of the column
;                            containing the data, starting with 1
;
; KEYWORD PARAMETERS:
;
;   COORD_IN= 1-character scalar, describing the input data coordinate system:
;                either 'C' or 'Q' : Celestial2000 = eQuatorial,
;                       'E'        : Ecliptic,
;                       'G'        : Galactic 
;             If set, it will over-ride the coordinates read from the FITS file header (when
;             applicable).  In absence of information, the input coordinates is
;             assumed to be celestial.
;             The data will be rotated so that the output coordinates are Celestial, as expected by GoogleSky
;
;   /HELP:  prints out this documentation
;
;   KML= filename of the KML file to be created (if the .kml suffix is missing,
;     it will be added automatically). Used only by Google Earth and Google Sky.
;     Default: hpx2googlesky.kml
;
;   PNG= filename of the PNG file to be created. Only to be used if you want the
;     filename to be different from the default (ie same as KML file, with a .png suffix instead
;     of .kml)
; 
;   Note: Pixels missing or unobserved in the input data will be
;   totally 'transparent' in the output file.
;
;   RESO_ARCMIN = pixel angular size in arcmin (at the equator) of the cartesian
;     map generated (default = 30.)
;
;   SUBTITLE = information on the data, will appear in KML file GroundOverlay
;     description field
;
;   TITLEPLOT = information on the data, will appear in KML file GroundOverlay
;     name field
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
;   creates a KML file and a PNG file
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
;   2007-??
;   2008-02-27: documentation, added COORD_IN
;-


;TODO : 
;       check orientation (Google coordinates: Celestial=eQuatorial)
;
;
; solve_kw_conflict,'png', 'gif',       k1=png,    k2=gif,   kout=img_file,
; /defined

routine = 'HPX2GS'
syntax = [routine+', File [, Select, COORD_IN=, /HELP, KML=, PNG=, RESO_ARCMIN=, $',$
          '         SUBTITLE=, TITLEPLOT=, $',$
          '       <+ most of CARTVIEW''s keywords> ]']

if (keyword_set(help)) then begin
    doc_library,routine
    return
endif

if (n_params() eq 0) then begin
    print,syntax,format='(a)'
    print,'Type: '+routine+', /help'
    print,'   for more help.'
    return
endif
;-----------------

if undefined(reso_arcmin) then reso_arcmin = 30.
pxsize = 360.*60./reso_arcmin
pysize = 180.*60./reso_arcmin

default_coord = size(file_in,/tname) eq 'STRING' ? ['C'] : ['C','C']
coord = (defined(coord_in)) ? [coord_in,'C'] : default_coord

; Do not allow cartview keywords that are set internally by this routine
if defined(other_kwds) then begin
    forbidden = ['CROP','PXSIZE','PYSIZE','ROT','TRANSPARENT','WINDOW','GIF']
    useless =   ['CHARSIZE','FITS','HXSIZE','NOBAR','NOLABELS','NOPOSITION' $
                 , 'PS','UNITS','XPOS','YPOS']
    provided = tag_names(other_kwds)
    for i=0, n_elements(provided)-1 do begin
        jforb = index_word(forbidden, provided(i))
        jusel = index_word(useless, provided(i))
        if jforb ge 0 then begin
            message,'Keyword '+provided(i)+' ('+forbidden(jforb)+')'+' not allowed.'
        endif
        if jusel ge 0 then begin
            message,/info,'Keyword '+provided(i)+' ('+useless(jusel)+')'+' will be ignored.'
        endif
    endfor
endif

; define KML and PNG filenames
kml_file = (size(kml,/tname) eq 'STRING') ? kml : 'hpx2googlesky.kml'
iend = strpos(strlowcase(kml_file),'.kml')
if (iend lt 0) then kml_file = kml_file + '.kml' ; add kml suffix

iend = strpos(strlowcase(kml_file),'.kml')
def_img_file = (iend gt 0) ? strmid(kml_file,0,iend) :  kml_file
img_file = (size(png,/tname) eq 'STRING') ? png : def_img_file+'.png'
iend = strpos(strlowcase(img_file),'.png')
if (iend lt 0) then img_file = img_file + '.png' ; add png suffix

; create a PNG file with cartesian projection of data
cartview, file_in, select_in, $
          /silent, $
          /CROP, $ ; crop image
          COORD = coord, $
          PXSIZE = pxsize, $
          PYSIZE = pysize, $
          ROT=[180,0,0], $ ; rotation constrained by Google Sky
          /TRANSPARENT, $  ; transparent missing pixels
          WINDOW = -1, $ ; virtual window
;;;          gif = gif, $
          png = img_file, $
          reso_arcmin = reso_arcmin, $
          _strict_extra = other_kwds


; write KML file
write_kml_file, img_file, kml_file, title=titleplot, subtitle=subtitle, silent=silent

return
end
