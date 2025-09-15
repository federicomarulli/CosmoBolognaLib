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
pro write_kml_file, img_file, kml_file, box=box, rotation=rotation, title=title, subtitle=subtitle

; this default box works for an image with RA increasing towards the left
if undefined(box) then box ={north:90.d0,south:-90.d0,east:180.d0,west:-180.d0}

if undefined(rotation) then rotation = 0.

;;;help,img_file,kml_file
openw,lun,kml_file,/get_lun

printf,lun,'<?xml version="1.0" encoding="UTF-8"?>'
printf,lun,'<kml xmlns="http://earth.google.com/kml/2.1" hint="target=sky">'
printf,lun,'<!--'
printf,lun,'      To open this KML Folder under Google Earth: do'
printf,lun,'           File > Open and select this file.
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
printf,lun,'      <!-- Swapping North and South (or East and West) numbers will flip the image accordingly-->
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

print,'KML file created in '+kml_file

return
end
;--------------------------------------------------
pro cartgoogle, file_in, select_in, $
                gif = gif, $
                kml = kml, $
                png = png, $
                _extra = other_kwds, $
                reso_arcmin = reso_arcmin, $
                subtitle = subtitle, $
                titleplot = titleplot
;+
; NAME:
;
;
;
; PURPOSE:
;    Turns a Healpix data set into a PNG or GIF file, and a KML folder file, usuable by Google Sky
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;   Cartgoogle, data [, select, GIF=, KML=, PNG=, RESO_ARCMIN=, SUBTITLE=,
;   TITLEPLOT=, 
;   + all CARTVIEW's keywords (COLT=, FACTOR=, GRATICULE=, IGRATICULE=, HBOUND=,
;   HIST_EQUAL=, LOG=, ...,MIN=, MAX=,  OUTLINE=, POLARIZATION=, PREVIEW= ....
;
;       *except*: CROP, PXSIZE, PYSIZE, ROT, TRANSPARENT, WINDOW) ]
;
;
; INPUTS:
;   data
;
;
; OPTIONAL INPUTS:
;   select
;
;
; KEYWORD PARAMETERS:
;   GIF=
;
;   /KML= if set, a KML file is created with the information necessary to use
;   the GIF or PNG file.
;
;   PNG= png file. Pixels missing or unobserved in the input data will be
;   totally 'transparent' in the output file.
;
;   RESO_ARCMIN = pixel angular size in arcmin (at the equator) of the cartesian
;   map generated
;
;   SUBTITLE = information on the data, will appear in KML file GroundOverlay
;   description field
;
;   TITLEPLOT = information on the data, will appear in KML file GroundOverlay
;   name field
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-
;    COLT=
;    COORD=
;    FACTOR=
;    FLIP=
;    GRATICULE=
;    HBOUND=
;    HELP=
;    HIST_EQUAL=
;    IGRATICULE=
;    LOG=
;    MAX=
;    MIN=
;    NESTED=
;    NO_DIPOLE = 
;    NO_MONOPLE=
;    OFFSET=
;    ONLINE=
;    OUTLINE=
;    POLARIZATION=
;    PREVIEW=
;    QUADCUBE= ,
;    SAVE=


;TODO : documentation (idl, latex)
;       check orientation (Google coordinates: Celestial)
;       transparency: Done
;
solve_kw_conflict,'png', 'gif',       k1=png,    k2=gif,   kout=img_file, /defined

if undefined(reso_arcmin) then reso_arcmin = 30.
pxsize = 360.*60./reso_arcmin
pysize = 180.*60./reso_arcmin

; Do not allow cartview keywords that are set internally by this routine
if defined(other_kwds) then begin
    forbidden = ['CROP','PXSIZE','PYSIZE','ROT','TRANSPARENT','WINDOW']
    useless =   ['CHARSIZE','FITS','HXSIZE','NOBAR','NOLABELS','NOPOSITION' &
    'PS','UNITS','XPOS','YPOS']
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

cartview, file_in, select_in, $
          /CROP, $ ; crop image
          PXSIZE = pxsize, $
          PYSIZE = pysize, $
          ROT=[180,0,0], $ ; rotation constrained by Google Sky
          /TRANSPARENT, $  ; transparent missing pixels
          WINDOW = -1, $ ; virtual window
          gif = gif, $
          png = png, $
          reso_arcmin = reso_arcmin, $
          _strict_extra = other_kwds

if keyword_set(kml) then begin
    if keyword_set(png) then suff = '.png'
    if keyword_set(gif) then suff = '.gif'
    iend = strpos(strlowcase(img_file),suff)
    if (iend gt 0) then kmlfile = strmid(img_file,0,iend) +'.kml' else kmlfile = img_file +'.kml'
    write_kml_file, img_file, kmlfile, title=titleplot, subtitle=subtitle
endif

return
end
