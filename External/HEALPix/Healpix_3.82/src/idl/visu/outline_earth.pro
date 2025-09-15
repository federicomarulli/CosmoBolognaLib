;+
; NAME:
;         outline_earth
;
; PURPOSE:
;        Creates or expands a HEALPix-compatible IDL structure containing
;        geographical data such as continental and/or countries boundaries,
;        (inner and outer) coastlines and rivers.
;        The data are in Celestial = eQuatorial coordinates
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;      geo = outline_earth([geo_in],
;                        [BEHAVIOR=, COASTS=, COLOR=, CONTINENTS=, COUNTRIES=, 
;                         HELP=, HIRES=, LINESTYLE=, RIVERS=, THICK=])
;
; INPUTS:
;      
;
; OPTIONAL INPUTS:
;       geo_in : structure to be expanded into geo.
;                If absent, a new structure is created in geo
;
;
; KEYWORD PARAMETERS:
;
;       BEHAVIOR= either 'GDL', 'IDL' or absent: 
;            changes the origin of the Earth data (see below)
;
;       COASTS= if set, adds coastlines, islands, and lakes information in geo
;                  with the features COLOR, LINESTYLE and THICK
;
;       COLOR= color index to be given to geographical data in final plot
;                 in [0, 255]
;
;       CONTINENTS= if set, adds continental boundaries.
;               Can be combined with COASTS, COUNTRIES and RIVERS
;
;       COUNTRIES= if set, adds political boundaries (as of 1993)
;
;       HELP= if set, this documentation is printed out and the routine exits
;
;       HIRES= if set >0 use high resolution information (if it is available) 
;              instead of default low resolution one.
;              Either 0 or 1 (and more) in IDL mode;
;              in {0,1,2,3,4} in GDL mode.
;             Beware that Hires data are voluminous, longer to process
;              and will result in large PNG or PDF files when plotted with eg mollview
;
;       LINESTYLE=  linestyle of the geographical data in mollview
;
;       RIVERS= if set, adds rivers
;               Can be combined with COASTS, CONTINENTS and COUNTRIES
;
;       THICK= thickness of the geographical lines (default=1.0)
;
; OUTPUTS:
;       geo: IDL structure containing geographical information, which
;         then can be used with eg 'mollview, outline = geo'
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;       None
;
; SIDE EFFECTS:
;       Read geographical data
;
; RESTRICTIONS:
;
;        If run under IDL (or if BEHAVIOR='IDL') the data are read from
;         ${IDL_DIR}/resource/maps/*/*.*
;         with 'low' (HIRES=0) and if available 'high' (HIRES>=1) resolution.
;
;        If run under GDL (or if BEHAVIOR='GDL') the GSHHS data 
;         (version 2.2 or more, available at 
;         https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-bin-*.zip)
;         are read from ${GSHHS_DATA_DIR}/*.b 
;         with HIRES=0,1,2,3,4 standing for 'coarse', 'low', 'intermediate', 
;         'high' and 'full' resolution respectively.
;         In that case, Coasts and Continents are degenerate.
;
;        Under FL, BEHAVIOR must be set to either 'IDL' or 'GDL', and
;         the corresponding data must be available.
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;      geo = outline_earth(/continents,/hires,thick=0.8) ; high res continents with thick lines
;      geo = outline_earth(geo, /rivers,      thick=0.1,col=20) ; add low res rivers with thiner lines and a different color
;      mollview, replicate(!healpix.bad_value,12), outline=geo, /flip, coord=['C','C']
;
; MODIFICATION HISTORY:
;      2018-12: creation
;      2019-01: add support for THICKness,
;               faster in GDL and FL by using HASH instead of STRUCTURE
;               
;-
;=====================================
function hash_support
; ---------------------------------------
; hash(es) are supported in IDL >= 8.0 and in GDL and FL.
; However they improve performance over structure(s)
; - significantly for GDL and FL,
; - barely for IDL >= 8.5
; - negatively for IDL 8.0 to 8.4
; ---------------------------------------
if is_idl() then begin
    dh = 0
    ;if (float(!version.release) gt 8.0) then dh = 1
endif else begin
    dh = 1
endelse
;print,'Hash support: ',dh
return, dh
end
; ---------------------------------------
function merge_structures, st1, st2, tag=tag
; ---------------------------------------
;  merge 2 structures into 1 with an optional tag for the second one
;  if the first structure does not exist, a new strucure is created
;  
; 2018-12-05: creation
; ---------------------------------------
do_tag = defined(tag)
if size(st1,/tname) eq 'STRUCT' then begin ; increment existing structure
    stout = do_tag ? create_struct(st1, tag, st2) : create_struct(st1, st2)
endif else begin ; create new structure
    stout = do_tag ? create_struct(     tag, st2) : st2
endelse

return, stout
end
; -----------------------------------------------------------------------------
pro hpx_c_gdl, fnames, name, hires, stc=stc, color=color, linestyle=linestyle, thick=thick
; ---------------------------------------
; read and parse GSHHS data for GDL
; adapted from D. W. Fanning's http://www.idlcoyote.com/programs/cgmap_gshhs.pro
; ---------------------------------------

dir = getenv('GSHHS_DATA_DIR')
print,'GSHHS_DATA_DIR = '+dir

; gshhs_?.? : continents + lakes + rivers (heavy)
; wdb_borders_?.?: countries    (light)
; wdb_rivers_?.?: rivers

filename = filepath(fnames[hires], root_dir=dir)
;print,filename
IF ~File_Test(filename, /Read) THEN Message, 'Cannot find or read the file ' + filename + '.'


;   ; Define the polygon header. This is for versions of the GSHHS software of 1.3 and earlier.
;   header = { id: 0L, $        ; A unique polygon ID number, starting at 0.
;              npoints: 0L, $   ; The number of points in this polygon.
;              polygonLevel: 0L, $ ; 1 land, 2 lake, 3 island-in-lake, 4 pond-in-island.
;              west: 0L, $      ; West extent of polygon boundary in micro-degrees.
;              east: 0L, $      ; East extent of polygon boundary in micro-degrees.
;              south: 0L, $     ; South extent of polygon boundary in micro-degrees.
;              north: 0L, $     ; North extent of polygon boundary in micro-degrees.
;              area: 0L, $      ; The area of polygon in 1/10 km^2.
;              version: 0L, $   ; Polygon version, always set to 3 in this version.
;              greenwich: 0S, $ ; Set to 1 if Greenwich median is crossed by polygon.
;              source: 0S }     ; Database source: 0 WDB, 1 WVS.

   ; Define the polygon header, for GSHHS software 1.4 through 1.11, which uses a 40 byte
   ; header structure. For example, gshhs_i.b from the gshhs_1.10.zip file.
;   header = { id: 0L, $        ; A unique polygon ID number, starting at 0.
;              npoints: 0L, $   ; The number of points in this polygon.
;              flag: 0L, $      ; Contains polygonlevel, version, greenwich, and source values.
;              west: 0L, $      ; West extent of polygon boundary in micro-degrees.
;              east: 0L, $      ; East extent of polygon boundary in micro-degrees.
;              south: 0L, $     ; South extent of polygon boundary in micro-degrees.
;              north: 0L, $     ; North extent of polygon boundary in micro-degrees.
;              area: 0L, $      ; Database source: 0 WDB, 1 WVS.
;              junk:bytarr(8)}  ; Eight bytes of junk to pad header.     

   ; Define the polygon header, for GSHHS software 1.4 through 1.11, which uses a 32 byte
   ; header structure. For example, gshhs_h.b from the gshhs_1.11.zip.
;   header = { id: 0L, $        ; A unique polygon ID number, starting at 0.
;              npoints: 0L, $   ; The number of points in this polygon.
;              flag: 0L, $      ; Contains polygonlevel, version, greenwich, and source values.
;              west: 0L, $      ; West extent of polygon boundary in micro-degrees.
;              east: 0L, $      ; East extent of polygon boundary in micro-degrees.
;              south: 0L, $     ; South extent of polygon boundary in micro-degrees.
;              north: 0L, $     ; North extent of polygon boundary in micro-degrees.
;              area: 0L}        ; Database source: 0 WDB, 1 WVS.
              
   ; Define the polygon header, for GSHHS software 2.0, which uses a 44 byte
   ; header structure. For example, gshhs_h.b from the gshhs_2.0.zip.
;   header = { id: 0L, $        ; A unique polygon ID number, starting at 0.
;              npoints: 0L, $   ; The number of points in this polygon.
;              flag: 0L, $      ; Contains polygon level, version, greenwich, source, and river values.
;              west: 0L, $      ; West extent of polygon boundary in micro-degrees.
;              east: 0L, $      ; East extent of polygon boundary in micro-degrees.
;              south: 0L, $     ; South extent of polygon boundary in micro-degrees.
;              north: 0L, $     ; North extent of polygon boundary in micro-degrees.
;              area: 0L, $      ; Area of polygon in 1/10 km^2.
;              area_full: 0L, $ ; Area of origiinal full-resolution polygon in 1/10 km^2.
;              container: 0L, $ ; ID of container polygon that encloses this polygon (-1 if "none").
;              ancestor: 0L }   ; ID of ancestor polygon in the full resolution set that was the source
;                               ; of this polygon (-1 of "none").
;
;  I don't have a polygon header for version 2.1 of the GSHHS software but I believe the 2.0 header will
;  work. This version of GSHHS was pretty buggy and not used very long. Better to upgrade to 2.2.
;  
   ; Define the polygon header, for GSHHS software 2.2, which uses a 44 byte
   ; header structure. For example, gshhs_h.b from the gshhg_2.2.zip.
header = { id: 0L, $            ; A unique polygon ID number, starting at 0.
           npoints: 0L, $       ; The number of points in this polygon.
           flag: 0L, $          ; Bytes defined as:
           ;    1st byte:    level = flag & 255: Values: 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake
              ;    2nd byte:    version = (flag >> 8) & 255: Values: Should be 9 for GSHHS release 9
              ;    3rd byte:    greenwich = (flag >> 16) & 3: Values: 0 if Greenwich nor Dateline are crossed,
              ;                 1 if Greenwich is crossed, 2 if Dateline is crossed, 3 if both is crossed.
              ;    4th byte:    source = (flag >> 24) & 1: Values: 0 = CIA WDBII, 1 = WVS
              ;    5th byte:    river = (flag >> 25) & 1: Values: 0 = not set, 1 = river-lake and GSHHS level = 2 (or WDBII class 0)
           ;    6th byte:    area magnitude scale p (as in 10^p) = flag >> 26.  We divide area by 10^p.
           west: 0L, $          ; West extent of polygon boundary in micro-degrees.
           east: 0L, $        ; East extent of polygon boundary in micro-degrees.
           south: 0L, $       ; South extent of polygon boundary in micro-degrees.
           north: 0L, $       ; North extent of polygon boundary in micro-degrees.
           area: 0L, $        ; Area of polygon in area/10^p = km^2.
           area_full: 0L, $   ; Area of original full-resolution polygon in area_full/10^p = km2.
           container: 0L, $   ; ID of container polygon that encloses this polygon (-1 if "none").
           ancestor: 0L }     ; ID of ancestor polygon in the full resolution set that was the source
                               ; of this polygon (-1 of "none").

; default values in structure
my_color = defined(color)     ? color     : 0
my_style = defined(linestyle) ? linestyle : 0
my_coord = defined(coord)     ? coord     : 'C' ; Celestial = eQuatorial
my_thick = defined(thick)     ? thick     : 1.0

t0 = systime(1)
; Open and read the GSHHS data file.
OPENR, lun, filename, /Get_Lun, /Swap_If_Little_Endian
npmin = ([1, 8, 16, 20, 20])[hires]
; count contours
count = 0L
tot = 0L
while (eof(lun) ne 1) do begin
    readu,lun, header
    polygon = LonArr(2, header.npoints, /NoZero)
    readu, lun, polygon
    tot += 1
    if (header.npoints ge npmin) then count += 1
endwhile
point_lun, lun, 0 ; rewind
; read data
i = 0L
stc_loc = 0
; in IDL use STRUCTURE, in GDL and FL use faster HASH
do_hash = hash_support()
buffer = do_hash ? handle_hash(/create) : 0
nb = 0L
nbmax = round(sqrt(count*1.0))
;print, tot, count,nbmax,npmin
t1 = systime(1)
WHILE (EOF(lun) NE 1) DO BEGIN
    READU, lun, header

   ; Get the polygon coordinates. Convert to lat/lon.
   polygon = LonArr(2, header.npoints, /NoZero)
   ;print, i, header.id, header.npoints
   READU, lun, polygon
   if (header.npoints ge npmin) then begin
       ; increment buffer
       sname = name+'_'+string(i,form='(i6.6)')
       st1 = {ra:   Reform(polygon[0,*] * 1.0e-6), $
              dec:  Reform(polygon[1,*] * 1.0e-6), $
              coord:my_coord, $
              linestyle:my_style, $
              color: my_color, $
              thick:my_thick}
       if do_hash then begin
           buffer[sname] = st1
       endif else begin
           buffer = merge_structures(buffer, st1, tag=sname)
       endelse
       nb += 1L
       if (nb ge nbmax && ~do_hash) then begin ; empty buffer
           stc_loc = merge_structures(stc_loc, buffer)
           buffer = 0
           nb = 0L
       endif
       i += 1L
   endif
ENDWHILE
t2 = systime(1)
if (nb gt 0) then begin; empty residual buffer
    if do_hash then begin
        stc_loc = handle_hash(buffer,/tostructure)
    endif else begin
        stc_loc = merge_structures(stc_loc, buffer)
    endelse
    buffer = 0
    nb = 0L
endif
Free_Lun, lun
t3 = systime(1)

; merge structures, or copy it
stc = merge_structures(stc, stc_loc)
t4 = systime(1)
;;print,t1-t0,t2-t1,t3-t2,t4-t3,t4-t0

return
end

; -----------------------------------------------------------------------------
FUNCTION hpx_map_getindex,indx, error
;
; Used to read in the index file for
; each map data file.  On successful completion, error is set to 0.
;
; Hacked from map_getindex in IDL 7.1
;
;COMPILE_OPT hidden

openr, lun, indx, /xdr, /get_lun, error = error
if error ne 0 then return, 0		;File not there or unreadable
segments=0L 
readu, lun, segments
dx_map=replicate({ fptr:0L, npts:0L,latmax:0.,latmin:0.,$
                   lonmax:0.,lonmin:0. }, segments )
readu, lun, dx_map
free_lun,lun
return, dx_map
END


; -----------------------------------------------------------------------------
pro hpx_continents_segments, fnames, name, hires, stc=stc, color=color, linestyle=linestyle, thick=thick
; create or expand a structure with Earth features
;
; fnames = [lowresname, hiresname]
; name = description for error message (boundaries, rivers, etc.)
; hires = 0 for low res, 1 for hires

;COMPILE_OPT hidden

;********************  find the files and read the index file **************
lun = -1
sub = (['low', 'high'])[hires]
;if strtrim(!dir,2) eq '' then !dir=getenv('IDL_DIR')
root_dir = strtrim(!dir,2) eq '' ? getenv('IDL_DIR') : !dir
fndx = FILEPATH(fnames[hires]+'.ndx', SUBDIR=['resource','maps',sub],root_dir=root_dir)
dat =  FILEPATH(fnames[hires]+'.dat', SUBDIR=['resource','maps',sub],root_dir=root_dir)
;print,fndx,dat
ndx = hpx_map_getindex(fndx, error)		;OPEN it
if error eq 0 then openr, lun, dat,/xdr,/stream, /get, error = error
if (error ne 0) and hires then begin 	;Try low res as a fallback
    message, 'High Res Map File: '+name+' not found, trying low res.', /INFO
    fndx = FILEPATH(fnames[0]+'.ndx', SUBDIR=['resource','maps','low'],root_dir=root_dir)
    dat =  FILEPATH(fnames[0]+'.dat', SUBDIR=['resource','maps','low'],root_dir=root_dir)
    ndx = hpx_map_getindex(fndx, error)		;OPEN it
endif                           ;Hires
if lun lt 0 then openr, lun, dat,/xdr,/stream, /get, error = error
if error ne 0 then message, 'Map file:'+fnames[hires]+' not found'
; print, fndx, dat




; ************** read the data file and fill structures ************
;

count = n_elements(ndx)
if count gt 0 then begin
    stc_loc = 0
    ; in IDL use STRUCTURE, in GDL and FL use faster HASH
    do_hash = hash_support()
    buffer = do_hash ? handle_hash(/create) : 0
    nb = 0L
    nbmax = round(sqrt(count*1.0))
    ;print,count,nbmax
; window,/free
; plot,/nodata,[-180,180],[-90,90]

    my_color = defined(color)     ? color     : 0
    my_style = defined(linestyle) ? linestyle : 0
    my_coord = defined(coord)     ? coord     : 'C' ; Celestial = eQuatorial
    my_thick = defined(thick)     ? thick     : 1.0
    ;print,'thick: ',my_thick

    for i=0L, count-1 do begin
        point_lun, lun, ndx[i].fptr
        if ndx[i].npts ge 2 then begin
            xy=fltarr(2,ndx[i].npts, /NOZERO)
            readu,lun,xy
            ; increment buffer
            sname  = name+'_'+string(i,form='(i6.6)')
            st1 = {ra:   reform(xy[1,*]), $
                   dec:  reform(xy[0,*]), $
                   coord:my_coord, $
                   linestyle:my_style, $
                   color: my_color, $
                   thick: my_thick}
            if do_hash then begin
                buffer[sname] = st1
            endif else begin
                buffer = merge_structures(buffer, st1, tag=sname)
            endelse
            nb += 1L
            if (nb ge nbmax && ~do_hash) then begin ; empty buffer
                stc_loc = merge_structures(stc_loc, buffer)
                buffer = 0
                nb = 0L
            endif
        endif
    endfor
    if (nb gt 0) then begin ; empty residual buffer
        if do_hash then begin
            stc_loc = handle_hash(buffer,/tostructure)
        endif else begin
            stc_loc = merge_structures(stc_loc, buffer)
        endelse
        buffer = 0
        nb = 0L
    endif
    ; merge structures (faster than incrementing in the loop)
    stc = merge_structures(stc, stc_loc)
endif

FREE_LUN, lun


return
end

; -----------------------------------------------------------------------------
function outline_earth, stc, $
  behavior=behavior, $
  COASTS=kcoasts, $
  color=color, $
  CONTINENTS = kcont, $
  COUNTRIES=kcountries, $
  HELP=help, $
  HIRES=khires, $
  linestyle=linestyle, $
  RIVERS=krivers, $
  thick=thick
; -----------------------------------------
; Hacked from map_continents in IDL 7.1
; -----------------------------------------

routine='outline_earth'
syntax = ["geo = outline_earth([geo_in],",$
          "              [BEHAVIOR=, COASTS=, COLOR=, CONTINENTS=, COUNTRIES=, ", $
          "               HELP=, HIRES=, LINESTYLE=, RIVERS=, THICK=])"]

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

; 			Map_continents keyword defaults:
cont   = keyword_set(kcont)
rivers = keyword_set(krivers)
coasts = keyword_set(kcoasts)
countries = keyword_set(kcountries)
;hires = keyword_set(khires)
hires = defined(khires) ? khires : 0

; ; continents by default
; if (countries+coasts+rivers eq 0) and n_elements(kcont) eq 0 then cont = 1

if ~(cont || rivers || coasts || countries) then begin
    print,syntax, form='(a)'
    print,'Type '
    print,routine+'(/help)'
    print,'for more information.'
    return,-1
endif

mode = defined(behavior) ? strupcase(behavior) : ''
if (mode ne 'IDL' && mode ne 'GDL') then begin
    if is_fl()  then mode = 'FL'
    if is_gdl() then mode = 'GDL'
    if is_idl() then mode = 'IDL'
endif


; if is_gdl() then begin
;     message,/info,'WARNING: '+routine+' is currently not supported in GDL'
;     return,stc
; endif
if (mode eq 'FL') then begin
    message,/info,'WARNING: Earth data currently not availale in FL '
    message,/info,'      set BEHAVIOR to IDL or GDL '
    message,/info,'      if those data are available '
    return,stc
endif

;ON_ERROR, 2



hires >= 0 ; set negative values to 0
if (mode eq 'IDL') then begin
    names  = ['low', 'high']
    hires <= (n_elements(names)-1)
    if (rivers) then begin	;Rivers
        hpx_continents_segments, 'r'+names, 'Rivers', hires, $
          stc=stc, color=color, linestyle=linestyle, thick=thick
    endif
    
    if (countries) then begin	;*** Countries ****
        hpx_continents_segments, 'b'+names, 'Boundaries', hires, $
          stc=stc, color=color, linestyle=linestyle, thick=thick
    endif
    
    if (coasts) then begin      ;Coasts
        hpx_continents_segments, 'c'+names, 'Coasts', hires, $
          stc=stc, color=color, linestyle=linestyle, thick=thick
    endif
    
    if (cont) then begin        ;Continents
        hpx_continents_segments, 'p'+names, 'Continents', hires, $
          stc=stc, color=color, linestyle=linestyle, thick=thick
    endif
endif
if (mode eq 'GDL') then begin
; resolution: among c(oarse), l(ow), i(ntermediate), h(igh), f(ull)
    names  = ['c','l','i','h','f']
    hires <= (n_elements(names)-1)
    if (coasts || cont) then begin ; Coasts or continents
        hpx_c_gdl, 'gshhs_'+names+'.b', 'Continents', hires, $
          stc=stc, color=color, linestyle=linestyle, thick=thick
    endif

    if (rivers) then begin ; Rivers
        hpx_c_gdl, 'wdb_rivers_'+names+'.b', 'Rivers', hires, $
          stc=stc, color=color, linestyle=linestyle, thick=thick
    endif

    if (countries) then begin ; Countries
        hpx_c_gdl, 'wdb_borders_'+names+'.b', 'Boundaries', hires, $
          stc=stc, color=color, linestyle=linestyle, thick=thick
    endif
endif


return, stc
end
