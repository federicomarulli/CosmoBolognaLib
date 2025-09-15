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
pro iprocess_mask, mask_in, distance_map $
              , binpath=binpath $
              , filled_mask=filled_mask $   
              , help=help $
              , hole_arcmin2=hole_arcmin2 $
              , hole_pixels=hole_pixels $
              , keep_tmp_files=keep_tmp_files $
              , nested=nested $
              , ordering=ordering $
              , ring=ring $
              , silent=silent $
              , tmpdir=tmpdir

;+
; NAME:
;    iprocess_mask
;
; PURPOSE:
;    interface to 'process_mask' F90 facility
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    iprocess_mask, mask_in, [distance_map , binpath=, filled_mask=, /help,
;    hole_arcmin2=, hole_pixels=, /keep_tmp_files, 
;       /nested, ordering=, /ring=,/silent=, tmpdir=]
;
; INPUT:
;   mask_in: input mask, can be a FITS file, or a memory array containing the
;        binary mask to process
;
; OUTPUT:
;   distance_map: (OPTIONAL) output distance map, can be a FITS file, or a
;   memory array. Will have the same ordering as the input mask
;
;       
; KEYWORD PARAMETERS:
;
;   binpath=: full path to back-end routine [default: $HEXE/process_mask, then $HEALPIX/bin/process_mask]
;                 a binpath starting with / (or \), ~ or $ is interpreted as absolute
;                 a binpath starting with ./ is interpreted as relative to current directory
;                 all other binpaths are relative to $HEALPIX
;
;   filled_mask=: OUTPUT mask with holes smaller than ... filled in.
;     Will have the same ordering as the input mask
;
;   /help:      if set, prints extended help
;
;   hole_arcmin2=: Minimal size (in arcmin^2) of invalid regions to be kept 
;    (can be used together with hole_pixels below, 
;     the result will be the largest of the two). Default=0.0
;
;   hole_pixels=: Minimal size (in pixels) of invalid regions to be kept 
;    (can be used together with hole_arcmin2 above, 
;     the result will be the largest of the two). Default=0
;
;   /keep_tmp_files: if set, temporary files are not discarded at the end of the
;                     run
;
;   /nested: if set, signals that *all* masks read online are in
;      NESTED scheme (does not apply to FITS file), see also /ring and Ordering
;
;   ordering=: either 'RING' or 'NESTED', ordering of online masks,
;    see /ring and /ordering
;
;   /ring: see /nested and ordering above
;
;   /silent:    if set, works silently
;
;   tmpdir=:      directory in which are written temporary files 
;         [default: IDL_TMPDIR (see IDL documentation about IDL_TMPDIR)]
;
;
; COMMON BLOCKS:
;    hxp_xface_com
;
; SIDE EFFECTS:
;    writes temporary parameter and data files, spawns external processes
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;     2012-02-24: 1st version, adapted from ismoothing
;     2012-03-01: retrofitted to run with GDL
;-

local = {routine: 'iprocess_mask', exe: 'process_mask', double: keyword_set(double)}
syntax = [local.routine+', mask_in  [, distance_map  ', $
'         , binpath=, filled_mask=, /help,', $
'            hole_arcmin2=, hole_pixels=, /keep_tmp_files, ', $
'            /nested, ordering=, /ring=,/silent=, tmpdir=]']

if keyword_set(help) then begin
    doc_library,local.routine
    return
endif

if (n_params() eq 0 or n_params() gt 2) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some input mask'
    return
endif
with_distance_map = (n_params() eq 2)
with_filled_mask  = defined(filled_mask) || arg_present(filled_mask)
if (~with_distance_map && ~with_filled_mask) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some output: Distance_map or FILLED_MASK='
    return
endif


;-------------------
hpx_xface_generic, fullpath, tmp_par_file, binpath, init=local, tmpdir=tmpdir

NoFile = " '' "
;;; if (~arg_present(distance_map)) then distance_map = NoFile

; deal with online data
tmp_mask_in       = hpx_mem2file(set_parameter(mask_in,    NoFile, /ifempty), /map,  /in, ring=ring, nested=nested, ordering=ordering)
tmp_distance_map  = hpx_mem2file(with_distance_map ? (defined(distance_map)?distance_map :-1) : NoFile, /out)
tmp_filled_mask   = hpx_mem2file(with_filled_mask  ? (defined(filled_mask )?filled_mask  :-1) : NoFile, /out)

; writes parameter file
openw,lunit,tmp_par_file, /get_lun
printf,lunit,'# parameter file for IDL interface to '+fullpath
printf,lunit,'# written: '+systime()+' by '+local.routine
printf,lunit,' '
;
printf,lunit,hpx_add_parameter('mask_file',     tmp_mask_in,      /expand)
;
printf,lunit,hpx_add_parameter('hole_min_size',       hole_pixels,           /skip_if_not_set)
printf,lunit,hpx_add_parameter('hole_min_surf_arcmin2', hole_arcmin2,        /skip_if_not_set) ; let code choose default
;
printf,lunit,hpx_add_parameter('distance_file', tmp_distance_map, def=NoFile,/ifempty,/expand,/force)
printf,lunit,hpx_add_parameter('filled_file',   tmp_filled_mask,  def=NoFile,/ifempty,/expand,/force)
free_lun, lunit

; execute command
hpx_xface_generic, /run, fullpath, tmp_par_file, silent=silent

; deal with online data
if (with_distance_map) then hpx_file2mem, tmp_distance_map, distance_map,/map
if (with_filled_mask)  then hpx_file2mem, tmp_filled_mask,  filled_mask, /map

; to_remove
hpx_xface_generic, clean = ~keyword_set(keep_tmp_files)

return
end
