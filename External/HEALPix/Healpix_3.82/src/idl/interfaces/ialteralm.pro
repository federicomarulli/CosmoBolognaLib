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
pro ialteralm, alm_in, alm_out $
              , beam_file_in   =  beam_file_in $
              , beam_file_out  =  beam_file_out $
              , binpath        = binpath $
              , coord_in       =  coord_in     $
              , coord_out      =  coord_out    $
              , epoch_in       =  epoch_in     $
              , epoch_out      =  epoch_out    $
              , fwhm_arcmin_in =  fwhm_arcmin_in  $
              , fwhm_arcmin_out=  fwhm_arcmin_out $
              , help           = help $
              , keep_tmp_files = keep_tmp_files $
              , lmax_out       =  lmax_out    $
              , nlmax_out      =  nlmax_out   $
              , nside_in       =  nside_in     $
              , nside_out      =  nside_out    $
              , nsmax_in       =  nsmax_in     $
              , nsmax_out      =  nsmax_out    $
              , silent         = silent $
              , tmpdir         = tmpdir $
              , windowfile_in  = windowfile_in   $
              , winfiledir_in  = winfiledir_in   $
              , windowfile_out = windowfile_out  $
              , winfiledir_out = winfiledir_out 

;+
; NAME:
;    ialteralm
;
; PURPOSE:
;    interface to 'alteralm' F90 facility
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    ialteralm, alm_in, alm_out, [beam_file_in=, beam_file_out=, 
;       binpath=, coord_in= , coord_out=, epoch_in=, epoch_out=,
;       fwhm_arcmin_in=, fwhm_arcmin_out=, /help,
;       /keep_tmp_files, lmax_out=, nlmax_out=, 
;       nside_in=, nside_out=, nsmax_in=, nsmax_out=
;       /silent=, tmpdir=, 
;       windowfile_in=,  winfiledir_in=, 
;       windowfile_out=, winfiledir_out=]
;
; INPUT:
;   alm_in: input alm, must be a FITS file
;
; OUTPUT:
;   alm_out: output alm, must be a FITS file
;
;       
; KEYWORD PARAMETERS:
;
;   binpath=   : full path to back-end routine [default: $HEXE/alteralm, then $HEALPIX/bin/alteralm]
;                 a binpath starting with / (or \), ~ or $ is interpreted as absolute
;                 a binpath starting with ./ is interpreted as relative to current directory
;                 all other binpaths are relative to $HEALPIX
;
;   beam_file_in= : beam window function of input alm
;                 either a FITS file or an array
;                 If defined and valid, it will override fwhm_arcmin_in
;                 Default: value of BEAM_LEG  keyword read from input alm FITS file
;
;   beam_file_out=: beam window function of output alm, 
;                 either a FITS file or an array. 
;                 If defined and valid, it will override fwhm_arcmin_out
;                 Default: '' (empty string, no beam window applied)
;
;   coord_in=     : Astrophysical coordinates system used to compute input alm.
;                 Case-insensitive single letter code.
;                 Valid choices are 'g','G'=Galactic, 'e','E'=Ecliptic,
;                 'c','q','C','Q'=Celestial/eQuatorial.
;                 Default: value of COORDSYS keyword read from input alm FITS file
;
;   coord_out=     : Astrophysical coordinates system of output alm.
;                 Default: coord_in
;
;   epoch_in=     : astronomical epoch of input coordinates (coord_in)
;                 Default: 2000.0
;
;   epoch_out=     : astronomical epoch of output coordinates (coord_out)
;                 Default: epoch_in
;
;   fwhm_arcmin_in= : Full Width Half-Maximum in arcmin of Gaussian beam applied
;                 to map from which are obtained input alm.
;                 Default: value of FWHM keyword in input alm FITS file
;
;   fwhm_arcmin_out=: FWHM in arcmin to be applied to output alm
;                 Default: fwhm_arcmin_in
;
;   /help:      if set, prints extended help
;
;   /keep_tmp_files: if set, temporary files are not discarded at the end of the
;                     run
;
;   lmax_out=, nlmax_out=  : maximum multipole of output alm
;
;   nside_in=, nsmax_in= : HEALPix resolution parameter of map
;                  from which were computed input alm
;                   Default: determined from input alm file
;
;   nside_out=, nsmax_out= : HEALPix resolution parameter Nside whose
;                  window function will be applied to output alm
;                  Could be set to 0 for infinitely small pixels (no window)
;                   Default: same as input nsmax_in
;
;   /silent     :if set, works silently
;
;   tmpdir=     :directory in which are written temporary files 
;         [default: IDL_TMPDIR (see IDL documentation about IDL_TMPDIR)]
;
;    windowfile_in=:    FITS file containing pixel window for nside_in
;        [default: determined automatically by back-end routine].
;      Do not set this keyword unless you really know what you are doing   
;
;    winfiledir_in=:   directory where windowfile_in is to be found 
;        [default: determined automatically by back-end routine].
;      Do not set this keyword unless you really know what you are doing
;
;    windowfile_out=:    FITS file containing pixel window for nside_out
;        [default: determined automatically by back-end routine].
;      Do not set this keyword unless you really know what you are doing   
;
;    winfiledir_out=:   directory where windowfile_out is to be found 
;        [default: determined automatically by back-end routine].
;      Do not set this keyword unless you really know what you are doing
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
;ialteralm,!healpix.path.test+'alm.fits','/tmp/alm_equat.fits',coord_in='g',coord_out='q'
;isynfast, 0, alm_in='/tmp/alm_equat.fits', '/tmp/map_equat.fits'
;mollview,'/tmp/map_equat.fits',1
;mollview,'/tmp/map_equat.fits',2
;
;   This example script reads the test (polarised) a_{lm} located in 
; $HEALPIX/test/alm.fits and rotates them from Galactic to Equatorial
;   coordinates, it then synthetizes a map out of those,
;   and finally plots its I and Q Stokes components (in Equatorial coordinates)
;
; MODIFICATION HISTORY:
;     2013-07-30: 1st version, adapted from iprocess_mask
;-

local = {routine: 'ialteralm', exe: 'alteralm', double: keyword_set(double)}
syntax = [local.routine+', alm_in, alm_out,  [beam_file_in=, beam_file_out=, ',$
'       binpath=, coord_in= , coord_out=, epoch_in=, epoch_out=,',$
'       fwhm_arcmin_in=, fwhm_arcmin_out=, /help,',$
'       /keep_tmp_files, lmax_out=, nlmax_out=,',$
'       nside_in=, nside_out=, nsmax_in=, nsmax_out=, ',$
'       /silent=, tmpdir=, ',$
'       windowfile_in=,  winfiledir_in=, ',$ 
'       windowfile_out=, winfiledir_out=]']

if keyword_set(help) then begin
    doc_library,local.routine
    return
endif

if (n_params() ne 2) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide input and output alm'
    return
endif

; check that no keyword is duplicated
solve_kw_conflict,'lmax_out',  'nlmax_out', k1=lmax_out, k2=nlmax_out, kout=nlmax_out, /defined
solve_kw_conflict,'nside_out', 'nsmax_out', k1=nside_out,k2=nsmax_out, kout=nsmax_out, /defined
solve_kw_conflict,'nside_in',  'nsmax_in',  k1=nside_in, k2=nsmax_in,  kout=nsmax_in,  /defined

;-------------------
hpx_xface_generic, fullpath, tmp_par_file, binpath, init=local, tmpdir=tmpdir

NoFile = " '' "

; deal with online data
tmp_beam_file_in  = hpx_mem2file(set_parameter(beam_file_in,  NoFile, /ifempty), /beam, /in)
tmp_beam_file_out = hpx_mem2file(set_parameter(beam_file_out, NoFile, /ifempty), /beam, /in)

; write parameter file
openw,lunit,tmp_par_file, /get_lun
printf,lunit,'# parameter file for IDL interface to '+fullpath
printf,lunit,'# written: '+systime()+' by '+local.routine
printf,lunit,' '
;
printf,lunit,hpx_add_parameter('infile_alms',     alm_in,      /expand)
printf,lunit,hpx_add_parameter('outfile_alms',    alm_out,     /expand, /force)
;
printf,lunit,hpx_add_parameter('beam_file_in',     tmp_beam_file_in,      /expand)
printf,lunit,hpx_add_parameter('beam_file_out',    tmp_beam_file_out,     /expand)
;
printf,lunit,hpx_add_parameter('coord_in',       coord_in,   /skip_if_not_set)
printf,lunit,hpx_add_parameter('coord_out',      coord_out,  /skip_if_not_set)
printf,lunit,hpx_add_parameter('epoch_in',       epoch_in,   /skip_if_not_set)
printf,lunit,hpx_add_parameter('epoch_out',      epoch_out,  /skip_if_not_set)
printf,lunit,hpx_add_parameter('fwhm_arcmin_in', fwhm_arcmin_in,   /skip_if_not_set)
printf,lunit,hpx_add_parameter('fwhm_arcmin_out',fwhm_arcmin_out,  /skip_if_not_set)
printf,lunit,hpx_add_parameter('nsmax_in',       nsmax_in,   /skip_if_not_set)
printf,lunit,hpx_add_parameter('nsmax_out',      nsmax_out,  /skip_if_not_set)
printf,lunit,hpx_add_parameter('nlmax_out',      nlmax_out,  /skip_if_not_set)
;
printf,lunit,hpx_add_parameter('winfiledir_in', winfiledir_in, /expand, /skip_if_not_set)
printf,lunit,hpx_add_parameter('windowfile_in', windowfile_in,          /skip_if_not_set)
printf,lunit,hpx_add_parameter('winfiledir_out', winfiledir_out, /expand, /skip_if_not_set)
printf,lunit,hpx_add_parameter('windowfile_out', windowfile_out,          /skip_if_not_set)
;
free_lun, lunit

; execute command
hpx_xface_generic, /run, fullpath, tmp_par_file, silent=silent

; to_remove
hpx_xface_generic, clean = ~keyword_set(keep_tmp_files)

return
end
