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
pro isynfast, cl_in, map_out $
              , alm_in=alm_in $   ;,   almsfile=almsfile $
              , alm_out=alm_out $ ;, outfile_alms=outfile_alms $
              , apply_windows=apply_windows $
              , beam_file=beam_file $
              , binpath=binpath $
              , double=double $   
              , fwhm_arcmin=fwhm_arcmin $
              , help=help $
              , iseed=iseed $
              , keep_tmp_files=keep_tmp_files $
              , lmax=lmax,       nlmax=nlmax $
              , nside=nside,     nsmax=nsmax $
              , plmfile=plmfile $
              , simul_type=simul_type $
              , silent=silent $
              , tmpdir=tmpdir $
              , windowfile=windowfile $
              , winfiledir=winfiledir 

;+
; NAME:
;    isynfast
;
; PURPOSE:
;    interface to 'synfast' F90 facility
;      can be used to    C(l) -> map [ & alm]
;                        C(l) -> alm
;                        alm  -> map
;                   C(l)+alm  -> map
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    isynfast, cl_in [, map_out,
;      alm_in=,  alm_out=, apply_windows=, beam_file=, binpath=, double=, fwhm_arcmin=, help=,
;      iseed=, keep_tmp_files=, lmax=, nlmax=, nside=, nsmax=, plmfile=,
;      simul_type=, silent=, tmpdir=, windowfile=, winfiledir=]
;
;
; INPUTS:
;   cl_in: input power spectrum, can be a FITS file, or a memory array containing the
;        C(l), used to generate a map or a set of gaussian alm.
;   If empty quotes ('') or a zero (0) are provided, it will be interpreted as "No input C(l)", in
;   which case some input alm's (alm_in) are required.
;
; OPTIONAL OUTPUTS:
;    map_out : map synthetised from the power spectrum or from constraining alm
;              either a FITS file or an array
;
; KEYWORD PARAMETERS:
;
;   alm_in=:    optional input (constraining) alm           [default: no alm]
;               must be a FITS file
;
;   alm_out=:   contains on output the effective alm
;               must be a FITS file
;
;   /apply_windows: if set, beam and pixel windows are applied to input alm_in
;   (if any)
;
;   beam_file=: beam window function, either a FITS file or an array
;
;   binpath=: full path to back-end routine [default: $HEXE/synfast, then $HEALPIX/bin/synfast]
;                 a binpath starting with / (or \), ~ or $ is interpreted as absolute
;                 a binpath starting with ./ is interpreted as relative to current directory
;                 all other binpaths are relative to $HEALPIX
;
;   /double:    if set, I/O is done in double precision [default: single
;       precision I/O]
;
;   fwhm_arcmin=: gaussian beam FWHM in arcmin [default = 0]
;
;   /help:      if set, prints extended help
;
;   iseed=: integer seed of radom sequence [default = 0]
;
;   /keep_tmp_files: if set, temporary files are not discarded at the end of the
;                     run
;
;   lmax=,nlmax= : maximum multipole of simulation [default: 2*Nside]
;
;   nside=,nsmax= : Healpix resolution parameter
;
;   plmfile=: FITS file containing precomputed Spherical Harmonics [default = no file]
;
;   simul_type=: 
;        1) Temperature only
;        2) Temperature + polarisation [default]
;        3) Temperature + 1st derivatives
;        4) Temperature + 1st & 2nd derivatives
;        5) T+P + 1st derivatives
;        6) T+P + 1st & 2nd derivates
;
;    /silent:    if set, works silently
;
;    tmpdir=:      directory in which are written temporary files 
;         [default: IDL_TMPDIR (see IDL documentation about IDL_TMPDIR)]
;
;    windowfile=:    FITS file containing pixel window 
;        [default: determined automatically by back-end routine].
;      Do not set this keyword unless you really know what you are doing   
;
;    winfiledir=:     directory where the pixel windows are to be found 
;        [default: determined automatically by back-end routine].
;      Do not set this keyword unless you really know what you are doing
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
;    isynfast, '$HEALPIX/test/cl.fits', map, fwhm=30, nside=256, /silent
;    mollview, map, 1, title='I'
;    mollview, map, 2, title='Q'
;
; will synthetize and plot I and Q  maps constistent with WMAP-1yr best fit
; power spectrum and observed with a circular gaussian 30 arcmin beam.
;
; MODIFICATION HISTORY:
;
;	July 2007, Eric Hivon, v1.0
;	Aug  2007, Eric Hivon, v1.1
;       Oct  2007, v1.1.1
;       Feb  2008, v1.1.2 : accept Cl array
;       May  2008, v1.1.3 (post beta release): bug correction for map_out
;       2008-08-27, v1.1.3 (post 2.10 release): do not try to read map_out in
;                 alm-only generation mode
;       2012-03-01: retrofitted to run with GDL
;-
local = {routine: 'isynfast', exe: 'synfast', double: keyword_set(double)}
syntax = [local.routine+', cl_in [, map_out,  ', $
'  alm_in=,  alm_out=, apply_windows=, beam_file=, binpath=, double=, fwhm_arcmin=, help=, $', $
'  iseed=, keep_tmp_files=, lmax=, nlmax=, nside=, nsmax=, plmfile=, $', $
'  simul_type=, silent=, tmpdir=, windowfile=, winfiledir= ]']

if keyword_set(help) then begin
    doc_library,local.routine
    return
endif

if (n_params() eq 0 || n_params() gt 2) then begin
    print,syntax,form="(a)"
    return
endif

if (undefined(cl_in) && undefined(alm_in)) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some input : cl_in or alm_in'
    return
endif
with_map_out = (n_params() eq 2)
if (~with_map_out  && undefined(alm_out)) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some output: map_out or alm_out'
    return
endif

; check that no keyword is duplicated
;solve_kw_conflict,'alm_in', 'almsfile',     k1=alm_in,  k2=almsfile, kout=alm_in, /defined
solve_kw_conflict,'lmax',   'nlmax',         k1=lmax,    k2=nlmax,    kout=lmax, /defined
solve_kw_conflict,'nside',  'nsmax',         k1=nside,   k2=nsmax,    kout=nside, /defined

;-------------------
hpx_xface_generic, fullpath, tmp_par_file, binpath, init=local, tmpdir=tmpdir

NoFile = " '' "
;if (~arg_present(map_out)) then map_out = NoFile

; deal with online data
x_cl_in = set_parameter(cl_in,     NoFile, /ifempty, /ifzero)
tmp_cl_in     = hpx_mem2file(x_cl_in, /cl,   /in)
tmp_beam_file = hpx_mem2file(set_parameter(beam_file, NoFile, /ifempty), /beam, /in)
tmp_map_out   = hpx_mem2file(with_map_out ? (defined(map_out)?map_out :-1) : NoFile,  /out)

no_cl = undefined(x_cl_in) || (size(/tname,x_cl_in) eq 'STRING' && x_cl_in eq NoFile)
if (no_cl && undefined(alm_in)) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some input : cl_in or alm_in'
    return
endif

; writes parameter file
openw,lunit,tmp_par_file, /get_lun
printf,lunit,'# parameter file for IDL interface to '+fullpath
printf,lunit,'# written: '+systime()+' by '+local.routine
printf,lunit,' '
tmp_nside                     = set_parameter(nside,       128)
printf,lunit,hpx_add_parameter('nsmax',       tmp_nside)
printf,lunit,hpx_add_parameter('simul_type', simul_type,  def=2,           /ifempty)
printf,lunit,hpx_add_parameter('nlmax',       lmax,        def=2*tmp_nside, /ifempty)
printf,lunit,hpx_add_parameter('iseed',       iseed,       def=0,           /ifempty)
printf,lunit,hpx_add_parameter('fwhm_arcmin',fwhm_arcmin, def=0,           /ifempty)
printf,lunit,hpx_add_parameter('apply_windows',apply_windows, def='F',    /ifempty)
;
printf,lunit,hpx_add_parameter('winfiledir', winfiledir, /expand, /skip_if_not_set)
printf,lunit,hpx_add_parameter('windowfile', windowfile,          /skip_if_not_set)
;
printf,lunit,hpx_add_parameter('infile',     tmp_cl_in,      /expand)
printf,lunit,hpx_add_parameter('beam_file',  tmp_beam_file, /expand)
             
printf,lunit,hpx_add_parameter('almsfile',  alm_in,   def=NoFile, /ifempty, /expand)
printf,lunit,hpx_add_parameter('plmfile',   plmfile,  def=NoFile, /ifempty, /expand)
;            
printf,lunit,hpx_add_parameter('outfile',      tmp_map_out,             /expand, /force)
printf,lunit,hpx_add_parameter('outfile_alms', alm_out,    def=NoFile, /expand, /force, /ifempty)
free_lun, lunit

; execute command
hpx_xface_generic, /run, fullpath, tmp_par_file, silent=silent

; deal with online data
if (with_map_out) then hpx_file2mem, tmp_map_out, map_out,/map ; 2008-08-27

; to_remove
hpx_xface_generic, clean = ~keyword_set(keep_tmp_files)

return
end
