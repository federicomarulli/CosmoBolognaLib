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
pro ismoothing, map_in, map_out $
              , beam_file=beam_file $
              , binpath=binpath $
              , double=double $   
              , fwhm_arcmin=fwhm_arcmin $
              , help=help $
              , iter_order = iter_order $
              , keep_tmp_files=keep_tmp_files $
              , lmax=lmax,      nlmax=nlmax $
              , nested=nested $
              , ordering=ordering $
              , plmfile=plmfile $
              , regression = regression $
              , ring=ring $
              , simul_type=simul_type $
              , silent=silent $
              , theta_cut_deg = theta_cut_deg $
              , tmpdir=tmpdir $
              , won = won $
              , w8file = w8file $
              , w8dir = w8dir  ; w8filedir

;+
; NAME:
;    ismoothing
;
; PURPOSE:
;    interface to 'smoothing' F90 facility
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    ismoothing, map_in, map_out [, beam_file=, binpath=, double=, fwhm_arcmin=, help=, iter_order=, keep_tmp_files=, 
;       lmax=, nlmax=, nested=, ordering=, plmfile=, regression=, ring=, 
;       simul_type=, silent=, theta_cut_deg=, tmpdir=, won=, w8file=, w8dir=]
;
; INPUT:
;   map_in: input map, can be a FITS file, or a memory array containing the
;        map to smooth
;
; OUTPUT:
;   map_out: output smoothed map, can be a FITS file, or a memory array
;   Will have the same orderin as the input map
;       
; KEYWORD PARAMETERS:
;
;   beam_file=: beam window function, either a FITS file or an array
;
;   binpath=: full path to back-end routine [default: $HEXE/smoothing, then $HEALPIX/bin/smoothing]
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
;   iter_order=: order of iteration in the analysis [default = 0]
;
;   /keep_tmp_files: if set, temporary files are not discarded at the end of the
;                     run
;
;   lmax=,nlmax= : maximum multipole of smoothing [default: determined by back-end routine (ie, smoothing)]
;
;   /nested: if set, signals that *all* maps and mask read online are in
;      NESTED scheme (does not apply to FITS file), see also /ring and Ordering
;
;   ordering=: either 'RING' or 'NESTED', ordering of online maps and masks,
;    see /ring and /ordering
;
;   plmfile=: FITS file containing precomputed Spherical Harmonics [default = no file]
;
;   regression= 0, 1 or 2, regress out best fit monopole and/or dipole before
;       alm analysis
;     [default = 0, analyze raw map] 
;
;   /ring: see /nested and ordering above
;
;   simul_type=: 1 or 2, smooths temperature only or temperature + polarization
;
;   /silent:    if set, works silently
;
;   theta_cut_deg: cut around the equatorial plane 
;
;   tmpdir=:      directory in which are written temporary files 
;         [default: IDL_TMPDIR (see IDL documentation about IDL_TMPDIR)]
;
;   won=:     if 0, no weighting, 
;              if set to 1, a ring-based  quadrature weighting scheme is applied,
;              if set to 2, a pixel-based quadrature weighting scheme is applied.
;       [default: 1 apply ring-based weighting]
;
;   w8file=:    FITS file containing weights 
;        [default: determined automatically by back-end routine] 
;      do not set this keyword unless you really know what you are doing             
;
;   w8dir=:     directory where the weights are to be found 
;        [default: determined automatically by back-end routine]
;      do not set this keyword unless you really know what you are doing             
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
;     whitenoise = randomn(seed, nside2npix(256))
;     ismoothing, whitenoise, rednoise, fwhm=120, /ring, simul=1,/silent
;     mollview, whitenoise,title='White noise'
;     mollview, rednoise,  title='Smoothed white Noise'
;
; will generate and plot a white noise map and its smoothed version
;
; MODIFICATION HISTORY:
;
;	July 2007, Eric Hivon, v1.0
;	Aug  2007, Eric Hivon, v1.1
;       Oct  2007, v1.1.1
;       2009-09-07:  w8filedir -> w8dir *EVERYWHERE*
;       2018-05-23: deal with pixel-based quadrature weights
;-

local = {routine: 'ismoothing', exe: 'smoothing', double: keyword_set(double)}
syntax = [local.routine+', map_in, map_out  ', $
'[ ,beam_file=, binpath=, double=, fwhm_arcmin=, help=, iter_order=, keep_tmp_files=, ', $
'  lmax=, nlmax=, nested=, ordering=, plmfile=, regression=, ring=, ', $
'  simul_type=, silent=, theta_cut_deg=, tmpdir=, won=, w8file=, w8dir=]']

if keyword_set(help) then begin
    doc_library,local.routine
    return
endif

if (n_params() ne 2) then begin
    print,syntax,form='(a)'
    print,local.routine+': Should provide some input and output maps'
    return
endif

; check that no keyword is duplicated
solve_kw_conflict,'lmax',   'nlmax',         k1=lmax,    k2=nlmax,    kout=lmax, /defined

;-------------------
hpx_xface_generic, fullpath, tmp_par_file, binpath, init=local, tmpdir=tmpdir

NoFile = " '' "
;;; if (~arg_present(map_out)) then map_out = NoFile

; deal with online data
tmp_map_in    = hpx_mem2file(set_parameter(map_in,    NoFile, /ifempty), /map,  /in, ring=ring, nested=nested, ordering=ordering)
tmp_beam_file = hpx_mem2file(set_parameter(beam_file, NoFile, /ifempty), /beam, /in)
tmp_map_out   = hpx_mem2file(map_out,                                           /out)

; writes parameter file
openw,lunit,tmp_par_file, /get_lun
printf,lunit,'# parameter file for IDL interface to '+fullpath
printf,lunit,'# written: '+systime()+' by '+local.routine
printf,lunit,' '
printf,lunit,hpx_add_parameter('simul_type', simul_type,      def=2,   /ifempty)
printf,lunit,hpx_add_parameter('fwhm_arcmin',fwhm_arcmin,     def=0,   /ifempty)
printf,lunit,hpx_add_parameter('iter_order',   iter_order,    def=0,   /ifempty)
printf,lunit,hpx_add_parameter('theta_cut_deg',theta_cut_deg, def=0.0, /ifempty)
printf,lunit,hpx_add_parameter('won',           won,           def=1,   /ifempty)
printf,lunit,hpx_add_parameter('plmfile',   plmfile,           def=NoFile, /ifempty, /expand)
;
printf,lunit,hpx_add_parameter('infile',     tmp_map_in,      /expand)
printf,lunit,hpx_add_parameter('beam_file',  tmp_beam_file, /expand)
;
printf,lunit,hpx_add_parameter('nlmax',       lmax,             /skip_if_not_set)
printf,lunit,hpx_add_parameter('w8filedir', w8dir,     /expand, /skip_if_not_set)
printf,lunit,hpx_add_parameter('w8file',    w8file,             /skip_if_not_set)
printf,lunit,hpx_add_parameter('regression', regression,        /skip_if_not_set) ; let code choose default
;
printf,lunit,hpx_add_parameter('outfile',      tmp_map_out,             /expand, /force)
free_lun, lunit

; execute command
hpx_xface_generic, /run, fullpath, tmp_par_file, silent=silent

; deal with online data
hpx_file2mem, tmp_map_out, map_out,/map

; to_remove
hpx_xface_generic, clean = ~keyword_set(keep_tmp_files)

return
end
