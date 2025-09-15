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
pro ianafast, map1_in, cl_out $
              , alm1_out = alm1_out $
              , alm2_out = alm2_out $
              , binpath=binpath $
              , cxx=cxx $
              , double = double $
              , healpix_data=healpix_data $
              , help=help $
              , iter_order = iter_order $
              , keep_tmp_files=keep_tmp_files $
              , map2_in = map2_in $
              , maskfile = maskfile $
              , nested=nested $
              , nlmax = nlmax $
              , nmmax = nmmax $
              , ordering=ordering $
              , plmfile = plmfile $
              , polarisation = polarisation $
              , regression = regression $
              , ring=ring $
              , show_cl=show_cl $
              , simul_type = simul_type $
              , silent=silent $
              , theta_cut_deg = theta_cut_deg $
              , tmpdir=tmpdir $
              , weighted=weighted $
              , won = won $
              , w8file = w8file $
              , w8dir = w8dir  ; w8filedir

;+
; NAME:
;    ianafast
;
; PURPOSE:
;    interface to 'anafast' F90 and 'anafast_cxx' C++ facilities
;
; CATEGORY:
;    IDL Interface to external facility
;
; CALLING SEQUENCE:
;    ianafast, map1_in [, cl_out,
;       alm1_out=, alm2_out=, binpath=, cxx=, double=, help=, healpix_data=, iter_order=, keep_tmp_files=, 
;       map2_in=, maskfile=, nested=, nlmax=, nmmax=, ordering=, plmfile=, polarisation=, 
;       regression=, ring=, show_cl=, simul_type=, silent=, theta_cut_deg=, tmpdir=, 
;       weighted=, won=, w8file=, w8dir= ]
;
;
; INPUTS:
;   map1_in: 1st input map, can be a FITS file, or a memory array containing the
;        map to analyze
;
; OPTIONAL OUTPUTS:
;    Cl_out : output auto or cross power spectrum
;
; KEYWORD PARAMETERS:
;
;   alm1_out=:   output alm of 1st map          [default: alm not kept]
;                must be a FITS file
;
;   alm2_out=:   output alm of 2nd map (if any) [default: alm not kept]
;                must be a FITS file
;
;   binpath=: full path to back-end routine [default: $HEXE/anafast, then $HEALPIX/bin/anafast 
;                   or $HEALPIX/bin/anafast_cxx if cxx is set]
;                 a binpath starting with / (or \), ~ or $ is interpreted as absolute
;                 a binpath starting with ./ is interpreted as relative to current directory
;                 all other binpaths are relative to $HEALPIX
;
;   /cxx: if set, the C++ back-end anafast_cxx is invoked instead of F90 anafast, 
;              AND the parameter file is written accordingly
;
;   /double:    if set, I/O is done in double precision [default: single precision I/O]
;
;   /help:      if set, prints extended help
;
;   healpix_data=: same as w8dir
;
;   iter_order=: order of iteration in the analysis [default = 0]
;
;   /keep_tmp_files: if set, temporary files are not discarded at the end of the
;                     run
;
;    map2_in=: 2nd input map (FITS file or array), if provided, Cl_out will
;     contain the cross power spectra of the 2 maps [default: no 2nd map]
;
;    maskfile=: pixel mask (FITS file or array), F90 only   [default: no mask]
;
;    /nested=: if set, signals that *all* maps and mask read online are in
;      NESTED scheme (does not apply to FITS file), see also /ring and Ordering
;
;    nlmax=:   maximum multipole of analysis, *required* for C++ anafast_cxx,
;         optional for F90 anafast
;
;    nmmax=:   maximum degree m, only valid for C++ anafast_cxx [default: nlmax]
;
;    ordering=: either 'RING' or 'NESTED', ordering of online maps and masks,
;    see /ring and /ordering
;
;    plmfile=: FITS file containing precomputed Spherical Harmonics [default = no file]
;
;    /polarisation: if set analyze temperature + polarization (same as simul_type = 2)
;
;    regression= 0, 1 or 2, regress out best fit monopole and/or dipole before
;       alm analysis
;     [default = 0, analyze raw map] 
;
;    /ring: see /nested and ordering above
;
;    /show_cl: if set, and cl_out is defined, l*(l+1) C(l)/2Pi is plotted
;
;    simul_type=: 1 or 2, analyze temperature only or temperature + polarization
;
;    /silent:    if set, works silently
;
;    theta_cut_deg=: cut around the equatorial plane, only in F90
;
;    tmpdir=:      directory in which are written temporary files 
;         [default: IDL_TMPDIR (see IDL documentation about IDL_TMPDIR)]
;
;    weighted=:     same as won
;       [default: apply weighting]
;
;    won=:     if 0, no weighting, 
;              if set to 1, a ring-based  quadrature weighting scheme is applied,
;              if set to 2, a pixel-based quadrature weighting scheme is applied.
;       [default: 1 apply ring-based weighting]
;
;    w8file=:    FITS file containing weights (applicable when won or weighted non-zero)
;       In F90: [default: determined automatically by back-end routine] 
;          do not set this keyword unless you really know what you are doing  
;       In C++: must be set to full path of weight file, consistent with value of won (or weighted)       
;
;    w8dir=:     directory where the weights are to be found, in F90 only 
;        [default: determined automatically by back-end routine]
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
;    whitenoise = randomn(seed, nside2npix(256))
;    ianafast, whitenoise, cl, /ring, /silent
;    plot, cl[*,0]
;
; will plot the power spectrum of a white noise map
;
; MODIFICATION HISTORY:
;
;	July 2007, Eric Hivon, v1.0
;	Aug  2007, Eric Hivon, v1.1
;       Oct  2007, v1.1.1
;       Oct  2007, v1.1.2: addition of /show_cl
;       July 2008, v1.2:   output Cl file really optional
;       2008-08-27: cosmetics (a closing quote was missing)
;       2009-04-30: tmpdir= is not ignored anymore, and IDL_TMPDIR is now used
;       as default temporary directory.
;       2009-09-07:  w8filedir -> w8dir *EVERYWHERE*
;       2009-09-09:  use !healpix.path.data instead of !healpix.directory+'/data'
;       2010-02-22: SILENT forwarded to hpx_file2mem
;       2012-03-01: retrofitted to run with GDL
;       2018-05-23: deal with pixel-based quadrature weights
;
;-
local = {routine: 'ianafast', exe: 'anafast', exe_cxx: 'anafast_cxx', double: keyword_set(double)}
syntax = [local.routine+', map1_in [, cl_out, ',$
'  alm1_out=, alm2_out=, binpath=, cxx=, double=, help=, healpix_data=, iter_order=, keep_tmp_files=, ', $
'  map2_in=, maskfile=, nested=, nlmax=, nmmax=, ordering=, plmfile=, polarisation=, ', $
'  regression=, ring=, show_cl=, simul_type=, silent=, theta_cut_deg=, tmpdir=, ', $
'  weighted=, won=, w8file=, w8dir= ]' ]

if keyword_set(help) then begin
    doc_library,local.routine
    return
endif

if (n_params() eq 0) then begin
    print,syntax,form='(a)'
    print,local.routine+',/help     for extended help'
    return
endif
if (n_params() gt 2) then begin
    print,syntax,form='(a)'
    message,'wrong number of arguments'
endif
if (keyword_set(cxx)) then begin
    if undefined(nlmax) then message,'nlmax is REQUIRED when invoking C++ code (/cxx option).'
endif
with_cl_out = (n_params() ge 2)

; check that no keyword is duplicated
solve_kw_conflict,'w8dir', 'healpix_data',       k1=w8dir,    k2=healpix_data,   kout=w8dir, /defined
solve_kw_conflict,'won',     'weighted',          k1=won,    k2=weighted,         kout=won, /defined
solve_kw_conflict,'simul_type', 'polarisation', k1=simul_type, k2=polarisation, kout=sim_final,/defined
if (undefined(sim_final)) then sim_final = 1 else sim_final = defined(polarisation) ?  polarisation + 1 : simul_type
pol_bool = (['false','true'])[(sim_final-1)<1>0]

;-------------------
hpx_xface_generic, fullpath, tmp_par_file, binpath, init=local, cxx=cxx, tmpdir=tmpdir


NoFile = keyword_set(cxx) ? " " : " '' "

if (keyword_set(cxx) && keyword_set(won) && undefined(w8file)) then begin
    message,'w8file must be defined when setting won or weighted with C++ code (/cxx option).'
endif

if (keyword_set(cxx) && defined(theta_cut_deg)) then begin
    message,/info,'WARNING: theta_cut_deg will be ignored by anafast_cxx'
endif
if (keyword_set(cxx) && defined(maskfile)) then begin
    message,/info,'WARNING: maskfile will be ignored by anafast_cxx'
endif

; deal with online data
tmp_cl_out   = hpx_mem2file(with_cl_out ? (defined(cl_out)?cl_out :-1) : NoFile, /out)
tmp_map1_in  = hpx_mem2file(set_parameter(map1_in, NoFile, /ifempty),  /in, /map, ring=ring, nested=nested, ordering=ordering)
tmp_map2_in  = hpx_mem2file(set_parameter(map2_in, NoFile, /ifempty),  /in, /map, ring=ring, nested=nested, ordering=ordering)
tmp_maskfile = hpx_mem2file(set_parameter(maskfile, NoFile,/ifempty), /in, /map, ring=ring, nested=nested, ordering=ordering)

; writes parameter file
openw,lunit,tmp_par_file, /get_lun
printf,lunit,'# parameter file for IDL interface to '+fullpath
printf,lunit,'# written: '+systime()+' by '+local.routine
printf,lunit,' '
; ------- F90 and C++ ------
printf,lunit,hpx_add_parameter('infile',    tmp_map1_in ,/expand)
printf,lunit,hpx_add_parameter('nlmax',       nlmax,      /skip_if_not_set) ; let code choose default
printf,lunit,hpx_add_parameter('iter_order',   iter_order,    def=0, /ifempty)
;printf,lunit,'outfile       = '+add_prefix(tmp_cl_out,'!',except=NoFile)
printf,lunit,hpx_add_parameter('outfile',      tmp_cl_out,             /expand, /force)
printf,lunit,hpx_add_parameter('outfile_alms', alm1_out,    def=NoFile, /expand, /force, /ifempty)
if (~keyword_set(cxx)) then begin
; ------- F90 only --------
    printf,lunit,hpx_add_parameter('infile2',   tmp_map2_in ,/expand)
    printf,lunit,' '
    printf,lunit,hpx_add_parameter('maskfile',  tmp_maskfile,                       /expand)
    printf,lunit,hpx_add_parameter('plmfile',     plmfile,    def=NoFile, /ifempty, /expand)
    printf,lunit,hpx_add_parameter('simul_type', sim_final ) ; F90 only

    printf,lunit,hpx_add_parameter('regression', regression, /skip_if_not_set) ; let code choose default
;
    printf,lunit,hpx_add_parameter('w8filedir', w8dir, /expand,     /skip_if_not_set) ; F90 only
    printf,lunit,hpx_add_parameter('w8file',    w8file,             /skip_if_not_set) ; F90 only
;
    printf,lunit,hpx_add_parameter('theta_cut_deg',theta_cut_deg, def=0.0, /ifempty) ; F90 only
    printf,lunit,hpx_add_parameter('won',           won,           def=1, /ifempty) ; F90 only
    printf,lunit,hpx_add_parameter('outfile_alms2', alm2_out,    def=NoFile, /expand, /force, /ifempty)
    printf,lunit,' '
endif else begin
; ------- C++ only -------
    printf,lunit,hpx_add_parameter('polarisation', pol_bool)
    if (pol_bool eq 'true') then begin
        printf,lunit,hpx_add_parameter('full_powerspectrum','true' )
    endif
    if (defined(regression) && regression ge 1) then begin
        printf,lunit,hpx_add_parameter('remove_monopole','true' )
    endif
    printf,lunit,hpx_add_parameter('double_precision',keyword_set(double)?'true':'false')
    printf,lunit,hpx_add_parameter('nmmax',       nmmax,      /skip_if_not_set) ; let code choose default
;    printf,lunit,hpx_add_parameter('weighted',keyword_set(won)?'true':'false')
    if (keyword_set(won)) then begin
;        printf,lunit,hpx_add_parameter('healpix_data', !healpix.directory+'/data', /expand, /skip_if_not_set)
        printf,lunit,hpx_add_parameter('healpix_data', !healpix.path.data, /expand, /skip_if_not_set)
        if (won eq 1) then printf,lunit,hpx_add_parameter('ringweights',  w8file, /skip_if_not_set)
        if (won eq 2) then printf,lunit,hpx_add_parameter('pixelweights', w8file, /skip_if_not_set)
    endif
endelse
; ------------------------
free_lun, lunit

; execute command
hpx_xface_generic, /run, fullpath, tmp_par_file, silent=silent

; deal with online data
if (with_cl_out) then hpx_file2mem, tmp_cl_out, cl_out, /cl, show_cl = show_cl, silent=silent

; to_remove
hpx_xface_generic, clean = ~keyword_set(keep_tmp_files)



return
end


