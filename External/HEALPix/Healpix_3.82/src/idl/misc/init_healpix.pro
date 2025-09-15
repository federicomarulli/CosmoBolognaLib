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
function defined_sysvar, sysvar
defsysv, sysvar, exist=is_defined
return, is_defined
end

function set_default_hpx_path, envar, dirlist

                                ; use env variable $HEALPIX and concatenate list of subdirectories
sep = path_sep()
path = getenv(envar)
if (strtrim(path,2) eq '') then begin
    print,'WARNING: environment variable '+envar+' not found. '
    path = ''
endif
if defined(dirlist) then begin
    for i=0, n_elements(dirlist)-1 do begin
        path = path + sep + dirlist[i]
    endfor
endif
                                ; end with separator
path = path + sep

return, path
end


pro init_healpix, verbose=verbose
;+
; defines the (structure) system variable Healpix
;
; 2006-05-22: print caution message in version 5.5a which has touchy expand_path
; 2006-oct : v 2.10, enabled nside > 8192
; 2008-oct: v2.11, use getenv
; 2009-09-09: v2.12: introduced !hpx* variables used by HFI-L2
; 2009-09-14: v2.12e, correct typo
; 2009-10-07: v2.12f, add !healpix.path.doc.(*) sub-structure
; 2015-03-17: v2.20,  add !healpix.path.src sub-structure
; 2019-02-28: v2.30,  adapted to Healpix 3.60 new C++ structure
; 2020-06-09: v2.40,  added path to doc/epub
;-

; system variable name
healpix_sysvar = '!HEALPIX'

; Healpix version
version = '3.82'

; release data
date = '2022-07-28'


; Healpix directory
directory = getenv('HEALPIX')
hexe      = getenv('HEXE')

hpx_path_data     = defined_sysvar('!hpx_path_data')    ? !hpx_path_data    : set_default_hpx_path('HEALPIX',['data'])
hpx_path_src      = defined_sysvar('!hpx_path_src')     ? !hpx_path_src     : set_default_hpx_path('HEALPIX',['src'])
hpx_path_test     = defined_sysvar('!hpx_path_test')    ? !hpx_path_test    : set_default_hpx_path('HEALPIX',['test'])
hpx_path_doc_html = defined_sysvar('!hpx_path_doc_html')? !hpx_path_doc_html: set_default_hpx_path('HEALPIX',['doc','html'])
hpx_path_doc_pdf  = defined_sysvar('!hpx_path_doc_pdf') ? !hpx_path_doc_pdf : set_default_hpx_path('HEALPIX',['doc','pdf'])
hpx_path_doc_epub = defined_sysvar('!hpx_path_doc_epub')? !hpx_path_doc_epub: set_default_hpx_path('HEALPIX',['doc','epub'])

hpx_path_bin_f90 = ''
if defined_sysvar('!hpx_path_bin_f90') then begin
    hpx_path_bin_f90 = !hpx_path_bin_f90
endif else begin
    if strtrim(hexe,2) ne '' then begin
        hpx_path_bin_f90 = set_default_hpx_path('HEXE')
    endif else begin
        if (strtrim(directory,2) ne '') then begin
            hpx_path_bin_f90 = set_default_hpx_path('HEALPIX',['bin'])
        endif ;else begin
;             print,'WARNING: environment variable HEALPIX not found '
;             directory = ''
;         endelse
    endelse
endelse

hpx_path_bin_cxx = ''
if defined_sysvar('!hpx_path_bin_cxx') then begin
    hpx_path_bin_cxx = !hpx_path_bin_cxx
endif else begin
    if (strtrim(directory,2) ne '') then begin
        if (float(version) le 3.501) then begin
            target    = getenv('HEALPIX_TARGET')
            if strtrim(target,2) eq '' then target = 'generic_gcc'
            hpx_path_bin_cxx = set_default_hpx_path('HEALPIX',['src','cxx',target,'bin']) ;up to 3.50
        endif else begin
            hpx_path_bin_cxx = set_default_hpx_path('HEALPIX',['bin']) ; 3.60 and above
        endelse
    endif
endelse

; give non blanck value to directory if HEALPIX not defined by !hpx_path_bin_f90
; is defined
if (strtrim(directory,2) eq '' && defined_sysvar('!hpx_path_bin_f90')) then begin
    directory = file_dirname(!hpx_path_bin_f90)+path_sep()
endif

stc_bin  = {cxx:hpx_path_bin_cxx,   f90:hpx_path_bin_f90}
stc_doc  = {html:hpx_path_doc_html, pdf:hpx_path_doc_pdf, epub:hpx_path_doc_epub}
stc_path = {bin:stc_bin, data: hpx_path_data, doc:stc_doc, src:hpx_path_src, test: hpx_path_test}


; list of possible Nside's
nside = 2L^lindgen(30)  ; 1, 2, 4, 8, ..., 8192, ..., 2^29 = 0.54e9

; flag for missing values
bad_value = -1.6375e30

comment = ['This system variable contains some information on Healpix :', $
           healpix_sysvar+'.VERSION   = current version number,', $
           healpix_sysvar+'.DATE      = date of release,',$
           healpix_sysvar+'.DIRECTORY = directory containing Healpix package,',$
           healpix_sysvar+'.PATH      = structure containing:',$
           healpix_sysvar+'.PATH.BIN  = structure containing binary path :',$
           healpix_sysvar+'.PATH.BIN.CXX  =     C++',$
           healpix_sysvar+'.PATH.BIN.F90  =     Fortran90',$
           healpix_sysvar+'.PATH.DATA = path to data subdirectory,',$
           healpix_sysvar+'.PATH.DOC  = path to doc subdirectories (.html, .pdf, .epub),',$
           healpix_sysvar+'.PATH.SRC  = path to src subdirectory,',$
           healpix_sysvar+'.PATH.TEST = path to test subdirectory,',$
           healpix_sysvar+'.NSIDE     = list of all valid values of Nside parameter,',$
           healpix_sysvar+'.BAD_VALUE = value of flag given to missing pixels in FITS files,',$
           healpix_sysvar+'.COMMENT   = this description.']

; create structure
stc = {version:version, date:date, directory:directory, path:stc_path, nside:nside, bad_value:bad_value, comment:comment}

; fill variable out
defsysv, healpix_sysvar, exists = exists
if (exists) then begin
    !Healpix = stc
endif else begin
    defsysv, healpix_sysvar, stc
endelse

if (keyword_set(verbose)) then begin

    print,'Initializing '+healpix_sysvar+' system variable'
    print
    print,comment,form='(a)'
    print
;     help,/st,healpix_sysvar
;     print
endif



return
end

