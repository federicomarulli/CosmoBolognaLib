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
pro alm_t2i, alm_tab, index, alm_list, HELP=help, MFIRST=mfirst
;+
; NAME:
;     alm_t2i
;
; PURPOSE:
;  turn an alm COMPLEX or REAL array of the form (respectively)
;        a[l,m,s] or a[l,m,0:1,s] (default)
;     or a[m,l,s] or a[m,l,0:1,s] (if MFIRST is set)
;  with s being the Stokes index (generally either 0 or {0,1,2})
;
;  into an index list (with index = l*(l+1)+m+1 and abs(m) <= l)
;  and a 2D REAL array alm_list
;    where alm_list[*,0,s] is the real part of alm
;    and   alm_list[*,1,s] is the imaginary part of alm
;
;
; index and alm_list can then be used by alm2fits
;
; internally Healpix use complex arrays of the form a(l,m)
; (where l is the first index ie, the fastest varying one)
;
; CALLING SEQUENCE:
;     alm_t2i, alm_tab, index, alm_list
; 
; INPUTS:
;     alm_tab: real or complex array, containing all the alm for l
;         in [0,lmax] and m in [0,mmax] (and s in [0,smax] if applicable)
;       if REAL    it has 3 (or 4) dimensions
;       if COMPLEX is has 2 (or 3) dimensions
;
; KEYWORD PARAMETERS
;
;    HELP: if set, prints out this help header
;
;    MFIRST: if set, the input array is a(m,l) instead of a(l,m)
;        where abs(m) <= l
;
; OUTPUTS:
;     index: Integer vector of size nl containing the index the 
;            of alm coefficients, related to {l,m} by 
;             i = l^2 + l + m + 1
;
;     alm_list: array of alm coefficients, with dimension (nl, 2 [, ns])
;            -- corresponding to
;                  nl: number of {l,m} indices
;                  2: for real and imaginary parts of alm coefficients
;                  ns: number of Stokes parameters (if > 1)
;
; RESTRICTIONS:
;     Unphysical values (m > l) are silently ignored
;
; PROCEDURE:
;
; EXAMPLE:
;    ; combining two different sets of alm for (l,m) <= 100
;    fits2alm, i1, a1, 'alm1.fits'                       ; read first set of alm from a FITS file
;    ac1 = alm_i2t(i1, a1, /complex, lmax=100, mmax=100) ; make an array out of those with (l,m) <= 100
;
;    fits2alm, i2, a2, 'alm2.fits'                       ; read second set of alm
;    ac2 = alm_i2t(i2, a2, /complex, lmax=100, mmax=100) ; make an array out of those with (l,m) <= 100
;
;    ac = 0.9*ac1 + 0.1*ac2          ; weighted sum the 2 alm sets 
;
;    alm_t2i, ac, i, a               ; makes an index list of the new alms
;    alm2fits, i, a, 'almsum.fits'   ; save the new alms into a FITS file
;
; MODIFICATION HISTORY:
;     2007-10-04: creation
;     2011-10-21: addition of MFIRST
;     2012-02-23: can deal with Stokes alm
;
;-

routine='alm_t2i'
syntax = strupcase(routine)+', alm_tab, index, alm_list [, /HELP, MFIRST=] '
if (keyword_set(help)) then begin
    doc_library,routine
    return
endif

if (n_params() eq 0) then begin
    print,syntax
    return
endif
if (n_params() ne 3) then begin
    print,syntax
    message,'Abort'
endif

; generates l and m arrays
sz = size(alm_tab)
if (sz[0] lt 1) then begin
    print,syntax
    print,'alm_tab must have at least 1 dimension'
    message,'Abort.'
endif
if (keyword_set(mfirst)) then begin ; m first index
    nm = sz[1]
    nl = (sz[0] gt 1) ? sz[2] : 1L
    mtab = lindgen(nm) # replicate(1,nl)
    ltab = replicate(1,nm) # lindgen(nl)
endif else begin
    nl = sz[1]
    nm = (sz[0] gt 1) ? sz[2] : 1L
    ltab = lindgen(nl) # replicate(1,nm)
    mtab = replicate(1,nl) # lindgen(nm)
endelse

if (nm gt nl) then begin
    print,syntax
    print,'It seems that Max m ('+strtrim(nm-1,2)+') is larger than Max l ('+strtrim(nl-1,2)+'), '
    print,'while we expect abs(m) <= l'
    message, 'Abort.'
endif

sample = alm_tab[0,0]*0 ; either (0.,0.) or (0.d0, 0.d0) or (0.) or (0.d0)
complex = (size(/tname, sample) eq 'COMPLEX' || size(/tname, sample) eq 'DCOMPLEX') 

; extra dimension (for eg, T, E, B alms)
de = (complex) ? 3 : 4 ; 3rd dimension if complex array, 4th dimension otherwize
nstokes = (sz[0] eq de) ? sz[de] : 1
if sz[0] gt de then begin
    print,syntax
    message,'Too many dimensions in input array'
endif

; find points where m <= l
k = where(mtab le ltab,nk)

; define index = l*(l+1)+m+1 for those points
index = ((ltab+1L)*ltab + mtab)[k] + 1L

; make sure that output alm are sorted by increasing index 
; (in which m is fastest varying)
if ~keyword_set(mfirst) then begin
    ii = sort(index)
    k = k[ ii ] 
    index = index [ ii ]
    ii = 0
endif

shift = nm*nl
if complex then begin
; fill the alm list with the real and imaginary parts
    alm_list = replicate(real_part(sample),nk,2, nstokes)
    for js=0, nstokes-1 do begin
        alm_list[0,0,js] = real_part(alm_tab[k + shift*js]) ; R(alm_tab[*,*,js])
        alm_list[0,1,js] = imaginary(alm_tab[k + shift*js]) ; I(alm_tab[*,*,js])
    endfor
endif else begin
; fill the alm list with the real and imaginary parts
    alm_list = replicate(sample,nk,2, nstokes)
    for js=0, nstokes-1 do begin
        alm_list[0,0,js] = alm_tab[k + shift*(2*js)  ] ; alm_tab[*,*,0,js]
        alm_list[0,1,js] = alm_tab[k + shift*(2*js+1)] ; alm_tab[*,*,1,js]
    endfor
endelse
if nstokes eq 1 then begin
    alm_list=reform(alm_list, nk, 2, /Overwrite)
endif



return
end

