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
pro bin_llcl, llcl_in, ubin, l_out, llcl_out, dllcl, uniform=uniform, deltal=deltal, flatten=flatten, help=help
;+
; NAME:
;   bin_llcl
;
; PURPOSE:
;   turns llcl_in (= continuous l*(l+1)*cl/2Pi) into a binned version
;    with a bin BIN
;
; CALLING SEQUENCE:
;   bin_llcl, llcl_in, bin, l_out, llcl_out, [dllcl, deltal=, /Flatten, /Help, /Uniform]
;
;
; INPUTS
;   llcl_in : input l*(l+1)*Cl/2Pi, 1D vector, defined for each l from l=0
;   bin : can be either a scalar = dl
;     or a vector defining the bins boundaries : 
;     [low0, low1, low2, ...,low(n-2), low(n-1)+1]
;
; OUTPUTS
;     l_out: center of bins
;     llcl_out: binned l*(l+1)*Cl/2Pi
;
; OPTIONAL OUTPUTS
;   dllcl : returns on output the rms of C(l) for a full sky
;   observation
;     = C(l) * sqrt( 2/ 2l+1 / dl)
;
;  deltal = width of each bin
;
; KEYWORDS
;   help: if set, print this documentation
;
;   flatten: if set, the input C(l) is multiplied internally by l*(l+1)/2Pi before being
;     binned. By default, the input C(l) is binned as is.
;
;   uniform: if set, each l is given the same weight in the bin.
;     By default, a weighting propto (2*l+1) (inverse cosmic
;     variance) is applied to each l.
;     In any cases, the l_out is the same
;
; HISTORY
;  ?? : Eric Hivon, first version
;  2006-oct: shortens bins list on the fly
;  2007-04-05: added deltal optional output
;  2007-05-31: added /flatten and /help
;  2008-03-03: ignore bins beyond data lmax, truncate last valid bin if necessary.
;  2013-03-15: correctly deals with optional Dllcl argument
;  2014-07-21: correction of a mistyped call to MESSAGE mostly showing up under GDL
;  2016-11-23: typo corrections in this header
;-


routine = 'bin_llcl'
uroutine = strupcase(routine)
usage = ['usage: '+uroutine+', llcl_in, bin, l_out, llcl_out, [dllcl, Uniform=, Deltal=, Flatten=, Help=]']

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if n_params() eq 0 then begin
    print,usage,form='(a)'
    print,'Type '+routine+', /help              for extended help'
    return
endif
if (n_params() lt 4 || n_params() gt 5) then begin
    print,'Number of arguments found: ',n_params()
    print,'Expected 4 or 5'
    message,usage,/noname
endif
if (min(ubin) lt 0) then begin
    print,'bins must be positive'
    message,usage,/noname
endif

bin = ubin
nb = n_elements(bin)
lmax_in = n_elements(llcl_in[*,0])-1
if (nb gt 1) then begin
    k = where(bin le (lmax_in+1),nk)
    if (nk eq (nb-1)) then bin = bin < (lmax_in+1) ; shorten last bin
    if (nk lt (nb-1)) then begin ; shorten last valid bins, and drop the ones beyond lmax
        bin = [bin[k], lmax_in+1]
        nb = nk + 1
;        print,'BIN_LLCL: too many bins'
;        print,nk,' bins available',nb,' bins given'
;        return
    endif
endif

if nb eq 1 then begin
    ; regular binning
    nbins = lmax_in/long(bin)
    lmax = nbins * bin-1
    l = findgen(lmax+1)
    w = 2*l + 1
    if keyword_set(uniform) then w = replicate(1.,lmax+1)
    y = llcl_in[0:lmax, 0]
    if keyword_set(flatten) then y = y * (l*(l+1.)/(2*!pi))
    w1 = reform(w,    bin, nbins)
    y1 = reform(y*w,  bin, nbins)
    l1 = reform(l,    bin, nbins)
    n1 = replicate(1, bin, nbins)
    
    llcl_out = total(y1,1)/total(w1,1)
    l_out    = total(l1,1)/total(n1,1)
    dl       = replicate(bin,nbins)
endif else begin
    ; irregular binning
    lmax = max(bin)-1
    nbins = nb-1
    good = where(bin lt lmax, ng)
    if (ng eq 0)  then begin
        message,/info,'l-range of binning does not intersect that of data'
        return
    endif
    l = findgen(lmax+1)
    n1 = replicate(1,lmax+1)
    w = 2*l + 1
    if keyword_set(uniform) then w = replicate(1.,lmax+1)
    y = llcl_in[0:lmax, 0]
    if keyword_set(flatten) then y = y * (l*(l+1.)/(2*!pi))
    l_out    = fltarr(nbins)
    llcl_out = fltarr(nbins)
    dl       = intarr(nbins)
    for i=0,nbins-1 do begin
        l_out[i] = mean( l[bin[i]:bin[i+1]-1])
        dl[i]    = bin[i+1]-bin[i]
        llcl_out[i] = total( (y*w) [bin[i]:bin[i+1]-1] ) $
                    / total(    w  [bin[i]:bin[i+1]-1] )
    endfor
endelse

dllcl = llcl_out * sqrt(2/(2*l_out+1.)/dl)

deltal = dl

return
end

