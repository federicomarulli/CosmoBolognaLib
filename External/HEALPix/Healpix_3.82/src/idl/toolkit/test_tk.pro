pro detect_tk_error, nside, pixel, pixel2, error, legend=legend

if total(abs(pixel2-pixel)) ne 0 then begin
    np = n_elements(pixel)
    bad = where(pixel2 ne pixel, nbad)
    print, 'Nside = ',nside
    if (defined(legend)) then print, legend
    print, nbad,' / ',np,' pixels : '
    nbad <= 5
    print, 'In:   ',pixel[bad[0:nbad-1]]
    print, 'Out:  ',pixel2[bad[0:nbad-1]]
    print, 'Diff: ',pixel2[bad[0:nbad-1]]-pixel[bad[0:nbad-1]]
    error = error + 1
endif
return
end
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
pro test_tk, nside, upix, help=help, random=random
;+
; test_tk [,nside ,upix, help=, random= ]
;
;  tests the consistency of Healpix toolkit routines (pix2*, ang2*, vec2*, ang2vec, vec2ang)
;
;  by default, all pixels in 0, Npix-1 are tested (Npix=12*nside*nside)
;  unless upix or random are defined
;
;  OPTIONAL INPUT
;   nside (default = 32)
;   upix , integer scalar or vector, list of hand-picked pixels
;
;  KEYWORD
;   HELP=   : if set, this information header is printed
; 
;   RANDOM= : float scalar, some pixels (Npix * random + 1) are picked randomly
;   in [0, Npix-1]
;
;
;  EXAMPLES:
;   test_tk
;   test_tk, 1024
;   test_tk, 8192,  random=1.d-3
;   test_tk, 2L^29, random=1.d-12
;
;  HISTORY:
;    see test_tk2 for a slightly longer set of tests
;  2008-03-17: enabled Nside > 8192
;  2011-08:    use detect_tk_error
;  2017-06:    cosmetic editions
;-

t0 = systime(1)
routine = 'test_tk'
syntax = routine+' [nside, upix, HELP=, RANDOM=]'

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if undefined(nside) then nside = 32
lnside = long(nside)

npix = nside2npix(lnside,err=err_nside)
snpix = strtrim(string(npix,form='(i20)'),2)
snpix1 = strtrim(string(npix-1,form='(i20)'),2)
snside = strtrim(string(nside,form='(i9)'),2)

if (err_nside ne 0) then begin
    print, syntax
    print,'Invalid Nside'
    return
endif

if defined(random) then begin
    nr = min([npix, long(npix*random)+1])
    print,'Nside = '+snside
    print,nr,' pixels are picked randomly in [0, '+snpix1+']'
    if (nside gt 8192) then begin
        pixel =long64( randomu(seed,nr, /double) * npix )
    endif else begin
        pixel =long( randomu(seed,nr, /double) * npix )
    endelse
    print,min(pixel),max(pixel)
endif
if defined(upix) then begin
    if min(upix) ge 0 and max(upix) lt npix then begin
        pixel = upix
        print,'Nside = '+snside
        print,'test selected pixels in ',min(upix), max(upix)
    endif else begin
        print,'invalid choice of pixels'
        return
    endelse
endif
if undefined(pixel) then begin
    pixel = lindgen(npix)
    print,'Nside = '+snside
    print,'test all pixels in [0, '+snpix1+']'
endif

error = 0

;--------------------------
pix2ang_ring, nside, pixel, theta, phi
ang2pix_ring, nside, theta, phi, pixel2
detect_tk_error, nside, pixel, pixel2, error, legend='ERROR: pix <-> ang RING'


pix2ang_nest, nside, pixel, theta, phi
ang2pix_nest, nside, theta, phi, pixel2
detect_tk_error, nside, pixel, pixel2, error, legend='ERROR: pix <-> ang NEST'
;---------------------

pix2vec_ring, nside, pixel, vec
vec2pix_ring, nside, vec*!DPI, pixel2
detect_tk_error, nside, pixel, pixel2, error, legend='ERROR: pix <-> vec RING'

pix2vec_nest, nside, pixel, vec
vec2pix_nest, nside, vec*5.d0, pixel2
detect_tk_error, nside, pixel, pixel2, error, legend='ERROR: pix <-> vec NEST'
;---------------------

ring2nest, nside, pixel, pix_n
nest2ring, nside, pix_n, pixel2
detect_tk_error, nside, pixel, pixel2, error, legend='ERROR:   NEST <-> RING'
;---------------------


pix2ang_ring, nside, pixel, theta, phi
ang2vec,             theta, phi, vec
vec2pix_nest, nside, vec, pixel1
nest2ring,    nside, pixel1, pixel2
detect_tk_error, nside, pixel, pixel2, error, legend='ERROR:   R -> A -> V -> N -> R'
;---------------------

pix2vec_ring, nside, pixel, vec
vec2ang,             vec, theta, phi
ang2pix_nest, nside, theta, phi, pixel1
nest2ring,    nside, pixel1, pixel2
detect_tk_error, nside, pixel, pixel2, error, legend='ERROR:   R -> V -> A -> N -> R'
;---------------------
t1 = systime(1)

if (error eq 0) then begin
    print, ' -----------------------'
    print, '   TEST PASSED ('+snside+')'
    print, ' -----------------------'
    print, 'Time [s]: ', t1-t0
endif

return
end
