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
pro test_tk2, nside, upix, random=random
;+
; test_tk2 [,nside ,upix, random= ]
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
;   random : float scalar, some pixels (random*Npix +1) are picked randomly in [0, Npix-1]
;
;  EXAMPLES:
;   test_tk2
;   test_tk2, 1024
;   test_tk2, 8192,  random=1.d-3
;   test_tk2, 2L^29, random=1.d-12
;
;  HISTORY:
;   based on test_tk, with a few extra tests on vectors
;  2011-08: cosmetic editions
;  2017-06: more cosmetic editions
;-

t0 = systime(1)
if undefined(nside) then nside = 32
;lnside = long(nside)

trigger = (machar(/double)).eps * 8 ; 1.77e-15

;npix = nside2npix(lnside,err=err_nside)
npix = nside2npix(nside,err=err_nside)
snpix = strtrim(string(npix,form='(i20)'),2)
snpix1 = strtrim(string(npix-1,form='(i20)'),2)
snside = strtrim(string(nside,form='(i9)'),2)

if (err_nside ne 0) then begin
    print,'Invalid Nside'
    return
endif

if defined(random) then begin
    nr = min([npix, long(random*npix)+1])
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
np = n_elements(pixel)

error = 0

;--------------------------
pix2ang_ring, nside, pixel, theta, phi
ang2pix_ring, nside, theta, phi, pixel2
if total(abs(pixel2-pixel)) ne 0 then begin
    print,'error pix <-> ang ring', nside
    print,'out pixel range: ',minmax(pixel2)
    error = error + 1
endif

pix2ang_nest, nside, pixel, theta, phi
ang2pix_nest, nside, theta, phi, pixel2
if total(abs(pixel2-pixel)) ne 0 then begin
    print,'error pix <-> ang nest', nside
    print,'out pixel range: ',minmax(pixel2)
    error = error + 1
endif
;---------------------

pix2vec_ring, nside, pixel, vec
vec2pix_ring, nside, vec, pixel2
if total(abs(pixel2-pixel)) ne 0 then begin
    print,'error pix <-> vec ring', nside
    error = error + 1
    help,pixel, pixel2
    bad = where( (pixel2-pixel) ne 0, nbad)
    print,'#bad:',nbad
    print,'in:',pixel[bad]
    print,'out:',pixel2[bad]
    print,'diff:',(pixel2-pixel)[bad]
endif

pix2vec_nest, nside, pixel, vec
vec2pix_nest, nside, vec, pixel2
if total(abs(pixel2-pixel)) ne 0 then begin
    print,'error pix <-> vec nest', nside
    print,'out pixel range: ',minmax(pixel2)
    error = error + 1
    bad = where( (pixel2-pixel) ne 0, nbad)
    print,'#bad:',nbad
;     print,pixel
;     print,pixel-pixel2
endif
;---------------------

ring2nest, nside, pixel, pix_n
nest2ring, nside, pix_n, pixel2
if total(abs(pixel2-pixel)) ne 0 then begin
    print,'error ring <-> nest', nside
    error = error + 1
endif
;---------------------


pix2ang_ring, nside, pixel, theta, phi
ang2vec,             theta, phi, vec
vec2pix_nest, nside, vec, pixel1
nest2ring,    nside, pixel1, pixel2
if total(abs(pixel2-pixel)) ne 0 then begin
    ;print,'error loop1', nside
    print,'R -> A -> V -> N -> R', nside
    error = error + 1
endif
;---------------------

pix2vec_ring, nside, pixel, vec
vec2ang,             vec, theta, phi
ang2pix_nest, nside, theta, phi, pixel1
nest2ring,    nside, pixel1, pixel2
if total(abs(pixel2-pixel)) ne 0 then begin
    ;print,'error loop2', nside
    print,'R -> V -> A -> N -> R', nside
    error = error + 1
endif
;---------------------

pix2vec_ring, nside, pixel, vec, vv
ring2nest, nside, pixel, pixnest
pix2vec_nest, nside, pixnest, vec2, vv2
t1 = total(abs(vec-vec2))/np
m1 = max(abs(vec-vec2))
if  (t1 gt trigger || m1 gt trigger) then begin
    ;print,'error loop3 vector', nside
    print,'R -> V  vs  R -> N -> V,  vector', nside
    print,t1,m1,n_elements(vec)
    print,total(vec,1)/np
    error = error + 1
endif
t2 = total(abs(vv-vv2))/np
m2 = max(abs(vv-vv2))
if (t2 gt trigger || m2 gt trigger) then begin
    ;print,'error loop3 vertex', nside
    print,'R -> V  vs  R -> N -> V,  vertex', nside
    print,t2,m2,n_elements(vv)
    print,total(total(vv,1),2)/np
    error = error + 1
endif

t1 = systime(1)
if (error eq 0) then begin
    print, ' -----------------'
    print, '   TEST2 PASSED ('+snside+')'
    print, ' -----------------'
    print, 'Time [s]: ', t1-t0
endif


return
end
