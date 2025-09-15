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
pro coadd_raw_obs
;+
; coadd_raw_obs : coadd raw observation into a signal and n_obs per
; pixel maps, put them in a FITS file, see them in Mollweide
; projection
;
; Nov 2006: replaced GIF output by PNG, 
;           use ang2vec
;           use more functionalities of write_fits_cut4
;-

; reads the table of pointings in Equat. 2000 coord.
file_in = 'raw_obs.txt'
npoints = 3000l
tmp = fltarr(3,npoints)
openr,1,file_in
readf,1,tmp
close,1
ratab = tmp(0,*)
dectab = tmp(1,*)
measure = tmp(2,*)


; array of unit vector in Equat. 2000 coordinate
ang2vec, dectab, ratab, vec_eq, /astro

; conversion  eQuatorial (=celestial) 2000 --> Galactic
loadsky                         ; cgis package routine, define rotation matrices
outcoord = 'G'                  ; Galactic coord
vec_gal = SKYCONV(vec_eq,inco='Q',outco=outcoord)  ; cgis package function

; conversion (x,y,z) Galactic --> pixel number
nside = 32 ; resolution 120'/pixel
npix = nside2npix(nside)
VEC2PIX_NEST, nside, vec_gal, ipnest   ;  <<<<<<< nested scheme (see below)

; makes maps (counts the number of hits per pixels and coadd signal)
n_obs = lonarr(npix)
n_obs(0:npix-1) = ( HISTOGRAM( ipnest, MIN = 0l, MAX= npix+1, REV = rr) )(0:npix-1)
print,total(n_obs)

pixel = where(n_obs gt 1)

signal = fltarr(npix)
for i=0l, npoints-1 do signal(ipnest(i)) = signal(ipnest(i)) + measure(i)
signal = signal / ( n_obs > 1)

signal = signal[pixel]
n_obs = n_obs[pixel]

; ----------- output the hit map in a FITS file ------------------
file_fits = 'coadded.fits'


;-------------------------------------
; extension unit
info_xhdr = STRARR(1)
SXADDPAR,info_xhdr,'HISTORY','-----------------------------------------------------'
SXADDPAR,info_xhdr,'HISTORY','file created by COADD_RAW_OBS from '+file_in
SXADDPAR,info_xhdr,'HISTORY','-----------------------------------------------------'
SXADDPAR,info_xhdr,'BAD_DATA',-1.6375E30,' value for missing data',form='(e15.8)'
units='mV'
ordering = 'NESTED'  ; <<<<<< nested scheme (see above)

; creates the structure containing the binary extension header and
; data columns

write_fits_cut4, file_fits, pixel, signal, n_obs, xh=info_xhdr, units=units, nside=nside, coord=outcoord, ordering=ordering, ext=0

; look at the map in Mollweide projection and makes a png file of it
select='signal'
mollview,file_fits,PNG=select+'.png'

select='n_obs'
mollview,file_fits,PNG=select+'.png',select

return
end

