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
pro vec2ang, vector, theta, phi, astro = astro
;*****************************************************************************
;+
; NAME:
;       VEC2ANG
;
; PURPOSE:
;       converts  cartesian coordinates vector Vector 
;       into angular coordinates Theta, Phi
;
; CALLING SEQUENCE:
;       Vec2ang, Vector, Theta, Phi, [/Astro]
;
; INPUTS:
;       Vector : array of cartesian position of dimension (np, 3)
;        with north pole = (x=0,y=0,z=1) and origin at the center of
;        the sphere
;       and vector = (x(0), x(1), ... x(np-1), y(0), y(1), ... y(np-1),
;         z(0), z(1), ...z(np-1))
;
; KEYWORD PARAMETERS:
;       Astro : see below
;
; OUTPUTS:
;       Theta : vector of angular coordinate of size np
;       Phi   : vector of angular coordinate of size np
;
;
;        * if ASTRO is NOT set (default) : geometric coordinates
;       Theta : colatitude in RADIAN measured Southward from North pole
;       Phi   : longitude in RADIAN, increasing Eastward
;
;        * if ASTRO is set : astronomic coordinates
;       Theta : latitude in DEGREE measured Northward from Equator
;       Phi   : longitude in DEGREE, increasing Eastward
;
; MODIFICATION HISTORY:
;       March 6, 1999    Eric Hivon, Caltech, Version 1.0
;       March 22, 2002     use Message
;       Oct 4, 2007      Check dimensions of Vector
;       Aug 2011, IAP: more precise theta determination close to pole
;-
;*****************************************************************************

routine = 'VEC2ANG'
if N_params() lt 2 then begin
    message,' syntax = vec2ang, vec, theta, phi, astro=',/noprefix
endif

np1 = N_ELEMENTS(vector)
np = np1 / 3L
if (np1 NE 3 * np) then begin
    message,/info,'input Vector should have Npix*3 elements',/noprefix
    message,'It currently has: '+strtrim(np1,2),/noprefix
endif
sz = size(vector)
ndim = sz[0] ; number of dimension
if (sz[ndim] ne 1 and sz[ndim] ne 3 and sz[ndim] ne np1) then begin
; if (sz[ndim] ne 3) then begin
    dimensions = strtrim(sz[1:ndim],2)+[replicate('x',ndim-1),'']
    message,/info,'Input Vector has dimensions: ',/noprefix
    print,routine+':         ', dimensions
;    message,/info,'I expected: '+strtrim(np,2)+'x 3 or
;    '+strtrim(np1,2),/noprefix
    message,/info,'I expected: '+strtrim(np,2)+'x 3',/noprefix
    print,routine+' x0, x1, x2, ...'
    print,routine+' y0, y1, y2, ...'
    print,routine+' z0, z1, z2, ...'
    message,'some reordering of the array may be necessary',/noprefix
endif

twopi = 2. * !Pi
radeg = !RaDeg
zero = 0.
ninety = 90.
if (size(/TYPE,vector) EQ 5) then begin ; double precision
    twopi = 2.d0 *!DPi 
    radeg = 180.d0 / !DPi
    zero = 0.d0
    ninety = 90.d0
endif

;---------------

Vector = REFORM(vector, np, 3, /OVER) ; condition the input vector
;theta_rad = ACOS( vector(*,2)  /  SQRT(  TOTAL(vector^2, 2) ) )
theta_rad = ATAN( SQRT(vector[*,0]^2+vector[*,1]^2)  , vector[*,2]) ; in [0,Pi]
phi_rad   = ATAN( vector[*,1], vector[*,0] )  ; in [-Pi,Pi]
vector = reform(vector,sz[1:ndim], /over) ; revert to original shape

phi_rad = phi_rad + twopi * (phi_rad LT zero)

IF KEYWORD_SET(astro) THEN BEGIN
    theta = ninety - theta_rad * RaDeg
    phi   = phi_rad * RaDeg
ENDIF ELSE BEGIN
    theta = theta_rad
    phi = phi_rad
ENDELSE

theta_rad = 0
phi_rad = 0


return
end

