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
pro ang2vec, theta, phi, vector, astro = astro
;*****************************************************************************
;+
; NAME:
;       ANG2VEC
;
; PURPOSE:
;       converts angular coordinates Theta, Phi into cartesian
;       coordinates vector Vector
;
; CALLING SEQUENCE:
;       Ang2vec, Theta, Phi, Vector, [/Astro]
;
; INPUTS:
;       Theta : vector or array of angular coordinate of size np
;       Phi   : vector or array of angular coordinate of size np
;
;        * if ASTRO is NOT set (default) : geometric coordinates
;       Theta : colatitude in RADIAN measured Southward from North pole
;       Phi   : longitude in RADIAN, increasing Eastward
;
;        * if ASTRO is set : astronomic coordinates
;       Theta : latitude (or Declination) in DEGREE measured Northward from Equator
;       Phi   : longitude (or Right Ascension) in DEGREE, increasing Eastward
;
; KEYWORD PARAMETERS:
;       Astro : see above
;
; OUTPUTS:
;       Vector : array of normalised cartesian position of dimension (np, 3)
;        with north pole = (x=0,y=0,z=1) and origin at the center of
;        the sphere
;       and vector = (x(0), x(1), ... x(np-1), y(0), y(1), ... y(np-1), z(0), z(1), ...z(np-1))
;
; MODIFICATION HISTORY:
;       March 6, 1999    Eric Hivon, Caltech, Version 1.0
;       March 22, 2002     use Message
;       Aug 2011, IAP: more precise (x,y) determination close to pole
;       May 2017: replaced obsolete DATATYPE() with IDL's SIZE(/TNAME)
;
;-
;*****************************************************************************

if N_params() ne 3 then begin
    message,' syntax = ang2vec, theta, phi, vec, astro=',/noprefix
endif

np = N_ELEMENTS(theta)
np1 = N_ELEMENTS(phi) 
if (np NE np1) then begin
    message,'inconsistent theta and phi in ang2vec',/noprefix
endif

dtor = !DtoR
if (size(theta,/type) EQ 5 || size(phi,/type) EQ 5) then begin ; double precision
    dtor = !DPi / 180.d0
endif

;---------------

IF KEYWORD_SET(astro) THEN BEGIN
    theta_rad = (90. - theta) * DtoR
    phi_rad   = phi * DtoR
ENDIF ELSE BEGIN
    theta_rad = theta
    phi_rad = phi
ENDELSE

max_theta = MAX(theta_rad, MIN = min_theta)
IF (max_theta GT !Pi   OR  min_theta LT 0) THEN  BEGIN
    message,' theta out of range in ang2vec',/noprefix
ENDIF

if (np gt 1) then theta_rad = REFORM(theta_rad,np,/OVER)
if (np gt 1) then phi_rad   = REFORM(phi_rad,np,/OVER)

cos_theta = COS(theta_rad)
;sin_theta = SQRT( (1.-cos_theta)*(1.+cos_theta) )
sin_theta = SIN(theta_rad)
theta_rad = 0

Vector = [ [sin_theta * COS(phi_rad)], [sin_theta*SIN(phi_rad)], [cos_theta] ]
phi_rad = 0


return
end

