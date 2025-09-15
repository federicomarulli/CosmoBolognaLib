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
PRO ang2pix_nest, nside, theta, phi, ipnest
;*******************************************************************
;+
; ANG2PIX_NEST, Nside, Theta, Phi, Ipix
;
;        renders the NESTED scheme pixel number Ipix for a pixel which, given the
;        map resolution parameter Nside, contains the point on the sphere
;        at angular coordinates Theta and Phi
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;       should be a power of 2 (not tested)
;	SCALAR
;    Theta : angle (along meridian), in [0,Pi], theta=0 : north pole,
;	can be an ARRAY
;    Phi   : angle (along parallel), in [0,2*Pi]
;	can be an ARRAY of same size as Theta
;
; OUTPUT
;    Ipix  : pixel number in the NESTED scheme HEALPIX pixelisation in [0,Npix-1]
;	is an ARRAY of same size as Theta and Phi
;
; SUBROUTINE
;    nside2npix
;
; HISTORY
;    June-October 1997,  Eric Hivon & Kris Gorski, TAC
;                  original ang_pix_nest
;    Feb 1999,           Eric Hivon, Caltech
;                  renamed ang2pix_nest
;    March 1999,  correction of a bug when nside is not LONG
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Aug  2004,  EH, use !PI as theta upper-bound instead of !DPI
;    Dec 2007, EH, enabled Nside > 8192
;    Aug 2011, EH, IAP: more accurate calculations close to pole
;
;-
;*****************************************************************************

  routine = 'ANG2PIX_NEST'
  if N_params() ne 4 then begin
      message,' syntax: '+routine+', nside, theta, phi, ipix'
  endif

  if (N_ELEMENTS(nside) GT 1) then message, 'Nside should be a scalar'
  npix = nside2npix(nside, error = error)
  if (error ne 0) then message, 'Invalid Nside: '+string(nside)

  np = N_ELEMENTS(theta)
  np1 = N_ELEMENTS(phi) 
  if (np NE np1) then begin
      message,'inconsistent theta and phi'
  endif
  if (MIN(theta) lt 0.) or (MAX(theta) gt !PI) then message, 'theta out of range'
;------------------------------------------------------------

  nl1 = LONG(nside[0]) ; 2018-06-06
  ns_max = nl1
  l64 = 0 ; all in 32 bit
  if (nl1 gt 8192) then begin
      ix = LONARR(np, /NoZero) 
      iy = LONARR(np, /Nozero)
      npface = nl1 * long64(nl1)
      l64 = (nl1 ge ishft(1L,29)) ; for Nside >= 2^29
  endif else begin
      ix = INTARR(np, /NoZero)
      iy = INTARR(np, /Nozero)
      npface = nl1 * (nl1)
  endelse
  face_num = BYTARR(np, /NoZero)

  common xy2pix, x2pix, y2pix
  if (n_elements(x2pix) eq 0) then init_xy2pix ; initiate x2pix and y2pix

  twopi  = 2.d0 * !DPI
  piover2 = !DPI / 2.d0
  z0 = 2.d0 /3.d0
  z = COS(DOUBLE(theta))
  phi_in = phi MOD twopi
  phi_in += (phi_in lt 0.0  )*twopi
  tt  = phi_in / piover2        ; in [0,4[
  phi_in = 0

  pix_eqt = WHERE(z LE z0 AND z GT -z0, n_eqt) ; equatorial strip
  IF (n_eqt GT 0) THEN BEGIN

;     (the index of edge lines increase when the longitude=phi goes up)
;       jp = LONG(ns_max*(0.5d0 + tt(pix_eqt) - z(pix_eqt)*0.75d0)) ;  ascend edge line index
;       jm = LONG(ns_max*(0.5d0 + tt(pix_eqt) + z(pix_eqt)*0.75d0)) ; descend edge line index
      jp = floor(ns_max*(0.5d0 + tt(pix_eqt) - z(pix_eqt)*0.75d0), L64=l64) ;  ascend edge line index
      jm = floor(ns_max*(0.5d0 + tt(pix_eqt) + z(pix_eqt)*0.75d0), L64=l64) ; descend edge line index

;     finds the face
      face_n = BYTARR(n_eqt)
      ifp = BYTE(jp / ns_max)   ; in {0,4}
      ifm = BYTE(jm / ns_max)
      p_np = WHERE(ifp EQ ifm, n_np)
      p_eq = WHERE(ifp LT ifm, n_eq)
      p_sp = WHERE(ifp GT ifm, n_sp)
      if (n_np GT 0) then face_n(p_np) = (ifp(p_np) MOD 4) + 4
      if (n_eq GT 0) then face_n(p_eq) = (ifp(p_eq) MOD 4)
      if (n_sp GT 0) then face_n(p_sp) = (ifm(p_sp) MOD 4) + 8

      face_num(pix_eqt) = face_n
      ix(pix_eqt) = jm MOD ns_max
      iy(pix_eqt) = ns_max - (jp MOD ns_max) - 1

      jp=0 & jm=0 & face_n=0 & ifp=0 & ifm=0 & face_n=0 & p_np=0 & p_eq=0 & p_sp=0
      pix_eqt = 0
  ENDIF

  pix_pol = WHERE(z GT z0 OR z LE -z0, n_pol) ; polar caps
  IF (n_pol GT 0) THEN BEGIN


      ;tmp = SQRT( 3.d0*(1.d0 - ABS(z(pix_pol))) ) ; in ]0,1]
      angle = double(theta[pix_pol])
      angle <= !DPI - angle ; either theta or Pi-theta, whichever is smaller
      tmp = sqrt(6.d0) * SIN(angle * 0.5d0) ; more accurate close to pole
      angle = 0

      ntt = FIX(tt(pix_pol)) < 3
      tp = tt(pix_pol) - ntt

;     (the index of edge lines increase when distance from the closest pole goes up)
      jp = LONG( ns_max * tp          * tmp ) ; line going toward the pole as phi increases
      jm = LONG( ns_max * (1.d0 - tp) * tmp ) ; that one goes away of the closest pole
      jp = jp < (ns_max-1)      ; for points too close to the boundary
      jm = jm < (ns_max-1)

;     finds the face and pixel's (x,y)
      p_sp = WHERE(z[pix_pol] lt 0., n_sp, complement=p_np, ncomplement=n_np)
      if (n_np GT 0) then begin
          face_num(pix_pol(p_np)) = ntt(p_np)
          ix(pix_pol(p_np)) = ns_max - jm(p_np) - 1
          iy(pix_pol(p_np)) = ns_max - jp(p_np) - 1
      endif
      if (n_sp GT 0) then begin
          face_num(pix_pol(p_sp)) = ntt(p_sp) + 8
          ix(pix_pol(p_sp)) = jp(p_sp)
          iy(pix_pol(p_sp)) = jm(p_sp)
      endif

      ntt=0 & tp=0 & tmp=0 & jp=0 & jm=0 & p_np=0 & p_sp=0
      pix_pol = 0
  ENDIF

  if (nl1 gt 8192) then begin
      smax = 4
      scale = 1LL
      ipnest = 0LL
      scale_factor = 128LL*128LL
  endif else begin
      smax = 1
      scale = 1L
      ipnest = 0L
      scale_factor = 128L*128L
  endelse      

  for i=0, smax-1 do begin
      ix_low = ix and 127 ; last 7 bits
      iy_low = iy and 127 ; last 7 bits
      ipnest += (x2pix[ix_low] + y2pix[iy_low]) * scale
      scale *= scale_factor
      ix /= 128                 ; truncate out last 7 bits
      iy /= 128 
  endfor
  ipnest += (x2pix[ix] + y2pix[iy]) * scale
  ix_low = 0 & ix = 0 & iy_low = 0 & iy = 0

  ipnest += face_num* npface ; in {0, 12*nside**2 - 1}


  RETURN
END                             ; ang2pix_nest
