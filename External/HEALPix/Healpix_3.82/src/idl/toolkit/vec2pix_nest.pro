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
PRO vec2pix_nest, nside, vec_in, ipnest
;*******************************************************************
;+
; VEC2PIX_NEST, Nside, Vec_in, Ipix
;
;        renders the NESTED scheme pixel number Ipix for a pixel which, given the
;        map resolution parameter Nside, contains the point on the sphere
;        at cartesian coordinates Vec_in
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;            should be a power of 2 (not tested)
;	     SCALAR
;    Vec_in    : (x,y,z) position unit vector(s) with North pole = (0,0,1)
;            stored as x(0), x(1), ..., y(0), y(1), ..., z(0), z(1) ..
;            ARRAY of dimension = (np,3)
;
; OUTPUT
;    Ipix  : pixel number in the RING scheme of HEALPIX pixelisation in [0,Npix-1]
;	is a VECTOR of dimension = (np)
;
; SUBROUTINE
;    nside2npix
;
; HISTORY
;    June-October 1997,  Eric Hivon & Kris Gorski, TAC
;    Feb 1999,           Eric Hivon,               Caltech
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Dec 2007, EH, enabled Nside > 8192
;    Aug 2011, EH, IAP: more accurate calculations close to pole
;
;
;-
;*****************************************************************************

  routine = 'VEC2PIX_NEST'
  if N_params() ne 3 then begin
      message,' syntax: '+routine+', nside, vec, ipix'
  endif

  if (N_ELEMENTS(nside) GT 1) then message,'Nside should be a scalar in '+routine
  npix = nside2npix(nside, error = error)
  if (error ne 0) then message,'Invalid Nside '+string(nside)


  np1 = N_ELEMENTS(vec_in)
  np = LONG(np1/3) 
  if (np1 NE np*3) then begin
      print,'inconsistent vec_in in '+routine
      print,routine+', nside, vec_in, ipix'
      message,'Abort'
  endif
  vec_in = reform(vec_in,np,3,/OVERWRITE)
;------------------------------------------------------------

  nl1 = LONG(nside[0]) ; 2018-06-06
  ns_max = nl1
  l64 = 0
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
  norm = SQRT( TOTAL(vec_in^2, 2, /double) )
  z    = vec_in[*,2] / norm
  sz   = SQRT( vec_in[*,0]^2 + vec_in[*,1]^2 ) / norm
  norm = 0
  phi_in = ATAN(vec_in(*,1), vec_in(*,0))
  phi_in += (phi_in lt 0.0  )*twopi
  tt  = phi_in / piover2        ; in [0,4[
  phi_in = 0

  pix_eqt = WHERE(z LE z0 AND z GT -z0, n_eqt, complement=pix_pol, ncomplement=n_pol) ; equatorial strip
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
      if (n_np GT 0) then face_n(p_np) = (ifp(p_np) MOD 4B) + 4B
      if (n_eq GT 0) then face_n(p_eq) = (ifp(p_eq) MOD 4B)
      if (n_sp GT 0) then face_n(p_sp) = (ifm(p_sp) MOD 4B) + 8B

      face_num(pix_eqt) = face_n
      ix(pix_eqt) = jm MOD ns_max
      iy(pix_eqt) = ns_max - (jp MOD ns_max) - 1

      jp=0 & jm=0 & face_n=0 & ifp=0 & ifm=0 & face_n=0 & p_np=0 & p_eq=0 & p_sp=0
      pix_eqt = 0
  ENDIF

;;;;;;  pix_pol = WHERE(z GT z0 OR z LE -z0, n_pol) ; polar caps
  IF (n_pol GT 0) THEN BEGIN

      zz = z(pix_pol)
      ntt = FIX(tt(pix_pol)) < 3
      tp = tt(pix_pol) - ntt
;       tmp = SQRT( 3.d0*(1.d0 - ABS(z(pix_pol))) ) ; in ]0,1]
      tmp = sz[pix_pol] * SQRT(3.d0 / (1.d0 + abs(z[pix_pol]))) ; more accurate

;     (the index of edge lines increase when distance from the closest pole goes up)
      jp = LONG( ns_max * tp          * tmp ) ; line going toward the pole as phi increases
      jm = LONG( ns_max * (1.d0 - tp) * tmp ) ; that one goes away of the closest pole
      jp = jp < (ns_max-1)      ; for points too close to the boundary
      jm = jm < (ns_max-1)

;     finds the face and pixel's (x,y)
      p_np = WHERE(zz gt 0., n_np)
      p_sp = WHERE(zz lt 0., n_sp)
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

  z = 0 & sz = 0

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
END                             ; vec2pix_nest

