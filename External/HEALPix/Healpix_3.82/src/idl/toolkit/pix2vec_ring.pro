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
PRO pix2vec_ring, nside, ipix, vec_out, vertex
;*******************************************************************************
;+
; PIX2VEC_RING, Nside, Ipix, Vec_out [,Vertex]
; 
;       renders cartesian coordinates Vec_out of the nominal pixel center
;       given the RING scheme pixel number Ipix and map resolution parameter Nside
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;	SCALAR
;    Ipix  : pixel number in the RING scheme of Healpix pixelisation in [0,Npix-1]
;	can be an ARRAY of size (np) 
;       pixels are numbered along parallels (ascending phi), 
;       and parallels are numbered from north pole to south pole (ascending theta)
;
; OUTPUT
;    Vec_out : (x,y,z) position unit vector(s) of pixel center with North pole = (0,0,1)
;       stored as x(0), x(1), ..., y(0), y(1), ..., z(0), z(1) ..
;       is an ARRAY of dimension (np,3)
;
; OPTIONAL OUTPUT
;    Vertex : (x,y,z) position unit vector of the 4 vertices of each pixel
;      stored as x_N(0),x_N(1), ... y_N(0), y_N(1), ... z_N(0),z_N(1), ...
;                x_W(0),x_W(1), ... y_W(0), y_W(1), ... z_W(0),z_W(1), ...
;                x_S(0),x_S(1), ... y_S(0), y_S(1), ... z_S(0),z_S(1), ...
;                x_E(0),x_E(1), ... y_E(0), y_E(1), ... z_E(0),z_E(1), ...
;        ( where N: northrn vertex, W: western vertex, ... )
;      is an ARRAY of dimension (np,3,4)
;
; SUBROUTINE
;    nside2npix
;    sub_compute_vertices
;
; HISTORY
;    June-October 1997,  Eric Hivon & Kris Gorski, TAC : pix_ang
;    Feb 1999,           Eric Hivon,               Caltech
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    Nov  2002,          EH
;           added vertex output
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Dec 2007, EH,  IAP, enabled nside > 8192
;    Aug 2008, EH, IAP: issues warning if ipix is not of integer type
;    Aug 2011, EH, IAP: uses CHEAP_ISQRT to compute ring index out of pixel index
;     more precise (x,y) determination close to pole (for center only, not vertex)
;
;-
;*******************************************************************************
  routine = 'PIX2VEC_RING'
  if N_params() lt 3 then begin
      message,' syntax: pix2vec_ring, nside, ipix, vec [,vertex]'
  endif

  if (N_ELEMENTS(nside) GT 1) then message,'Nside should be a scalar in '+routine
  npix = nside2npix(nside, error = error)
  if (error ne 0) then message,'Invalid Nside '+string(nside)

  assert_pixindex_type, ipix[0],/warning ; warning if ipix is not integer
  min_pix = MIN(ipix, MAX = max_pix)
  IF (min_pix LT 0) THEN BEGIN
      PRINT,'pixel index : ',min_pix,FORMAT='(A,I18)'
      PRINT,'is out of range : ',0,npix-1,FORMAT='(A,I2,I18)'
      message,'Abort'
  ENDIF
  IF (max_pix GT npix-1) THEN BEGIN
      PRINT,'pixel index : ',max_pix,FORMAT='(A,I18)'
      PRINT,'is out of range : ',0,npix-1,FORMAT='(A,I2,I18)'
      message,'Abort'
  ENDIF

  do_vertex = (N_params() eq 4)


  nl1 = LONG(nside[0]) ; 2018-06-06
  nl2 = 2L*nl1
  if (nl1 gt 8192) then begin
      nl4 = 4LL*nl1
      ncap = nl2*(nl1-1LL)
      nsup = nl2*(5LL*nl1+1LL)
      one = 1LL
      four = 4LL
      l64 = 1
  endif else begin
      nl4 = 4L*nl1
      ncap = nl2*(nl1-1L)
      nsup = nl2*(5L*nl1+1L)
      one = 1L
      four = 4L
      l64 = 0
  endelse
  fact1 = 1.5d0*nl1
  fact2 = (3.d0*nl1)*nl1
  ifact2 = 1.d0 / ((3.d0*nl1)*nl1)
  np = N_ELEMENTS(ipix)
  vec_out = DBLARR(np,3, /NoZero)
  if (do_vertex) then vertex = dblarr(np,3,4, /NoZero)
  halfpi = !dpi / 2.0d0

  pix_np = WHERE(ipix LT ncap,   n_np) ; north polar cap
  IF (n_np GT 0) THEN BEGIN     ; north polar cap ; ---------------------------------

      ip = ROUND(ipix[pix_np], L64=l64) + one
      ;;iring = LONG( SQRT( ip/2.d0 - SQRT(ip/2) ) ) + 1L ; counted from NORTH pole
      iring = (cheap_isqrt(2*ip) + 1)/2L ; counted from NORTH pole, starting at 1
      iphi  = ip - 2L*iring*(iring-one)

      phi   = (iphi - 0.5d0) * !DPI/(2.d0*iring)
      z = 1.d0 - double(iring)^2 * ifact2 
      ;sz = SQRT(1.d0 - z*z) ; gets bad close to pole
      sz = iring * sqrt( ifact2 * (1.d0+z) ) ; more accurate
      vec_out[pix_np,2] = z
      vec_out[pix_np,1] = sz * SIN( phi )
      vec_out[pix_np,0] = sz * COS( phi )

      if (do_vertex) then begin
          hdelta_phi = !DPI/(4.0d0*iring) ; half pixel width
          z_nv = 1.0d0 - (iring-1.d0)^2 * ifact2 
          z_sv = 1.0d0 - (iring+1.d0)^2 * ifact2 
          iphi_mod = (iphi-1) MOD iring ; in {0,1,... iring-1}
          iphi_rat = (iphi-1)  /  iring ; in {0,1,2,3}
          phi_nv = replicate(0.0d0, n_np)
          pgt1 = where(iring gt 1, npgt1)
          if (npgt1 gt 0) then begin
              phi_nv[pgt1] = HALFPI * (iphi_rat[pgt1] +  iphi_mod[pgt1]   /double(iring[pgt1]-1))
          endif
          phi_sv           = HALFPI * (iphi_rat + (iphi_mod+1)/double(iring+1))

          vertex[pix_np,0:2,0:3] = sub_compute_vertices (z, z_nv, z_sv, phi, phi_nv, phi_sv, hdelta_phi)
          z_nv = 0 & z_sv = 0 & phi_nv = 0 & phi_sv = 0 & hdelta_phi = 0
          iphi_mod = 0 & iphi_rat = 0 & pgt1 = 0
      endif

      ip = 0 & iphi = 0 & iring = 0 & phi = 0 & z = 0 & sz = 0 ; free memory
      pix_np = 0                ; free memory

  ENDIF                         ; ------------------------------------------------------------------------

  pix_eq = WHERE(ipix GE ncap AND ipix LT nsup,  n_eq) ; equatorial strip
  IF (n_eq GT 0) THEN BEGIN     ; equatorial strip ; ---------------------------------

      ip    = ROUND(ipix[pix_eq], L64=l64) - ncap
      iring = LONG( ip / nl4) + nl1 ; counted from NORTH pole
      iphi  = ( ip MOD nl4 )  + one

      fodd  = 0.5d0 * (1 + ((iring+nl1) MOD 2)) ; 1 if iring is odd, 1/2 otherwise

      phi   = (iphi - fodd) * !DPI/(2.d0*nl1)
      z = (nl2 - iring) / fact1 
      sz = SQRT(1.d0 - z*z)
      vec_out[pix_eq,2] = z
      vec_out[pix_eq,1] = sz * SIN( phi )
      vec_out[pix_eq,0] = sz * COS( phi )

      if (do_vertex) then begin
          hdelta_phi = !DPI/(4.0d0*nl1) ; half pixel width
          phi_nv = phi
          phi_sv = phi
          z_nv = (nl2 - iring +1) / fact1 
          z_sv = (nl2 - iring -1) / fact1
          k1 = where(iring eq   nl1, nk1) ; northern transition
          k3 = where(iring eq 3*nl1, nk3) ; southern transition
          if (nk1 gt 0) then begin 
              z_nv[k1] = 1.0d0 - (nl1-1.d0)^2 * ifact2 
              iphi_mod = (iphi[k1]-1) MOD nl1 ; in {0,1,... nside-1}
              iphi_rat = (iphi[k1]-1)  /  nl1 ; in {0,1,2,3}
              if (nl1 gt 1) then phi_nv[k1] = HALFPI * (iphi_rat +  iphi_mod   /double(nl1-1))
          endif
          if (nk3 gt 0) then begin
              z_sv[k3] = -1.0d0 + (nl1-1.d0)^2 * ifact2 
              iphi_mod = (iphi[k3]-1) MOD nl1 ; in {0,1,... iring-1}
              iphi_rat = (iphi[k3]-1)   / nl1 ; in {0,1,2,3}
              if (nl1 gt 1) then phi_sv[k3] = HALFPI * (iphi_rat +  iphi_mod   /double(nl1-1))
          endif
          vertex[pix_eq,0:2,0:3] = sub_compute_vertices (z, z_nv, z_sv, phi, phi_nv, phi_sv, hdelta_phi)
          z_nv = 0 & z_sv = 0 & phi_nv = 0 & phi_sv = 0 & hdelta_phi = 0
          iphi_mod = 0 & iphi_rat = 0 & k1 = 0 & k3 = 0
      endif

      ip = 0 & iphi = 0 & iring = 0 & phi = 0 & z = 0 & sz = 0 & fodd = 0 ; free memory
      pix_eq = 0                ; free memory

  ENDIF                         ; ------------------------------------------------------------------------

  pix_sp = WHERE(ipix GE nsup,   n_sp) ; south polar cap
  if ((n_np + n_eq + n_sp) NE np ) then message,'error in '+routine
  IF (n_sp GT 0) THEN BEGIN     ; south polar cap ; ---------------------------------

      ip =  npix - ROUND(ipix[pix_sp], L64=l64)
;       iring = LONG( SQRT( ip/2.d0 - SQRT(ip/2) ) ) + 1L ; counted from SOUTH pole
;       iphi  = one + four*iring  - (ip - 2L*iring*(iring-one))
      iring = (cheap_isqrt(2*ip) + 1)/2L                           ; counted from SOUTH pole, starting at 1
      iphi  = one  - ip + 2L*iring*(iring+one) ; in [1, 4*iring]
      
      phi   = (iphi - 0.5d0) * !DPI/(2.d0*iring)
      z = - 1.d0 + double(iring)^2 * ifact2 
      ;sz = SQRT(1.d0 - z*z) ; gets bad close to pole
      sz = iring * sqrt( ifact2 * (1.d0-z) ) ; more accurate
      vec_out[pix_sp,2] = z
      vec_out[pix_sp,1] = sz * SIN( phi )
      vec_out[pix_sp,0] = sz * COS( phi )

      if (do_vertex) then begin
          hdelta_phi = !DPI/(4.0d0*iring) ; half pixel width
          z_nv = -1.0d0 + (iring+1.d0)^2 * ifact2 
          z_sv = -1.0d0 + (iring-1.d0)^2 * ifact2 
          iphi_mod = (iphi-1) MOD iring ; in {0,1,... iring-1}
          iphi_rat = (iphi-1)  /  iring ; in {0,1,2,3}
          phi_sv = replicate(0.0d0, n_sp)
          pgt1 = where(iring gt 1, npgt1)
          if (npgt1 gt 0) then begin
              phi_sv[pgt1] = HALFPI * (iphi_rat[pgt1] +  iphi_mod[pgt1]   /double(iring[pgt1]-1))
          endif
          phi_nv           = HALFPI * (iphi_rat + (iphi_mod+1)/double(iring+1))
          vertex[pix_sp,0:2,0:3] = sub_compute_vertices (z, z_nv, z_sv, phi, phi_nv, phi_sv, hdelta_phi)
          z_nv = 0 & z_sv = 0 & phi_nv = 0 & phi_sv = 0 & hdelta_phi = 0
          iphi_mod = 0 & iphi_rat = 0 & pgt1 = 0
      endif

      ip = 0 & iphi = 0 & iring = 0 & phi = 0 & z = 0 & sz = 0 ; free memory
      pix_sp = 0                ; free memory

  ENDIF                         ; ------------------------------------------------------------------------


  RETURN
END
