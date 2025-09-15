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
PRO pix2vec_nest, nside, ipix, vec_out, vertex
;************************************************************************************
;+
; PIX2VEC_NEST, Nside, Ipix, Vec_out [, Vertex]
;
;       renders cartesian coordinates Vec_out of the nominal pixel center
;       given the NESTED scheme pixel number Ipix and map resolution parameter Nside
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;       should be a power of 2
;	SCALAR
;    Ipix  : pixel number in the NEST scheme of Healpix pixelisation in [0,Npix-1]
;	can be an ARRAY of size (np) 
;
; OUTPUT
;    Vec_out : (x,y,z) position unit vector(s) with North pole = (0,0,1)
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
;    March 1999,         EH
;         correction of an error when all pixels are in the same
;         regime (eg. equatorial)
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    Nov  2002,          EH
;           added vertex output
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Oct 2006, EH, IAP, enabled nside > 8192
;    Aug 2008, EH, IAP: issues warning if ipix is not of integer type
;    Aug 2011, EH, IAP: more precise (x,y) determination close to pole 
;       (for center only, not vertex)
;
;-
;*****************************************************************************
  routine = 'PIX2VEC_NEST'
  if N_params() lt 3 then begin
      message,' syntax: pix2vec_nest, nside, ipix, vec [,vertex]'
  endif

  if (N_ELEMENTS(nside) GT 1) then message,'Nside should be a scalar in '+routine
  npix = nside2npix(nside, error = error)
  if (error ne 0) then message,'Invalid Nside '+string(nside)

  assert_pixindex_type, ipix[0], /warning ; warning if ipix is not integer
  min_pix = MIN(ipix, max=max_pix)
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

;   coordinate of the lowest corner of each face
  jrll = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4] ; in unit of nside
  jpll = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7] ; in unit of nside/2

  nl1 = LONG(nside[0]) ; 2018-06-06
  nl2 = 2L*nl1
  nl3 = 3L*nl1
  if (nl1 gt 8192) then begin
      nl4 = 4LL*nl1
      npface = nl1 * LONG64(nl1)
      l64 = 1
  endif else begin
      nl4 = 4L*nl1
      npface = nl1 * nl1
      l64 = 0
  endelse
  np = N_ELEMENTS(ipix)
  fn = DOUBLE(nl1)
  fact1 = 1.d0/(3.d0*fn*fn)
  fact2 = 2.d0/(3.d0*fn)
  piover2 = !DPI / 2.d0
  nr     = LONARR(np, /NOZERO)
  kshift = BYTARR(np, /NOZERO)
  do_vertex = (N_params() eq 4)

;     initiates the array for the pixel number -> (x,y) mapping
  common pix2xy, pix2x, pix2y
  if (n_elements(pix2x) eq 0) then init_pix2xy ; initiate pix2x and pix2y

;     finds the face, and the number in the face
  ipf = ROUND(ipix, L64=l64)  ; make sure pixel number is integer
  face_num = BYTE(ipf/npface)  ; face number in {0,11}
  ipf MOD= npface         ; pixel number in the face {0,npface-1}

;     finds the x,y on the face (starting from the lowest corner)
;     from the pixel number
  scale = 1L
  ix = 0L & iy = 0L
  smax = (nl1 le 8192) ? 2 : 5
  for i=0, smax-1 do begin
      ip_low = ipf and 1023 ; last 10 bits
      ix += scale * pix2x[ip_low]
      iy += scale * pix2y[ip_low]
      scale *= 32L
      ipf /= 1024               ; truncate out last 10 bits
  endfor
  ix += scale * pix2x[ipf]
  iy += scale * pix2y[ipf]
  ip_low = 0 & ipf = 0 ; free memory
  
;     transforms this in (horizontal, vertical) coordinates
  jrt = ix + iy                 ; 'vertical' in {0,2*(nside-1)}
  jpt = ix - iy                 ; 'horizontal' in {-nside+1,nside-1}
  ix = 0 & iy = 0               ; free memory

;     computes the z coordinate on the sphere
  jr =  jrll(face_num)*nl1 - jrt - 1 ; ring number in {1,4*nside-1}
  jrt = 0                       ; free memory

  pix_eqt = WHERE( jr GE nl1 AND jr Le nl3, n_eqt) ; equatorial region
  IF (n_eqt GT 0) THEN BEGIN   ; equatorial region 
      nr(pix_eqt)     = nl1
      kshift(pix_eqt) = (jr(pix_eqt) - nl1) MOD 2
  ENDIF

  pix_npl = WHERE( jr LT nl1, n_npl) ; north pole
  IF (n_npl GT 0) THEN BEGIN   ; north pole region
      nr(pix_npl)     = jr(pix_npl)
      kshift[pix_npl] = 0B
  ENDIF

  pix_spl = WHERE( jr GT nl3,   n_spl) ; south pole
  IF (n_npl + n_spl + n_eqt) NE np THEN message,'error in '+routine
  IF (n_spl GT 0) THEN BEGIN   ; south pole region
      nr(pix_spl)     = nl4 - jr(pix_spl)
      kshift[pix_spl] = 0B
  ENDIF

;     computes the phi coordinate on the sphere, in [0,2Pi]
  jp = (jpll(face_num)*LONG(nr) + jpt + 1 + kshift)/2 ; 'phi' number in the ring in {1,4*nr}
  jpt = 0  & face_num = 0         ; free memory
  jp = jp - nl4 * (jp GT nl4)
  jp = jp + nl4 * (jp LT 1)

  if (do_vertex) then begin
      iphi_mod = (jp-1) MOD nr ; in {0,1,... nr-1}
      iphi_rat = (jp-1) / nr      ; in {0,1,2,3}
      phi_up = piover2 * (iphi_rat +  iphi_mod   /double((nr-1)>1))
      phi_dn = piover2 * (iphi_rat + (iphi_mod+1)/double(nr+1))
      iphi_rat = 0 & iphi_mod = 0
      vertex = dblarr(np, 3, 4, /NOZERO)
   endif
  vec_out = DBLARR(np, 3, /NOZERO)

  if (n_eqt GT 0) then begin
      phi = (jp[pix_eqt] - (kshift[pix_eqt]+1)*0.5d0) * (piover2 / nr[pix_eqt])
      z           = (nl2-jr(pix_eqt))*fact2
      sz          = SQRT(1.d0 - z^2)
      vec_out(pix_eqt,0) = sz*cos(phi) & vec_out(pix_eqt,1) = sz*sin(phi) & vec_out(pix_eqt,2) = z
      if (do_vertex) then begin
          z_nv = (nl2-jr[pix_eqt]+1)*fact2
          z_sv = (nl2-jr[pix_eqt]-1)*fact2
          k1 = where(jr[pix_eqt] eq   nl1, nk1) ; northern transition
          k3 = where(jr[pix_eqt] eq 3*nl1, nk3) ; southern transition
          hdelta_phi = !dpi/(4.0d0*nr[pix_eqt])
          vertex[pix_eqt,*,*] = sub_compute_vertices(z, z_nv, z_sv, phi, phi, phi, hdelta_phi)
          if (nk1 gt 0) then begin
              z_nv[k1] =  1.0d0 - (nl1-1.d0)^2 * fact1
              ; phi_nv = phi_up
              vertex[pix_eqt[k1],*,*] = sub_compute_vertices(z[k1],z_nv[k1],z_sv[k1],phi[k1],phi_up[pix_eqt[k1]],phi[k1], hdelta_phi[k1])
          endif
          if (nk3 gt 0) then begin
              z_sv[k3] = -1.0d0 + (nl1-1.d0)^2 * fact1
              ; phi_sv = phi_up
              vertex[pix_eqt[k3],*,*] = sub_compute_vertices(z[k3],z_nv[k3],z_sv[k3],phi[k3],phi[k3],phi_up[pix_eqt[k3]], hdelta_phi[k3])
          endif
          z_nv = 0 & z_sv = 0
          hdelta_phi = 0 & k1 = 0 & k3 = 0
      endif
      pix_eqt = 0 & sz = 0 & z = 0 & phi = 0 ; free  memory
  endif

  if (n_npl GT 0) then begin
      phi = (jp[pix_npl] - (kshift[pix_npl]+1)*0.5d0) * (piover2 / nr[pix_npl])
      z           = 1.d0 - double(nr[pix_npl])^2 * fact1
      ;sz          = SQRT(1.d0 - z^2)
      sz          = nr[pix_npl] * sqrt( fact1 * (1.d0+z) ) 
      vec_out(pix_npl,0) = sz*cos(phi) & vec_out(pix_npl,1) = sz*sin(phi) & vec_out(pix_npl,2) = z
      if (do_vertex) then begin
          z_nv = 1.0d0 - (nr[pix_npl]-1.d0)^2*fact1
          z_sv = 1.0d0 - (nr[pix_npl]+1.d0)^2*fact1
          vertex[pix_npl,*,*] = sub_compute_vertices(z, z_nv, z_sv, phi, phi_up[pix_npl], phi_dn[pix_npl], !dpi/(4.0d0*nr[pix_npl]))
          z_nv = 0 & z_sv = 0
      endif
      pix_npl = 0 & sz = 0 & z = 0 & phi = 0 ; free  memory
  endif

  if (n_spl GT 0) then begin
      phi = (jp[pix_spl] - (kshift[pix_spl]+1)*0.5d0) * (piover2 / nr[pix_spl])
      z           = -1.d0 + double(nr[pix_spl])^2 * fact1
      ; sz          = SQRT(1.d0 - z^2)
      sz          = nr[pix_spl] * sqrt( fact1 * (1.d0-z) ) 
      vec_out(pix_spl,0) = sz*cos(phi) & vec_out(pix_spl,1) = sz*sin(phi) & vec_out(pix_spl,2) = z
      if (do_vertex) then begin
          z_nv = - 1.0d0 + (nr[pix_spl]+1.d0)^2*fact1
          z_sv = - 1.0d0 + (nr[pix_spl]-1.d0)^2*fact1
          vertex[pix_spl,*,*] = sub_compute_vertices(z, z_nv, z_sv, phi, phi_dn[pix_spl], phi_up[pix_spl], !dpi/(4.0d0*nr[pix_spl]))
          z_nv = 0 & z_sv = 0
      endif
      pix_spl = 0 & sz = 0 & z = 0 & phi = 0 ; free  memory
  endif
  
  jp = 0 & kshift = 0 & nr = 0
  jr = 0

  return
end ; pix2vec_nest

