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
PRO nest2ring, nside, ipnest, ipring
;*******************************************************************************
;+
;   NEST2RING, Nside, Ipnest, Ipring
;
;         performs conversion from NESTED to RING pixel number in
;         the Healpix tesselation for a parameter Nside
;
;   INPUT :
;     Nside  : determines the resolution (Npix = 12* Nside^2)
;         should be a power of 2 (not tested)
;	  SCALAR
;     Ipnest : pixel number in the NEST scheme of Healpix pixelisation in [0,Npix-1]
;	can be an ARRAY
;
;   OUTPUT :
;     Ipring : pixel number in the RING scheme of Healpix pixelisation in [0,Npix-1]
;	is an ARRAY of same size as Ipnest
;
;
; HISTORY
;    Feb 1999,           Eric Hivon,               Caltech
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Dec 2007, EH, enabled Nside > 8192
;    Aug 2008, EH, IAP: issues warning if ipix is not of integer type
;
;-
;*******************************************************************************
  routine = 'NEST2RING'

  if N_params() ne 3 then begin
      message,' syntax: '+routine+', nside, ipnest, ipring'
  endif

  if (N_ELEMENTS(nside) GT 1) then begin
      message,'nside should be a scalar in '+routine
  endif
  npix = nside2npix(nside, error = error)
  if (error ne 0) then message, 'Invalid Nside: '+string(nside)

  assert_pixindex_type, ipnest[0], /warning ; warning if ipix is not integer
  min_pix = MIN(ipnest, MAX = max_pix)
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

  np = N_ELEMENTS(ipnest)
  nl1 = LONG(nside[0]) ; 2018-06-06
  nl3  = 3L*nl1
  if (nl1 gt 8192) then begin
      nl4  = 4LL*nl1
      npface = nl1*LONG64(nl1)
      ncap = 2LL*nl1*LONG64(nl1-1)
      n_before = LON64ARR(np, /NOZERO)
      nr       = LONARR(np, /NoZero)
      one      = 1LL
      l64 = 1
  endif else begin
      nl4  = 4L*nl1
      npface = nl1*nl1
      ncap = 2L*nl1*(nl1-1)
      n_before = LONARR(np, /NOZERO)
      nr       = INTARR(np, /NoZero)
      one      = 1L
      l64 = 0
  endelse
  kshift   = BYTARR(np, /NOZERO)


;     initiates the array for the pixel number -> (x,y) mapping
  common pix2xy, pix2x, pix2y
  if (n_elements(pix2x) eq 0) then init_pix2xy ; initiate pix2x and pix2y

;     finds the face, and the number in the face
  ipf = ROUND(ipnest, L64=l64)  ; make sure pixel number is integer
  face_num = BYTE(ipf/npface)      ; face number in {0,11}
  ipf MOD= npface       ; pixel number in the face {0,npface-1}

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
  jr =  jrll(face_num)*long(nl1) - jrt - 1 ; ring number in {1,4*nside-1}
  jrt = 0                       ; free memory

  pix_eqt = WHERE( jr GE nl1 AND jr Le nl3, n_eqt) ; equatorial region
  IF (n_eqt GT 0) THEN BEGIN   ; equatorial region 
      nr[pix_eqt]       = nl1
      n_before(pix_eqt) = ncap + nl4 * (jr(pix_eqt) - nl1)
      kshift(pix_eqt)   = (jr(pix_eqt) - nl1) MOD  2
      pix_eqt = 0 ; free memory
  ENDIF

  pix_npl = WHERE( jr LT nl1, n_npl) ; north pole
  IF (n_npl GT 0) THEN BEGIN    ; north pole region
      nr[pix_npl] = jr[pix_npl]
      n_before(pix_npl) = 2L * nr[pix_npl] * (nr[pix_npl] - one)
      kshift(pix_npl) = 0
      pix_npl = 0 ; free memory
  ENDIF

  pix_spl = WHERE( jr GT nl3,   n_spl) ; south pole
  IF (n_npl + n_spl + n_eqt) NE np THEN message,'error in nest2ring'
  IF (n_spl GT 0) THEN BEGIN    ; south pole region
      nr[pix_spl] = nl4 - jr(pix_spl)
      n_before(pix_spl) = npix - 2L * (nr[pix_spl] + one) * nr[pix_spl]
      kshift(pix_spl) = 0
      pix_spl = 0 ; free memory
  ENDIF

  jr = 0 ; free memory
;     computes the phi coordinate on the sphere, in [0,2Pi]
  if (nl1 gt 8192) then begin
      jp = (jpll(face_num)*LONG64(nr) + jpt + 1 + kshift)/2 ; 'phi' number in the ring in {1,4*nr}
  endif else begin
      jp = (jpll(face_num)*nr + jpt + 1 + kshift)/2 ; 'phi' number in the ring in {1,4*nr}
  endelse
  jpt = 0 ; free memory

  jp = jp - nl4*(jp GT nl4)
  jp = jp + nl4*(jp LT 1)

  ipring = n_before + jp - 1L ; in {0, npix-1}

  RETURN
END      ; nest2ring
