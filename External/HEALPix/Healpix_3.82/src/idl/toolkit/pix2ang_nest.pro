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
PRO pix2ang_nest, nside, ipix, theta, phi
;*****************************************************************************
;+
; PIX2ANG_NEST, nside, ipix, theta, phi
;
;       renders Theta and Phi coordinates of the nominal pixel center
;       given the NESTED scheme pixel number Ipix and map resolution parameter Nside
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2),
;       should be a power of 2
;	SCALAR
;    Ipix  : pixel number in the NESTED scheme of Healpix pixelisation in [0,Npix-1]
;	can be an ARRAY 
;
; OUTPUT
;    Theta : angle (along meridian = co-latitude), in [0,Pi], theta=0 : north pole,
;	is an ARRAY of same size as Ipix
;    Phi   : angle (along parallel = azimut), in [0,2*Pi]
;	is an ARRAY of same size as Ipix
;
; SUBROUTINE
;    nside2npix
;
; HISTORY
;     Feb 1999,           Eric Hivon,               Caltech
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Oct 2006, EH, IAP, enabled nside > 8192
;    Aug 2008, EH, IAP: issues warning if ipix is not of integer type
;    Aug 2011, EH, IAP: more precise theta determination close to pole;
;        sligtly faster by not using kshift
;    Sep 2011, EH, IAP: use HISTOGRAM instead of multiple WHEREs
;    Apr 2013, EH, IAP: works with scalar ipix
;
;-
;*****************************************************************************
  routine = 'PIX2ANG_NEST'
  if N_params() lt 3 then begin
      message,' syntax: '+routine+', nside, ipix, theta, phi'
  endif

  if (N_ELEMENTS(nside) GT 1) then message,'Nside should be a scalar in '+routine
  npix = nside2npix(nside, error = error)
  if (error ne 0) then message,'Invalid Nside '+string(nside)

  assert_pixindex_type, ipix[0], /warning ; warning if ipix is not integer
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

;   coordinate of the lowest corner of each face
  jrll = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4] ; in unit of nside
  jpll = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7] ; in unit of nside/2


  nl1 = LONG(nside[0]) ; 2018-06-06
  nl3 = 3L*nl1
  if (nl1 gt 8192) then begin
      nl4 = 4LL*nl1
      npface = nl1 * long64(nl1)
      l64 = 1
      if (nl1 ge 2LL^28) then jpll=long64(jpll) ; necessary to accomodate 8*Nside
  endif else begin
      nl4 = 4L*nl1
      npface = nl1 * nl1
      l64 = 0
  endelse
  np = N_ELEMENTS(ipix)
  fn = DOUBLE(nl1)
;  fact1 = 1.d0/(3.d0*fn*fn)
  sfact = 1.d0 / (sqrt(6.d0) * fn)
  fact2 = 2.d0/(3.d0*fn)
  piover4 = !DPI / 4.d0

;     initiates the array for the pixel number -> (x,y) mapping
  common pix2xy, pix2x, pix2y
  if (n_elements(pix2x) eq 0) then init_pix2xy ; initiate pix2x and pix2y

;     finds the face, and the number in the face
  ipf = ROUND(ipix, L64=l64)  ; make sure pixel number is integer
  face_num = BYTE(ipf/npface)        ; face number in {0,11}
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
  jr =  long(jrll(face_num)*nl1 - jrt - 1) ; ring number in {1,4*nside-1}
  jrt = 0                       ; free memory

  nr     = LONARR(np, /NOZERO)
  theta  = DBLARR(np, /NOZERO)

  if (n_elements(jr) eq 1) then jr=[jr] ; make sure jr is an array (for HISTOGRAM)
  hist = HISTOGRAM(jr, min=-nl1, max=nl4, binsize=2L*nl1, reverse=rev)
  if (hist[0]+hist[1]+hist[2]) NE np THEN message,'error in '+routine

  kc = 1 ; equatorial region: [nl1,nl3]
  IF (hist[kc] GT 0) THEN BEGIN
      pix_eqt = rev[rev[kc]:rev[kc+1]-1]
      nr[pix_eqt]     = nl1   ; equatorial region 
      theta[pix_eqt]  = ACOS(  (2*nl1-jr[pix_eqt])*fact2  )
      pix_eqt = 0  ; free memory
  ENDIF

  kc = 0 ; north pole: [1,nl1]
  IF (hist[kc] GT 0) THEN BEGIN
      pix_npl = rev[rev[kc]:rev[kc+1]-1]
      nr[pix_npl]     = jr[pix_npl]
      theta[pix_npl] = 2.d0 * ASIN( nr[pix_npl] * sfact )
      pix_npl = 0; free memory
  ENDIF

  kc = 2 ; south pole: [nl3,nl4-1]
  IF (hist[kc] GT 0) THEN BEGIN
      pix_spl = rev[rev[kc]:rev[kc+1]-1]
      nr[pix_spl]     = nl4 - jr[pix_spl]
      theta[pix_spl] = !DPI - 2.d0 * ASIN( nr[pix_spl] * sfact )
      pix_spl = 0; free memory
  ENDIF
  rev = 0 & hist=0

; ;     computes the phi coordinate on the sphere, in [0,2Pi]
  jp = jpll[face_num]*nr + jpt
  jpt = 0  & face_num = 0         ; free memory
  jp += (2*nl4) * (jp LT 0)    ; in [0,8Nside-1] because of 1/2 step stagger
  phi = jp * (piover4 / nr)

  jp = 0
  nr = 0
  
  return
end ; pix2ang_nest

