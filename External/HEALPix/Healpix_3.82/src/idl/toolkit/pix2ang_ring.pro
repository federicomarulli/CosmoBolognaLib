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
PRO pix2ang_ring, nside, ipix, theta, phi
;****************************************************************************************
;+
; PIX2ANG_RING, nside, ipix, theta, phi
; 
;       renders Theta and Phi coordinates of the nominal pixel center
;       given the RING scheme pixel number Ipix and map resolution parameter Nside
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;	SCALAR
;    Ipix  : pixel number in the RING scheme of Healpix pixelisation in [0,Npix-1]
;	can be an ARRAY 
;       pixels are numbered along parallels (ascending phi), 
;       and parallels are numbered from north pole to south pole (ascending theta)
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
;    June-October 1997,  Eric Hivon & Kris Gorski, TAC
;    Aug  1997 : treats correctly the case nside = 1
;    Feb 1999,           Eric Hivon,               Caltech
;         renamed pix2ang_ring
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Dec 2007, EH,  IAP, enabled nside > 8192
;    Aug 2008, EH, IAP: issues warning if ipix is not of integer type
;    Aug 2011, EH, IAP: uses CHEAP_ISQRT to compute ring index out of pixel
;      index ;  
;     more precise theta determination close to pole
;
;
;-
;****************************************************************************************
  routine = 'PIX2ANG_RING'
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
  ;fact2 = (3.d0*nl1)*nl1
  sfact = 1.d0 / (nl1 * sqrt(6.d0))
  np = N_ELEMENTS(ipix)
  theta = DBLARR(np, /NoZero)
  phi   = DBLARR(np, /NoZero)


  pix_np = WHERE(ipix LT ncap,   n_np) ; north polar cap
  IF (n_np GT 0) THEN BEGIN     ; north polar cap ; ---------------------------------

      ip    = ROUND(ipix[pix_np], L64=l64) + one ; ipix+1
      iring = (cheap_isqrt(2*ip) + 1)/2L ; counted from NORTH pole, starting at 1
      iphi  = ip - 2L*iring*(iring-one)

      theta[pix_np] = 2.d0 * ASIN( iring * sfact)
      phi[pix_np]   = (iphi - 0.5d0) * !DPI/(2.d0*iring)
      ip = 0 & iring =0 & iphi = 0 ; free memory
      pix_np = 0                ; free memory

  ENDIF                         ; ------------------------------------------------------------------------

  pix_eq = WHERE(ipix GE ncap AND ipix LT nsup,  n_eq) ; equatorial strip
  IF (n_eq GT 0) THEN BEGIN     ; equatorial strip ; ---------------------------------

      ip    = ROUND(ipix[pix_eq], L64=l64) - ncap
      iring = LONG( ip / nl4) + nl1 ; counted from NORTH pole in {Nside, 3*Nside}
      iphi  = ( ip MOD nl4 )  + one

      fodd  = 0.5d0 * (1 + ((iring+nl1) MOD 2)) ; 1 if iring is odd, 1/2 otherwise

      theta[pix_eq] = ACOS( (nl2 - iring) / fact1 )
      phi[pix_eq]   = (iphi - fodd) * !DPI/(2.d0*nl1)
      ;;print, ip, iring, iphi, fodd
      ip = 0 & iring =0 & iphi = 0 & fodd = 0; free memory
      pix_eq = 0                ; free memory
      
  ENDIF                         ; ------------------------------------------------------------------------

  pix_sp = WHERE(ipix GE nsup,   n_sp) ; south polar cap
  IF (n_np + n_sp + n_eq) NE np THEN message,'error in '+routine
  IF (n_sp GT 0) THEN BEGIN     ; south polar cap ; ---------------------------------

      ip    = npix - ROUND(ipix[pix_sp], L64=l64)            ; Npix - ipix
      iring = (cheap_isqrt(2*ip) + 1)/2L                           ; counted from SOUTH pole, starting at 1
      iphi  = one  - ip + 2L*iring*(iring+one) ; in [1, 4*iring]
      
      theta[pix_sp] = !DPI - 2.d0 * ASIN( iring * sfact)
      phi[pix_sp]   = (iphi - 0.5d0) * !DPI/(2.d0*iring)
      ip = 0 & iring =0 & iphi = 0 ; free memory
      pix_sp = 0                ; free memory


  ENDIF                         ; ------------------------------------------------------------------------


  RETURN
END
