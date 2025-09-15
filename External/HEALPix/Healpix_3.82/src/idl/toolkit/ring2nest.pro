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
PRO ring2nest, nside, ipring, ipnest
;*******************************************************************************
;+
;   RING2NEST, Nside, Ipring, Ipnest
;
;         performs conversion from RING to NESTED pixel number in
;         the Healpix tesselation for a parameter Nside
;
;   INPUT :
;     Nside  : determines the resolution (Npix = 12* Nside^2)
;         should be a power of 2 (not tested)
;	  SCALAR
;     Ipring : pixel number in the RING scheme of Healpix pixelisation in [0,Npix-1]
;	can be an ARRAY
;
;   OUTPUT :
;     Ipnest : pixel number in the NEST scheme of Healpix pixelisation in [0,Npix-1]
;	is an ARRAY of same size as Ipring
;
;
; HISTORY
;    Feb 1999,           Eric Hivon,               Caltech
;    April 1999, EH : correction of a bug on the face number in the
;    North pole regime, (produced only if the input pixel indices were not consecutive)
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Dec 2007, EH, enabled Nside > 8192
;    Aug 2008, EH, IAP: issues warning if ipix is not of integer type
;    Oct 2011, EH, IAP: uses CHEAP_ISQRT to compute ring index out of pixel index
;
;-
;*******************************************************************************
  routine = 'RING2NEST'

;   ns_max = 8192L

  if N_params() ne 3 then begin
      message,' syntax: '+routine+', nside, ipring, ipnest'
  endif

  if (N_ELEMENTS(nside) GT 1) then begin
      message,'nside should be a scalar in '+routine
  endif
  npix = nside2npix(nside, error = error)
  if (error ne 0) then message, 'Invalid Nside: '+string(nside)

  assert_pixindex_type, ipring[0], /warning ; warning if ipix is not integer
  min_pix = MIN(ipring, MAX = max_pix)
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

  np = N_ELEMENTS(ipring)
  nl1 = LONG(nside[0]) ; 2018-06-06
  nl2 = 2L*nl1
  if (nl1 gt 8192) then begin
      nl4 = 4LL*nl1
      nl8 = 8LL*nl1
      ncap = nl2*(nl1-1LL)
      nsup = nl2*(5LL*nl1+1LL)
      npface = nl1 * LONG64(nl1)
      iphi   = LON64ARR(np, /NoZero)
      iring  = LONARR(np, /NoZero)
      nr     = LONARR(np, /NoZero)
      one    = 1LL
      l64    = 1
  endif else begin
      nl4 = 4L*nl1
      nl8 = 8L*nl1
      ncap = nl2*(nl1-1)
      nsup = nl2*(5*nl1+1)
      npface = nl1 * (nl1)
      iphi   = LONARR(np, /NoZero)
      iring  = INTARR(np, /NoZero)
      nr     = INTARR(np, /NoZero)
      one    = 1L
      l64    = 0
  endelse
  face_num = BYTARR(np, /NoZero)
  kshift   = BYTARR(np, /NoZero)

  common xy2pix, x2pix, y2pix
  if (n_elements(x2pix) eq 0) then init_xy2pix ; initiate x2pix and y2pix


  pix_npl = WHERE(ipring LT ncap,   n_npl) ; north polar cap
  IF (n_npl GT 0) THEN BEGIN     ; north polar cap ; ---------------------------------

      ip = ROUND(ipring[pix_npl], L64=l64) + one
;;;      irn = LONG( SQRT( ip/2.d0 - SQRT(ip/2) ) ) + 1L ; counted from NORTH pole
      irn = (cheap_isqrt(2*ip) + 1)/2L ; counted from NORTH pole

      iring[pix_npl]  = irn
      iphi[pix_npl]   = ip - 2L*irn*(irn-one)
      kshift[pix_npl] = 0
      nr[pix_npl]     = irn ; 1/4 of the number of points on the current ring
      face_num[pix_npl] = (iphi[pix_npl]-1) / irn ; in {0,3}
      ip = 0 & irn = 0
      pix_npl = 0

  ENDIF                         ; ------------------------------------------------

  pix_eqt = WHERE(ipring GE ncap AND ipring LT nsup,  n_eqt) ; equatorial strip
  IF (n_eqt GT 0) THEN BEGIN     ; equatorial strip ; ---------------------------------

      ip    = ROUND(ipring[pix_eqt], L64=l64) - ncap

      iring[pix_eqt] = LONG( ip / nl4) + nl1 ; counted from NORTH pole
      iphi[pix_eqt]  = ( ip MOD nl4 )  + 1
      kshift[pix_eqt]  = (iring[pix_eqt]+nl1) MOD 2 ; 1 if iring+nside is odd, 0 otherwise
      nr[pix_eqt] = nl1
      ip =0

      ire =  iring[pix_eqt] - nl1 + 1 ; in {1, 2*nside +1}
      irm =  nl2 + 2 - ire
      ifm = (iphi[pix_eqt] - ire/2 + nl1 -1) / nl1 ; face boundary
      ifp = (iphi[pix_eqt] - irm/2 + nl1 -1) / nl1
      ire = 0 & irm = 0         ; free memory
      ifd = ifp - ifm
      ifm = 0                   ; free memory
      sub1 = ( (ifp MOD 4) + 4) * (ifd EQ 0) ; faces 4 to 7
      sub2 = (  ifp           ) * (ifd EQ -1) ; (half-)faces 0 to 3
      sub3 = (  ifp + 7       ) * (ifd EQ 1) ; (half-)faces 8 to 11

      face_num[pix_eqt] = sub1 + sub2 + sub3
      sub1 = 0 & sub2 = 0 & sub3 = 0 & ifd = 0 & ifp = 0
      pix_eqt = 0

  ENDIF                         ; -------------------------------------------------

  pix_spl = WHERE(ipring GE nsup,   n_spl) ; south polar cap
  IF (n_npl + n_spl + n_eqt) NE np THEN message, 'error in ring2nest'
  IF (n_spl GT 0) THEN BEGIN     ; south polar cap ; ---------------------------------
      
      ip =  npix - ROUND(ipring[pix_spl], L64=l64)
;;      irs = LONG( SQRT( ip/2.d0 - SQRT(ip/2) ) ) + 1 ; counted from SOUTH pole
      irs = (cheap_isqrt(2*ip) + 1)/2L ; counted from SOUTH pole

      iring[pix_spl]   = nl4 - irs
;;      iphi[pix_spl]  = 4*irs + one - (ip - 2L*irs*(irs-one))
      iphi[pix_spl]  = one - ip + 2L*irs*(irs+one)
      kshift[pix_spl] = 0
      nr[pix_spl] = irs
      face_num[pix_spl] = (iphi[pix_spl]-1) / irs + 8 ; in {8,11}
      ip = 0 & irs = 0
      pix_spl = 0

  ENDIF                         ; -------------------------------------------------

;     finds the (x,y) on the face
  irt = LONG(  iring  - jrll(face_num)*nl1 + 1 )      ; in {-nside+1,0}
  ipt = LONG( 2*iphi  - jpll(face_num)*nr - kshift - 1) ; in {-nside+1,nside-1}
  iring = 0 & kshift = 0 & nr = 0 ; free memory
  ipt = ipt - nl8 * (ipt GE nl2) ; for the face #4
  ix =  (ipt - irt ) / 2
  iy = -(ipt + irt ) / 2
  ipt = 0 & irt = 0 & iphi = 0  ; free memory

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
END   ; ring2nest

