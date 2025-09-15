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
PRO ang2pix_ring, nside, theta, phi, ipring
;*********************************************************************************
;+
; ANG2PIX_RING, Nside, Theta, Phi, Ipring
;
;        renders the RING scheme pixel number Ipring for a pixel which, given the
;        map resolution parameter Nside, contains the point on the sphere
;        at angular coordinates Theta and Phi
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;	SCALAR
;    Theta : angle (along meridian), in [0,Pi], theta=0 : north pole,
;	can be an ARRAY
;    Phi   : angle (along parallel), in [0,2*Pi]
;	can be an ARRAY of same size as theta
;
; OUTPUT
;    Ipring  : pixel number in the RING scheme of HEALPIX pixelisation in [0,Npix-1]
;	can be an ARRAY of same size as Theta and Phi
;    pixels are numbered along parallels (ascending phi), 
;    and parallels are numbered from north pole to south pole (ascending theta)
;
;
; SUBROUTINE
;    nside2npix
;
; HISTORY
;    June-October 1997,  Eric Hivon & Kris Gorski, TAC, 
;            original ang_pix
;    Feb 1999,           Eric Hivon,               Caltech
;            name changed to ang2pix_ring_new
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Aug  2004,  EH, use !PI as theta upper-bound instead of !DPI
;    Dec 2007, EH,  IAP, enabled nside > 8192
;    Mar 2011, EH, IAP: MOD -> AND, replaced repeated WHERE with HISTOGRAM
;    Apr 2011, EH, IAP, correction of a bug affecting phi out of [-2Pi, 2Pi] 
;    Aug 2011, EH, IAP: more accurate calculations close to pole
;    Apr 2013, EH, IAP: works with scalar theta and phi
;-
;
; in North polar cap
; z = cos(theta) = 1 - i^2/(3.*Nside^2)  with i in {1,Nside}
; sin(theta/2)   = i/Nside /sqrt(6)
; sin(theta)     = i/Nside * sqrt( (1+z)/3 )
;
;
;*********************************************************************************
  routine = 'ANG2PIX_RING'
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
  nl1  = LONG(nside[0]) ; 2018-06-06
  nl2   = 2L*nl1
  l64 = 0
  if (nl1 gt 8192) then begin
      nl4   = 4LL*nl1
      ncap  = nl2*(nl1-1LL)
      ipring = LON64ARR(np, /NoZero)
      one = 1LL
      l64 = (nl1 ge ishft(1L,29)) ; for Nside >= 2^29
  endif else begin
      nl4   = 4L*nl1
      ncap  = nl2*(nl1-1L)
      ipring = LONARR(np, /NoZero)
      one = 1L
  endelse
  pion2 = !DPI * 0.5d0
  twopi = !DPI * 2.d0

  cth_in = COS(DOUBLE(theta))
  phi_in = phi MOD twopi
  ;;;;;phi_in = phi + (phi LE 0.d0)*twopi
  phi_in += (phi LE 0.d0)*twopi ; bug correction 2011-04-28

  if (n_elements(cth_in) eq 1) then cth_in = [cth_in] ; make sure cth_in is an array (for HISTOGRAM)
  histo = histogram( cth_in * 1.5d0, min=-3.d0, max=3.d0, binsize=2.d0, reverse=rev)
  n_sp  = histo[0]
  n_eqt = histo[1]
  n_np  = histo[2]
  if ((n_sp+n_np+n_eqt) ne np) then message,'Error in '+routine
  IF (n_eqt GT 0) THEN BEGIN    ; equatorial strip ----------------
      pix_eqt = rev[rev[1]:rev[2]-1]
      tt = 0.5d0 + phi_in[pix_eqt] / pion2
      zz = cth_in[pix_eqt]*0.75d0

      jp = floor(nl1*( tt - zz), L64=l64) ; increasing edge line index
      jm = floor(nl1*( tt + zz), L64=l64) ; decreasing edge line index

      ir = (nl1 + 1) + jp - jm ; in {1,2n+1} (ring number counted from z=2/3)
      k =  ~(ir and 1)   ; k=1 if ir even, and 0 otherwise

      ;ip = LONG( ( jp+jm+k + (1-nl1) ) / 2 ) and (nl4-1) ; in {0,4n-1}
      ip = LONG( ( (one-nl1) + jp +jm + k  ) / 2 ) and (nl4-1) ; in {0,4n-1} ; intermediate values are 64bits
      
      ipring[pix_eqt] = ncap + nl4*(ir-one) + ip
      tt = 0 & zz = 0 & jp = 0 & jm = 0 & ir = 0 & k = 0 & ip = 0
      pix_eqt = 0
  ENDIF

  IF (n_np GT 0) THEN BEGIN     ; north polar caps ------------------------

      pix_np  = rev[rev[2]:rev[3]-1]
      tt = phi_in[pix_np] / pion2
      tp = tt MOD 1.d0
      ;tmp = SQRT( 3.d0*(1.d0 - ABS(cth_in[pix_np])) )
      tmp = SQRT(6.d0) * SIN( double(theta[pix_np]) * 0.5d0) ; more accurate for theta~0

      jp = LONG( nl1 * tp          * tmp ) ; increasing edge line index
      jm = LONG( nl1 * (1.d0 - tp) * tmp ) ; decreasing edge line index

      ir = jp + jm + 1          ; ring number counted from the closest pole
      ip = LONG( tt * ir ) + 1  ; in {1,4*ir}
      ir4 = 4*ir
      ip -= ir4*(ip GT ir4)

      ipring[pix_np] =        2*ir*(ir-one) + ip - one
      tt = 0 & tp = 0 & tmp =0 & jp = 0 & jm = 0 & ir = 0 & ip = 0 & ir4 = 0
      pix_np = 0
  ENDIF                         ; -------------------------------------------------------

  IF (n_sp GT 0) THEN BEGIN     ; south polar caps ------------------------

      pix_sp  = rev[rev[0]:rev[1]-1]
      tt = phi_in[pix_sp] / pion2
      tp = tt MOD 1.d0
      ;tmp = SQRT( 3.d0*(1.d0 - ABS(cth_in[pix_sp])) )
      tmp = SQRT(6.d0) * COS( double(theta[pix_sp]) * 0.5d0) ; more accurate for theta~Pi

      jp = LONG( nl1 * tp          * tmp ) ; increasing edge line index
      jm = LONG( nl1 * (1.d0 - tp) * tmp ) ; decreasing edge line index

      ir = jp + jm + 1          ; ring number counted from the closest pole
      ip = LONG( tt * ir ) + 1  ; in {1,4*ir}
      ir4 = 4*ir
      ip -= ir4*(ip GT ir4)

      ipring[pix_sp] = npix - 2*ir*(ir+one) + ip - one
      tt = 0 & tp = 0 & tmp = 0 & jp = 0 & jm = 0 & ir = 0 & ip = 0 & ir4 = 0
      pix_sp = 0
  ENDIF                         ; -------------------------------------------------------

  return
end                             ; ang2pix_ring

