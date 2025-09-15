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
PRO vec2pix_ring, nside, vec_in, ipring
;*******************************************************************
;+
; VEC2PIX_RING, Nside, Vec_in, Ipring
;
;        renders the RING scheme pixel number Ipring for a pixel which, given the
;        map resolution parameter Nside, contains the point on the sphere
;        at cartesian coordinates Vec_in
;
; INPUT
;    Nside     : determines the resolution (Npix = 12* Nside^2)
;	     SCALAR 
;    Vec_in    : (x,y,z) position unit vector(s) with North pole = (0,0,1)
;            stored as x(0), x(1), ..., y(0), y(1), ..., z(0), z(1) ..
;            ARRAY of dimension = (np,3)
;
; OUTPUT
;    Ipring  : pixel number in the RING scheme in Healpix pixelisation in [0,Npix-1]
;	is a VECTOR of dimension = (np)
;    pixels are numbered along parallels
;    and parallels are numbered from north pole to south pole
;
;
; SUBROUTINE
;    nside2npix
;
; HISTORY
;    June-October 1997,  Eric Hivon & Kris Gorski, TAC
;    Feb 1999,  Eric Hivon, Caltech
;    Sept 2000,          EH
;           free memory by discarding unused variables
;    June 2003,  EH, replaced STOPs by MESSAGEs
;    Dec 2007, EH,  IAP, enabled nside > 8192
;    Aug 2011, EH, IAP: more accurate calculations close to pole
;
;-
;*****************************************************************************
  routine = 'VEC2PIX_RING'
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
  nl1  = LONG(nside[0]) ; 2018-06-06
  nl2   = 2*nl1
  l64 = 0
  if (nl1 gt 8192) then begin
      nl4   = 4LL*nl1
      ncap  = nl2*(nl1-1LL)
      ipring = LON64ARR(np, /NoZero)
      one = 1LL
      l64 = (nl1 ge ishft(1L, 29)) ; for Nside >= 2^29
  endif else begin
      nl4   = 4L*nl1
      ncap  = nl2*(nl1-1L)
      ipring = LONARR(np, /NoZero)
      one = 1L
  endelse
  pion2 = !DPI * 0.5d0
  twopi = !DPI * 2.d0
  cth0 = 2.d0/3.d0

  norm = SQRT( TOTAL(vec_in^2, 2, /double) )
  cth_in = vec_in[*,2]/ norm
  sth_in = SQRT( vec_in[*,0]^2 + vec_in[*,1]^2 ) / norm
  norm = 0
  phi_in = ATAN(vec_in[*,1],vec_in[*,0])
  phi_in += (phi_in LE 0.d0)*twopi
  pix_eqt = WHERE(cth_in LE cth0 AND cth_in GT -cth0, n_eqt) ; equatorial strip
  IF (n_eqt GT 0) THEN BEGIN    ; equatorial strip ----------------
      tt = phi_in(pix_eqt) / pion2

;       jp = LONG(nl1*(0.5d0 + tt - cth_in(pix_eqt)*0.75d0)) ; increasing edge line index
;       jm = LONG(nl1*(0.5d0 + tt + cth_in(pix_eqt)*0.75d0)) ; decreasing edge line index
      jp = floor(nl1*(0.5d0 + tt - cth_in(pix_eqt)*0.75d0), L64=l64) ; increasing edge line index
      jm = floor(nl1*(0.5d0 + tt + cth_in(pix_eqt)*0.75d0), L64=l64) ; decreasing edge line index

      ir = (nl1 + 1) + jp - jm ; in {1,2n+1} (ring number counted from z=2/3)
      k =  ( (ir MOD 2) EQ 0)   ; k=1 if ir even, and 0 otherwise

      ;ip = LONG( ( jp+jm+k + (1-nl1) ) / 2 ) + 1 ; in {1,4n}
      ip = floor( ( one-nl1 + jp+jm+k  ) / 2 , L64=l64) + 1 ; in {1,4n}
      ip = ip - nl4*(ip GT nl4)

      ipring(pix_eqt) = ncap + nl4*(ir-one) + ip - one
      tt = 0 & jp = 0 & jm = 0 & ir = 0 & k = 0 & ip = 0
      pix_eqt =0
  ENDIF

  pix_np  = WHERE(cth_in GT cth0, n_np) ; north caps
  IF (n_np GT 0) THEN BEGIN     ; north polar caps ------------------------

      tt = phi_in(pix_np) / pion2
      tp = tt MOD 1.d0
;       tmp = SQRT( 3.d0*(1.d0 - ABS(cth_in(pix_np))) )
      tmp = sth_in[pix_np] * SQRT(3.d0 / (1.d0 + cth_in[pix_np])) ; more accurate

      jp = LONG( nl1 * tp          * tmp ) ; increasing edge line index
      jm = LONG( nl1 * (1.d0 - tp) * tmp ) ; decreasing edge line index

      ir = jp + jm + 1          ; ring number counted from the closest pole
      ip = LONG( tt * ir ) + 1  ; in {1,4*ir}
      ir4 = 4*ir
      ip = ip - ir4*(ip GT ir4)

      ipring(pix_np) =        2*ir*(ir-one) + ip - one
      tt = 0 & tp = 0 & tmp = 0 & jp =0 & jm = 0 & ir = 0 & ip =0 & ir4=0
      pix_np = 0
  ENDIF                         ; -------------------------------------------------------

  pix_sp  = WHERE(cth_in LE -cth0, n_sp) ; south pole
  if ((n_sp + n_np + n_eqt) ne np) then message,'Error in '+routine
  IF (n_sp GT 0) THEN BEGIN     ; south polar caps ------------------------

      tt = phi_in(pix_sp) / pion2
      tp = tt MOD 1.d0
;      tmp = SQRT( 3.d0*(1.d0 - ABS(cth_in(pix_sp))) )
      tmp = sth_in[pix_sp] * SQRT(3.d0 / (1.d0 - cth_in[pix_sp])) ; more accurate

      jp = LONG( nl1 * tp          * tmp ) ; increasing edge line index
      jm = LONG( nl1 * (1.d0 - tp) * tmp ) ; decreasing edge line index

      ir = jp + jm + 1          ; ring number counted from the closest pole
      ip = LONG( tt * ir ) + 1  ; in {1,4*ir}
      ir4 = 4*ir
      ip = ip - ir4*(ip GT ir4)

      ipring(pix_sp) = npix - 2*ir*(ir+one) + ip - one
      tt = 0 & tp = 0 & tmp = 0 & jp =0 & jm = 0 & ir = 0 & ip =0 & ir4=0
      pix_sp = 0
  ENDIF                         ; -------------------------------------------------------

  return
end

