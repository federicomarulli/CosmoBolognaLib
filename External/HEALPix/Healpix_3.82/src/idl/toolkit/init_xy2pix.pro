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
PRO init_xy2pix
;************************************************************************
;
;     sets the array giving the number of the pixel lying in (x,y)
;     x and y are in {1,128}
;     the pixel number is in {0,128**2-1}
;
;     if  i-1 = sum_p=0  b_p * 2^p
;     then ix = sum_p=0  b_p * 4^p
;          iy = 2*ix
;     ix + iy in {0, 128**2 -1}
;-
;************************************************************************

  common xy2pix, x2pix, y2pix

  x2pix = LONARR(128)
  y2pix = LONARR(128)

  for i = 0,127 do begin        ;for converting x,y into
      j  = i                    ;pixel numbers
      k  = 0
      ip = 1
      loop:
      if (j eq 0) then begin
          x2pix(i) = K
      endif else begin
          id = j mod 2
          j  = j/2
          k  = IP*ID+K
          ip = ip*4
          goto, loop
      endelse
  endfor
  y2pix = 2 * x2pix

  return
end

;*****************************************************************************
