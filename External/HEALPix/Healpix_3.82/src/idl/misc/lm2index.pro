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
PRO LM2INDEX, l, m, index

;+
; NAME:
;       LM2INDEX
;
; PURPOSE:
;       to convert l,m mode indices in the spherical harmonic basis
;       to an integer index value.
;       index = l(l+1) + m + 1
;
; CALLING SEQUENCE:
;       LM2INDEX, l, m, index
; 
; INPUTS:
;       l     = (long) integer array
;       m     = (long) integer array
;
; OUTPUTS:
;       index = long integer array
;
; MODIFICATION HISTORY:
;       May 1999: written by A.J. Banday (MPA)       
;       Feb 2000: now works for large l, EH (Caltech)
;       Aug 2002: added error message, EH (Caltech)
;       Sep 2014: cosmetic changes
;
;-

if n_params() ne 3 then begin
    message,'syntax : lm2index, l, m, index'
endif


; compute {l,m} to index relation
index = long(l)^2 + l + m  + 1

; Exit routine
return
end



