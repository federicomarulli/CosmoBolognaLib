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
pro init_astrolib
;+
; trivial routine to call ASTROLIB routine (ASTRON library)
;
; 
; prompted by the removal of !DEBUG in Jan 2009
; now the variable !TEXTOUT is tested before calling ASTROLIB
;
; 2009-01-19: v1.0, EH, IAP
;-

defsysv, '!TEXTOUT', EXISTS = i  ; check if astrolib variables have been set-up
if (i ne 1) then astrolib         ; if not, run astrolib to do so

return
end
