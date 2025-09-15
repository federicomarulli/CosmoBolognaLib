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
function decode_coord, coord_code, nlong=nlong, nlat=nlat, flong=flong, flat = flat, error = error
;+
; NAME:
;   decode_coord
;
; PURPOSE:
;   returns the meaning of the coordinate coding used for Cobe
;
; CATEGORY:
;
; CALLING SEQUENCE:
;  result = DECODE_COORD(Coord_code, [Nlong=, Nlat=, Flong=, Flat=, Error=])
; 
; INPUTS:
;  Coord_code = 1 character string coding for coordinate system
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
;  result = string with the full name of the coordinate system
;
; OPTIONAL OUTPUTS:
;  Nlong, Nlat : named variable containing on output
;     the longitude-like and latitude-like coordinates
;  Flong, Float : named variable containing on output
;     the Fits coding for longitude and latitude coordinates
;  Error : set to 1 on output if invalid input
;
;  Coord_code         result   Nlong   Nlat    Flong      Flat   Slong   Slat
;         'E'     'Ecliptic'  'Long'  'Lat'   'ELON-'   'ELAT-'
;         'Q'   'Equatorial'    'RA'  'Dec'   'RA---'   'DEC--'   'Alpha' 'Delta'
;         'C'    same as for 'Q'
;         'G'     'Galactic'  'Long'  'Lat'   'GLON-'   'GLAT-'       'l'        'b'
;       other      'Unknown'  '????'  '???'   '-----'   '-----'
;      
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
; MODIFICATION HISTORY:
;         March 1999,   EH Caltech, version 1.0
;         Feb   2000,  added error variable
;         Feb   2002,  added missing '-' in flong and flat (pointed
;            out by P.R.McCullough)
;         Jan 2013,      EH IAP, 
;               treat 'C' as 'Q'
;               defines all output values for unknown coordinates
;-


coordinate ='unknown'
error = 1
if N_params() ne 1 then begin
    print,'incorrect call to DECODE_COORD'
    print, 'result = DECODE_COORD(Coord_code, [Nlong=, Nlat=, Flong=, Flat=])'
    return, coordinate
ENDIF

code = STRMID(STRTRIM(STRUPCASE(coord_code),2),0,1) ; first letter, upper case

case 1 of
    (code eq 'G') : begin
        coordinate = 'Galactic'
        nlong = 'Long'  & nlat = 'Lat'
        flong = 'GLON-' & flat = 'GLAT-'
        error = 0
        end
    (code eq 'E') : begin
        coordinate = 'Ecliptic'
        nlong = 'Long'  & nlat = 'Lat'
        flong = 'ELON-' & flat = 'ELAT-'
        error = 0
        end
    (code eq 'Q' || code eq 'C') : begin
        coordinate = 'Equatorial'
        nlong = 'RA'    & nlat = 'Dec'
        flong = 'RA---' & flat = 'DEC--'
        error = 0
        end
    else : begin
        print,'coordinate code unknown'
        coordinate = 'Unknown'
        nlong = '????'  & nlat = '???'
        flong = '-----' & flat = '-----'
        error = 1
    end
endcase

return,coordinate

end

