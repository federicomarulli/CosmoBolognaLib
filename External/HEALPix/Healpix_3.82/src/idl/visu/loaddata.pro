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
;==================================================================================================
pro loaddata, file_in, select_in,$
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
  SAVE=save,ONLINE=online,NESTED=nested_online,UNITS=units,COORD=coord,ROT=rot_ang,QUADCUBE=quadcube,LOG=log, $
  ERROR=error, FLIP=flip, POLARIZATION=polarization, FACTOR=factor_u, OFFSET = offset_u


print,'Warning : The HEALPix routine LOADDATA is not in service anymore.'
print,'    It has been renamed LOADDATA_HEALPIX to avoid conflict with '
print,'      an existing routine in the coyote IDL library.'
message,/info,'Please, use LOADDATA_HEALPIX instead of LOADDATA.'

loaddata_healpix, file_in, select_in,$
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
  SAVE=save,ONLINE=online,NESTED=nested_online,UNITS=units,COORD=coord,ROT=rot_ang,QUADCUBE=quadcube,LOG=log, $
  ERROR=error, FLIP=flip, POLARIZATION=polarization, FACTOR=factor_u, OFFSET = offset_u


return
end
