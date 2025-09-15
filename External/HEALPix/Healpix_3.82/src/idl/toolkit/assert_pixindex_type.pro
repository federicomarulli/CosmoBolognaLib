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
pro assert_pixindex_type, pixel, warning=warning
;+
;
; assert_pixindex_type, pixel, [warning=]
;
;  returns an error or a warning if pixel is not of integer type
;
;
;  used by pix2ang_*, pix2vec_*, ring2nest, nest2ring, template_pixel_*
;
; EH, IAP, Aug 2008
;-

type = size(/type, pixel[0])

if (type ge 4 and type le 11) then begin
                                ; exclude float, double, (d)complex, string,
                                ; structure, pointer and object reference
    if keyword_set(warning) then begin
        message_patch,level=-1,/info,'***************** WARNING ************************'
    endif else begin
        message_patch,level=-1,/info,'**************************************************'
    endelse
    message_patch,level=-1,/info,'Pixel index must be a (1, 2,) 4 or 8-byte INTEGER,'
    message_patch,level=-1,/info,'ie, either (BYTE, INT,) LONG or LONG64,'
    message_patch,level=-1,/info,'any other type may create round-off errors.'
    message_patch,level=-1,/info,'It is currently: '+size(/tname,pixel[0])
    message_patch,level=-1,info=warning,'**************************************************'
    if (keyword_set(warning)) then begin
        if is_gdl() then begin
            help, call=trace
            print, trace,form='(a)'
        endif else begin
            help,/traceback
        endelse
    endif
endif



return
end
