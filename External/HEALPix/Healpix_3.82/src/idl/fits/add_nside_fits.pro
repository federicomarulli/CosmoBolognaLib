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
pro add_nside_fits, info_header, nside=nside, partial=partial, error=error, obs_npix = obs_npix
;+
; add_nside_fits
;
; add_nside_fits,info_header, nside=, partial=, error=
;
;  adds the NSIDE keyword defined by 'nside' in
;  the fits header 'info_header'
;
;  a full sky coverage is assumed unless partial is set
;
;  EH 2000-02
;  2000-11 : added INDXSCHM and GRAIN keywords
;  2001-07 : corrected bug
;  2005-02-08 : remove creation of header primer if no info_header provided
;               insure comments on grain correctly placed
;  2020-02-14 : removed GRAIN
;-

add_grain = 0
error=0

npix = nside2npix(nside,err=errpix)
if errpix gt 0 then begin
    print,'Invalid nside ',nside
    error=1
    return
endif

if (add_grain) then begin
    old_grain = sxpar(info_header,'GRAIN',count=cgrain)
endif

if keyword_set(partial) then begin
    grain = 1
    sxaddpar,info_header,'NSIDE',nside,' Healpix resolution parameter'
    sxdelpar,info_header,['NPIX','OBJECT','FIRSTPIX','LASTPIX']
    sxaddpar,info_header,'OBJECT','PARTIAL',' Sky coverage, either FULLSKY or PARTIAL'
    sxaddpar,info_header,'INDXSCHM','EXPLICIT',' indexing : IMPLICIT or EXPLICIT'
    if (add_grain) then begin
        sxaddpar,info_header,'GRAIN',grain,' GRAIN = 0: No index,'
        if (cgrain eq 0) then begin
            sxaddpar,info_header,'COMMENT','        GRAIN >1: 1 pixel index for Grain consecutive pixels',after='GRAIN'
            sxaddpar,info_header,'COMMENT','        GRAIN =1: 1 pixel index for each pixel,',after='GRAIN'
        endif
    endif
    if defined(obs_npix) then sxaddpar,info_header, 'OBS_NPIX', obs_npix, ' Number of pixel observed and recorded',after='OBJECT'
endif else begin
    grain = 0
    sxaddpar,info_header,'NSIDE',nside,' Healpix resolution parameter'
    sxaddpar,info_header,'NPIX',npix,' Total number of pixels',after='NSIDE'
    sxaddpar,info_header,'OBJECT','FULLSKY',' Sky coverage, either FULLSKY or PARTIAL'
    sxaddpar,info_header,'FIRSTPIX',0,' First pixel # (0 based)'
    sxaddpar,info_header,'LASTPIX',npix-1,' Last pixel # (zero based)'
    sxaddpar,info_header,'INDXSCHM','IMPLICIT',' indexing : IMPLICIT or EXPLICIT'
    if (add_grain) then begin
        sxaddpar,info_header,'GRAIN',grain,' GRAIN = 0: No index,'
        if (cgrain eq 0) then begin
            sxaddpar,info_header,'COMMENT','        GRAIN >1: 1 pixel index for Grain consecutive pixels',after='GRAIN'
            sxaddpar,info_header,'COMMENT','        GRAIN =1: 1 pixel index for each pixel,',after='GRAIN'
        endif
    endif
endelse

return
end

;--------------------------------
