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
pro preview_file, file, ps=ps, gif=gif, png=png, jpeg=jpeg, landscape=landscape, pdf=pdf
;+
; NAME:
;     preview_file
;
; PURPOSE:
;     previews Gif, of Png or Postscript files created by Healpix
;     visualisation routines
;
; CATEGORY:
;     visualisation
;
;
; CALLING SEQUENCE:
;     preview_file, file, GIF=, JPEG=,  LANDSCAPE=, PNG=, PS=, PDF=
;
; INPUTS:
;     file
;
; OPTIONAL INPUTS:
;
; 
; KEYWORD PARAMETERS:
;        GIF
;        JPEG
;        LANDSCAPE
;        PNG
;        PS
;
; OUTPUTS:
;        None
;
; OPTIONAL OUTPUTS:
;        None
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;        Spawn an external code
;
; RESTRICTIONS:
;       Requires the definition file idl_default_previewer
;
; PROCEDURE:
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;        Eric Hivon, Dec 2002, McM
;        May 2005, define 'factory' settings and issue reminder
;        Mar 2006: added comments before crash
;-


; ****************************************************************

; papersize : paper size used
; media : either paper size or blank (in which case BBox is used)
; ps_scom : short name of the PS viewing command
; ps_com  : full name of the PS viewing command
; gif_scom : short name of the PNG/GIF viewing command
; gif_com  : full name of the PNG/GIF viewing command

; default visualisation utilities
papersize = 'a4'
media = ' '
ps_scom = 'gv'
ps_com  = ps_scom
pdf_com = ps_scom
gif_scom = 'display'
gif_com = gif_scom
settings = 'factory' ; that's how it feels

comments=["-----------------------------------------------------------------------------",$
          "You can choose the facilities used to visualize PDF, Postscript, PNG, GIF, JPEG files",$
          "and the hard copy paper size,",$
          "by running the configure script in the main Healpix directory.",$
;          "           (no need to restart IDL ;-)",$
          "-----------------------------------------------------------------------------"]

; 2010-May-27
test_preview, found_preview
@idl_default_previewer.pro

; issue a friendly but annoying reminder if user did not customize settings
if (settings eq 'factory') then print,comments,form='(a)'

do_bbox = 0
if defined(media) then begin
    do_bbox = (strcompress(media,/remove_all) eq '')
endif

fullfile = file
blank=' '

; POSTSCRIPT
; ==========

case ps_scom of
    'gv' : begin
    ; gv
        ps_pps = '  -media '+papersize+blank
        ps_lsc = ' -seascape ' 
        ps_ext = ' ' ; extra options
    end
    'ghostview ' : begin
    ; ghostview
        ps_pps = '  -'+papersize+blank
        ps_lsc = ' -seascape' 
        ps_ext = ' '
    end
    'gs' : begin
    ; ghostscript
        ps_pps = ' -sPAPERSIZE='+papersize+blank
        ps_lsc = ' '
        ps_ext = ' '
    end
    else : begin
    ; unknown : neutral options
        ps_pps = ' '
        ps_lsc = ' '
        ps_ext = ' '
    end
endcase

if do_bbox then ps_pps = ' ' ; use default BBox
ps_opt = ps_ext + ps_pps
;;;print,ps_opt+ps_lsc


; PDF
; ===

case pdf_scom of
    else : begin
    ; unknown : neutral options
        pdf_pps = ' '
        pdf_lsc = ' '
        pdf_ext = ' '
    end
endcase

if do_bbox then pdf_pps = ' ' ; use default BBox
pdf_opt = pdf_ext + pdf_pps
;;;print,ps_opt+ps_lsc

; PNG and GIF
; ===========

; gif_com = 'xv '
; gif_com = 'gimp '
; gif_com = 'kpaint '
; gif_com = 'kview '
; gif_opt = ' '
if (keyword_set(gif) || keyword_set(png) || keyword_set(jpeg)) then begin
    if (gif_scom eq 'netscape' or gif_scom eq 'mozilla') then begin
        cd, current=cwd         ; find full path to current directory
        fullpath = filepath(file,root=cwd) ; build full path to file (OS dependent)
        if (1 - file_test(fullpath)) then fullpath=file ; file is already a full path
        fullfile = 'file:'+fullpath
        gif_opt = ' '
    endif else begin
        gif_opt = ' '
    endelse
endif


; background task
; ===============

os_family = !version.os_family
case os_family of
    'unix' : begin
        bg = ' & '
    end
    'vms' : begin
        bg = '   ' ; don't know
    end
    'windows' : begin
        bg = '   ' ; don't know
    end
    else : begin
        bg = '   '
    end
endcase

; **************************************************************
;   construct command
; **************************************************************

; POSTSCRIPT
if (keyword_set(ps)) then begin
    if ( ~ keyword_set(landscape)) then ps_lsc = ' '
    command = ps_com + ps_opt + ps_lsc + blank + fullfile + bg
endif

; PDF
if (keyword_set(pdf)) then begin
    command = pdf_com + pdf_opt + blank + fullfile + bg
endif

; PNG and GIF
if (keyword_set(gif) || keyword_set(png) || keyword_set(jpeg)) then begin
    command = gif_com + blank + gif_opt + blank + fullfile + bg
endif

; **************************************************************
;   execute command
; **************************************************************

if defined(command) then spawn, command

return
end

