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
;============================================================
pro planck_colors, option, get=get, help=help, show=show
;+
; NAME:
;    Planck_colors
;
; PURPOSE:
;    Creates Planck specific RGB color tables for CMB map or
;    frequency sky maps
;
; CATEGORY:
;    Data visualization
;
; CALLING SEQUENCE:
;    Planck_colors, option, [GET=, HELP=, SHOW=]
;
; INPUTS:
;    Option:   scalar integer in [1,2]
;      1 creates the 'parchment' Blue-red color table suitable for maps
;        dominated by Gaussian signal (eg, CMB)
;      2 creates a Blue-red-white color table suitable for maps with 
;        high dynamic signal (eg, Galactic foreground)
;
; KEYWORD PARAMETERS:
;    HELP:
;     if set, this documentation header is printed out and the routine exits
;
;    SHOW:
;     if set, the chosen color table is shown in a new window
;
; OPTIONAL OUTPUTS:
;    GET: 
;     on ouput, contains the newly created RGB color table
;     in a [256, 3] array
;
; EXAMPLES:
;
;   planck_colors, 1, /show ; defines and shows color table 1
;
;   mollview, cmb, col='-planck1'  ; applies reversed color table 1 to a CMB map
;
;   mollview, m217, col='planck2' ; applies color table 2 to a frequency map
;
; SIDE EFFECTS:
;   modify current color table using TVLCT
;
; MODIFICATION HISTORY:
;     Adapted from S. Colombi create_colortable1.
;     Replaced calls to MODIFYCT (which require write access 
;     to main IDL directories) with TVLCT which modify
;     the current color table
;
; v1.0: 2013-04-10, EH @ IAP
; v1.1: 2014-05-14: updated documentation header
; v1.2: 2014-05-16: added SHOW
;
; For more questions about the color tables, ask 
;   Stephane Colombi (colombiATiap.fr) for option 1
;   Krzysztof Gorski (krzysztof.m.gorskiATjpl.nasa.gov) for option 2
;-
;============================================================

routine = 'planck_colors'
syntax = routine+', option [,GET=, HELP=, SHOW=]'

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if (n_params() ne 1) then begin
    print,syntax
    return
endif

; define R, G, B on a small set of nodes
case option of
    1: begin
        Itab = [   0,  42,  85, 127, 170, 212, 255]
        Rtab = [   0,   0,   0, 255, 255, 255, 100]
        Gtab = [   0, 112, 221, 237, 180,  75,   0]
        Btab = [ 255, 255, 255, 217,   0,   0,   0]
    end
    2: begin
        Itab = [   0,  13,  26,  39,  52,  65,  76,  77,  88, $
                   101, 114, 127,   140, 153, 166, 179,   192, $
                   205, 218, 231, 255]
        Rtab = [   0,  10,  30,  80, 191, 228, 241, 241, 245, $
                   248, 249.9, 242.25,204, 165, 114, 127.5, 178.5, $
                   204, 229.5, 242.25, 252.45]
        Gtab = [   0,  20, 184, 235, 239, 240, 241, 241, 240, $
                   235, 204,   153,  76.5,  32,   0, 127.5, 178.5, $
                   204, 229.5, 242.25, 252.45]
        Btab = [ 255, 255, 255, 255, 250, 245, 212, 212, 175, $
                 130, 38.25, 12.75,  0,   32,  32, 153,   204,   $
                 229.5, 242.25, 249.9, 255]
    end
    else: begin
        print,syntax
        message,'option must be in [1,2]'
    end

endcase

; interpolate linearly between nodes

ncolors = 256 ; number of colors on output
ii = dindgen(ncolors)
R = byte(interpol(Rtab, Itab, ii))
G = byte(interpol(Gtab, Itab, ii))
B = byte(interpol(Btab, Itab, ii))

; implements new color table
tvlct, R, G, B

; put (R,G,B) in GET vector
get = [[R], [G], [B]]

; show color table
if keyword_set(show) then begin
    xsize= ncolors * 2	
    ysize= 100
    window, /free, xsize=xsize, ysize=ysize, title='Planck color table '+strtrim(option,2)
    tvscl, lindgen(xsize) #  replicate(1,ysize)
endif

return
end
