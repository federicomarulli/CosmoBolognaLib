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
function angulardistance, v1, v2, help=help
;+
; NAME:
;        AngularDistance
;
; PURPOSE:
;        Computes the angular distance (in RADIANS) between 2 3D vectors,
;        or lists of vectors.
;
; CATEGORY:
;        3D Geometry
;
; CALLING SEQUENCE:
;        distance = angulardistance(V, W [, HELP=])
;
; INPUTS:
;        V: single 3D vector (of shape (3) or (1,3)) 
;           or lists of 3D vectors (of shape (n,3))
;
;        W: single 3D vector (of shape (3) or (1,3)) 
;           or lists of 3D vectors (of shape (n,3))
;
;        It is not necessary for V and W vectors to be normalised to 1
;        upon calling the function.
;
;        Note: if V and W both are lists of vectors, 
;        they should be of the same length
;
;        If V (and/or W) has the form (n,3,4) 
;        (like the pixel corners returned by pix2vec\_*), 
;        it should be preprocessed with
;        V = reform( transpose(V, [0,2,1]), n_elements(V)/3, 3)
;        to take the form (n*4,3) before being passed to angulardistance.  
;
;
; OPTIONAL INPUTS:
;        None
;
; KEYWORD PARAMETERS:
;        Help= : returns this documentation header and exits
;
; OUTPUTS:
;        distance: angular distance in RADIANS
;        - if both V and W are single vectors, then 
;          distance = dist(V, W)
;        - if V (or W) is a single vector, while W (or V) is a list, then
;         distance[i] = dist(V[i], W)
;        - if both V and W are lists of vectors (of same length), then
;         distance[i] = dist(V[i], W[i])
;
; OPTIONAL OUTPUTS:
;       None
;
; COMMON BLOCKS:
;       None
;
; SIDE EFFECTS:
;       None
;
; RESTRICTIONS:
;
;     If V (and/or W) has the form (n,3,4) (like the pixel corners returned by
;     pix2vec_*), it should be preprocessed with
;     V = reform( transpose(V, [0,2,1]), n_elements(V)/3, 3)
;      before being passed to this routine.
;
; PROCEDURE:
;
;        V and W are internally renormalized, and then
;        distance = acos ( V . W )  in general
;        distance = 2 * asin ( |V-W|/2), when V and W are almost aligned
;
;
; EXAMPLE:
;      nside=8
;      pix2vec_ring, nside, lindgen(nside2npix(nside)), vpix
;      mollview, angulardistance( vpix, [1,1,1])
;
;    will plot the angular distance between the Healpix pixel center for
;    Nside=8, and the vector (x,y,z) = (1,1,1)/sqrt(3)
;
;
; MODIFICATION HISTORY:
;    Oct 2003: EH, first version of angdist routine
;    EH, corrected bug on scalar product boundary
;    2010-03-01: works for non-normalized input vectors
;    2011-04-12: turned in angulardistance function which accepts lists of
;    vectors
;
;-


routine = 'angulardistance'
syntax = 'dist = '+routine+'(V, W [,HELP=] )'

if (keyword_set(help)) then begin
    doc_library,routine
    return, -1
end

if (n_params() ne 2) then begin
    print,syntax
    return, -1
endif

; make sure vectors have shape (n, 3)
nn1 = n_elements(v1) & n1 = nn1/3 & sz1 = size(v1)
nn2 = n_elements(v2) & n2 = nn2/3 & sz2 = size(v2)
if (n1*3 ne nn1 || n2*3 ne nn2) then begin
    print,syntax
    message,'vectors must be of shape (3) or (n,3)'
endif
if (sz1[0] gt 2 || sz2[0] gt 2 || (sz1[0] eq 2 && sz1[2] ne 3) || (sz2[0] eq 2 && sz2[2] ne 3)) then begin
    print,syntax
    print, 'V has shape: ',sz1[1:sz1[0]]
    print, 'W has shape: ',sz2[1:sz2[0]]
    message,'vector lists must be of shape (n,3)'
endif
if (n1 gt 1 && n2 gt 1 && n1 ne n2) then begin
    print,syntax
    print, 'size of lists:', n1, n2
    message,'if V and W are both lists of vectors, they should be of the same length'
endif
one = replicate(1.d0, n1>n2)
r1 = (n1 eq 1) ? one # reform(v1, 3) : double(v1)
r2 = (n2 eq 1) ? one # reform(v2, 3) : double(v2)
; normalize both vectors
one = replicate(1.d0, 3)
r1 /= (sqrt(total(r1*r1, 2)) # one)
r2 /= (sqrt(total(r2*r2, 2)) # one)

; scalar product
sprod = total(r1*r2, 2, /double)
    
; standard calculation
dist = acos ( sprod )

cut = 0.9d0

; almost colinear vectors
p1 = where(sprod gt cut,    n_p1)
if (n_p1 gt 0) then begin
    vdiff = r1[p1,0:2] - r2[p1,0:2]
    diff = sqrt(total(vdiff*vdiff, 2)) ; norm of difference
    dist[p1] = 2.0d0 * asin(diff * 0.5d0)
endif

; almost anti-colinear vectors
m1 = where(sprod lt (-cut), n_m1)
if (n_m1 gt 0) then begin
    vdiff = r1[m1,0:2] + r2[m1,0:2]
    diff = sqrt(total(vdiff*vdiff, 2)) ; norm of sum
    dist[m1] = !DPI - 2.0d0 * asin(diff * 0.5d0)
endif

return, dist
end


