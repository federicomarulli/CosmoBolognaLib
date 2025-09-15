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
function rotate_coord, in_uvec, $
                       inco=in_coord, $
                       outco=out_coord, $
                       euler_matrix=eul_mat, $
                       stokes_parameters = stokes, $
                       delta_psi = dpsi, $
                       help=help, $
                       Free_Norm = Free_Norm         ;,junk=junk

;+
; NAME:
;   rotate_coord
;
;
; PURPOSE:
;   apply change of astronomical coordinates OR active rotation to
;   position vectors and Stokes parameters Q and U
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;    Result= rotate_coord(Vec, Delta_Psi=, Euler_Matrix=, Inco=, Outco=, Stokes_Parameters=, 
;                             /Free_Norm, /Help )
;
; 
; INPUTS:
;    Vec : a (N x 3) array of coordinates
;
;
; OPTIONAL INPUTS:
;    Delta_Psi: vector of size N containing on output the change in azimuth
;       in Radians resulting from the rotation
;       (measured with respect to the local meridian, from South to East) 
;       so that for a field of spin s the output Q',U' is related to the input Q,U via
;                  Q' = Q cos(s psi)  -  U sin (s psi)
;                  U' = U cos(s psi)  +  Q sin (s psi)
;       with s=2 for polarization Stokes parameters
;
;    Euler_matrix : 3x3 matrix of the active rotation
;      (see eg, Euler_matrix_new routine)
;
;    Inco : input coordinate system (must be used with Outco)
;     'G' :     Galactic
;     'E' :     Ecliptic
;     'C','Q' : Celestial, aka eQuatorial
;
;    Outco : output coordinate system (must be used with Inco)
;      one of 'G', 'E', 'C', 'Q'
;
;    Stokes_Parameters : (N x 2) array
;     * On input: contains the values of Stokes parameters Q and U 
;      at the location of the input vector vec;
;      where Q measures the polarisation along or across the meridian, 
;      and U the polarisation at 45deg of the meridian.
;     * On output: contains the values of Q and U resulting from the rotation
;
;          
;      
; KEYWORDS:
;     /Free_Norm: if set (and Stokes_Parameters and/or Delta_Psi are present)
;        the input (and output) coordinate vectors are *not* assumed to be normalized to 1.
;        Using this option is therefore safer, but 20 to 30% slower.
;        (3D vectors produced by ang2vec, pix2vec_nest and pix2vec_ring *are* properly normalized).
;        Ignored when Stokes_parameters and Delta_Psi are both absent.
;
;    /Help:  if set, prints the documentation header and exits
;
;      
;
; OUTPUTS:
;    Result : a (N x 3) array of coordinates, with the same norm as the input ones
;
;
; OPTIONAL OUTPUTS:
;   Stokes_Parameters is modified on output
;   Delta_Psi         can be output
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;    Only one of the rotations (astrophysical coordinate system or
;    active rotation) can be applied at a time
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;    ???: EH, creation
;    2014-08-25: added Help keyword
;    2018-07-20: addition of Delta_Psi and Free_Norm; fixed internal sign of psi
;-

routine = 'rotate_coord'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if (n_params() ne 1) then begin
    message,'Invalid number of arguments',/noprefix,/infor
    message,'result= '+routine+'(Vec, [Euler_Matrix=, HELP=, Inco=, Outco=, Stokes_parameters=, Delta_Psi=, Free_Norm=])',/noprefix,/noname

endif

if ((defined(in_coord) or defined(out_coord)) and defined(eul_mat)) then begin
    message,' can not apply the active rotation and the change of coordinate system at the same time',/noprefix
endif

vec0 = transpose([0.d0,0.d0,1.d0]) ; north pole of initial coordinate system

nv = (size(in_uvec,/dim))[0]
if (defined(in_coord) and defined(out_coord)) then begin
    ic = strupcase(strmid(in_coord, 0,1)) ; 1st letter, uppercase
    oc = strupcase(strmid(out_coord,0,1))

    loadsky
    if ic eq 'C' then ic = 'Q'
    if oc eq 'C' then oc = 'Q'
    vec = skyconv(in_uvec, inco=ic, outco=oc)

    vec0 = skyconv(vec0, inco=oc, outco=ic)

endif else begin
    vec = in_uvec

    vec0 = vec0

endelse

if (defined(eul_mat)) then begin
    sz = size(eul_mat)
    if ((sz[1] ne sz[2]) or (sz[0] ne 2) or (sz[1] ne 3)) then begin
        message,'The Euler matrix should be 3x3, it is : '+string(sz[0:2],form='(3(i5))'),/noprefix
    endif

    vec  = vec  # eul_mat
    vec0 = vec0 # transpose(eul_mat)

endif


do_psi    = arg_present(dpsi)
do_stokes = defined(stokes)

if (do_stokes || do_psi) then begin
    do_renorm = keyword_set(Free_Norm)
    if (do_stokes) then begin
        ;ns = n_elements(stokes[*,0])
        ns = (size(stokes,/dim))[0]
        if (ns ne nv) then begin
            message,'Number of Stokes parameter pairs (Q,U) should match number of 3D vectors',/noprefix,/inform
            message,string([ns,nv],form='(2(i10))'),/noprefix
        endif
    endif
    
    ; p = initial North Pole
    ; e_phi   =  p X r      / sqrt(x^2+y^2) : unit vector along parallel, pointing East
    ; e_theta = (p X r) X r / sqrt(x^2+y^2) : unit vector along meridian, pointing South
    ; R(e_phi), R(e_theta) : vectors after rotation
    ; R(e) . p = e . Rinv(p)
    ;
    ; R(r) = (x',y',z')
    ; Rinv(p) = (x0', y0', z0')
    ;  abs( p X r ) = sqrt(1-z*z)
    ;   Rinv(p) . r = p . R(r) = z'
    ;
    ; sin(psi) = (p X r) . Rinv(p)       * C
    ;          = -y*x0' + x*y0'          * C
    ; cos(psi) = ((p X r) X r). Rinv(-p) * C
    ;          = z0'(r.r) -z'z           * C
    ;
    ; with C = 1 / sqrt( (x^2+y^2) * (x'^2+y'^2) )

    ; not normalised SIN(psi) and COS(psi)
    if (do_renorm) then begin
        r2 = total(in_uvec^2, 2, /double) ; x^2 + y^2 +z^2
    endif else begin
        r2 = 1.d0
    endelse
    sin_psi = (vec0[0,1]*in_uvec[*,0] - vec0[0,0]*in_uvec[*,1])*sqrt(r2)
    cos_psi =  vec0[0,2]*r2           - vec [*,2]*in_uvec[*,2]
;     if arg_present(junk) then begin
;         print, 'compute junk'
;         fudge_i = sqrt(in_uvec[*,0]^2 + in_uvec[*,1]^2)
;         fudge_o = sqrt(    vec[*,0]^2 +     vec[*,1]^2)
;         junk = [[sin_psi/fudge_i/fudge_o], [cos_psi/fudge_i/fudge_o]]
;     endif

    if (do_psi) then begin
        dpsi = atan( sin_psi, cos_psi)
    endif

    if (do_stokes) then begin
        ; compute SIN(2*psi) and COS(2*psi)
        if (do_psi) then begin
;       --- slower method, more accurate ?, less memory
            s2psi = sin(2.d0 * dpsi)
            c2psi = cos(2.d0 * dpsi)
        endif else begin
;        --- faster method, less accurate ?, more memory
            norm    = sin_psi * sin_psi + cos_psi * cos_psi
            s2psi   =             2.0d0 * (sin_psi * cos_psi  / norm)
            c2psi   = (cos_psi + sin_psi)*(cos_psi - sin_psi) / norm
            cos_psi = 0.        ; free memory
            sin_psi = 0.        ; free memory
            norm    = 0.        ; free memory
        endelse

        ; new Q and U
        q = stokes[*,0] * c2psi - stokes[*,1] * s2psi
        u = stokes[*,1] * c2psi + stokes[*,0] * s2psi

        stokes = [ [q], [u] ]
        q = 0. & u = 0.         ; free memory
    endif

endif


return,vec
end

