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
function euler_matrix_new, a1, a2, a3, $
  DEG=deg, $
  HELP=help, $
  X=x, Y=y, ZYX=zyx

;+
; NAME:
;         EULER_MATRIX_NEW
;
; PURPOSE:
;         computes the Euler matrix of an arbitrary rotation described
;         by 3 Euler angles
;         correct bugs present in Euler_Matrix
;        
; CALLING SEQUENCE:
;         result = euler_matrix_new (a1, a2, a3 [, DEG=, HELP=, X=, Y=, ZYX= ])
; 
; INPUTS:
;         a1, a2, a3 = Euler angles, scalar
;	       (in radian by default, in degree if DEG is set)
;	       all the angles are measured counterclockwise
;      
;         correspond to x, y, zyx-conventions (see Goldstein)
;        the default is x
;
; KEYWORD PARAMETERS:
;       DEG : if set, the angles are measured in degrees
;
;       HELP: if set, the documentation header is printed
;
; 	X : 	rotation a1 around original Z 
;     		rotation a2 around interm   X 
;     		rotation a3 around final    Z
;	 DEFAULT,  classical mechanics convention
;           same as: rotation a3 around Z,
;                    rotation a2 around initial (unrotated) X,
;                    rotation a1 around initial (unrotated) Z.
;
;
; 	Y : 	rotation a1 around original Z
;     		rotation a2 around interm   Y
;     		rotation a3 around final    Z
;	 quantum mechanics convention
;           same as: rotation a3 around Z,
;                    rotation a2 around initial (unrotated) Y,
;                    rotation a1 around initial (unrotated) Z.
;
;
; 	ZYX : 	rotation a1 around original Z
;     		rotation a2 around interm   Y
;     		rotation a3 around final    X
;	 aeronautics convention
;           same as: rotation a3 around X,
;                    rotation a2 around initial (unrotated) Y,
;                    rotation a1 around initial (unrotated) Z.
;
;       * these last three keywords are obviously mutually exclusive *
;
; OUTPUTS:
;       result is a 3x3 matrix
;
; USAGE: if vec is an Nx3 array containing N 3D vectors,
;    vec2 = vec # euler_matrix_new(a1,a2,a3,/Y) will be the rotated vectors
;
;    alternatively
;   
;    vec2 = rotate_coord(vec, euler_matrix=euler_matrix_new(a1,a2,a3,/Y) )
;        will also produce the rotated vectors
;
;
; MODIFICATION HISTORY:
;	March 2002, EH, Caltech, rewritting of euler_matrix
;       Dec   2012, Paddy Leahy, Manchester: use double precision d2r if needed.
;       2014-08: EH, HELP keyword added
;       2016-11: EH, added description of rotations around initial axes
;
;  convention   euler_matrix_new           euler_matrix
;      X:       M_new(a,b,c,/X)  =  M_old(-a,-b,-c,/X) = Transpose( M_old(c, b, a,/X))
;      Y:       M_new(a,b,c,/Y)  =  M_old(-a, b,-c,/Y) = Transpose( M_old(c,-b, a,/Y))
;    ZYX:       M_new(a,b,c,/Z)  =  M_old(-a, b,-c,/Z)
;
;-

routine = 'euler_matrix_new'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if (n_params() ne 3) then begin
    message,'  Invalid number of arguments ',/noprefix,/inform
    message,' Syntax : result = '+routine+'(a1, a2, a3, [HELP=, DEG=, X=, Y=, ZYX=])',/noprefix,/noname
endif

t_k = 0
IF KEYWORD_SET(zyx) THEN t_k +=1
IF KEYWORD_SET(x  ) THEN t_k +=1
IF KEYWORD_SET(y  ) THEN t_k +=1
IF t_k GT 1 THEN BEGIN
	message, 'incompatible keywords choice',/noprefix
ENDIF

double = SIZE(a1,/TYPE) EQ 5 || SIZE(a2,/TYPE) EQ 5 || SIZE(a3,/TYPE) EQ 5
convert = 1.0
IF KEYWORD_SET(deg) THEN convert = double ? !dpi/180d : !DtoR

c1 = COS(a1*convert)
s1 = SIN(a1*convert)
c2 = COS(a2*convert)
s2 = SIN(a2*convert)
c3 = COS(a3*convert)
s3 = SIN(a3*convert)


if (keyword_set(ZYX)) then begin

    m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
    m2 = [[ c2,  0, s2],[  0,  1,  0],[-s2,  0, c2]] ; around   y
    m3 = [[  1,  0,  0],[  0, c3,-s3],[  0, s3, c3]] ; around   x

endif else if (keyword_set(Y)) then begin

    m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
    m2 = [[ c2,  0, s2],[  0,  1,  0],[-s2,  0, c2]] ; around   y
    m3 = [[ c3,-s3,  0],[ s3, c3,  0],[  0,  0,  1]] ; around   z

endif else begin

    m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
    m2 = [[  1,  0,  0],[  0, c2,-s2],[  0, s2, c2]] ; around   x
    m3 = [[ c3,-s3,  0],[ s3, c3,  0],[  0,  0,  1]] ; around   z

endelse

; A ## B = matrix product AB
; A # B  = matrix product BA
; vec # A = product A.vec  if A is 3x3 and vec is nx3


euler_matrix = m3 # ( m2 # m1) ; m1 m2 m3


return, euler_matrix
end

