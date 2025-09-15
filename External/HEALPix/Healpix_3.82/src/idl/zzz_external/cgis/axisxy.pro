;******************************************************************************
;
;******************************************************************************
PRO AXISXY,C,NFACE,X,Y
;C
;C	CONVERTS UNIT VECTOR C INTO NFACE NUMBER (0-5) AND X,Y
;C	IN RANGE 0-1
;C
;
; New algorithm developed by J.M Gales to eliminate CASE statement
; Apr 92
;
n = n_elements(c)/3

c0 = c(*,0)
c1 = c(*,1)
c2 = c(*,2)

abs_yx = ABS(c1/c0)
abs_zx = ABS(c2/c0)
abs_zy = ABS(c2/c1)

nface = 0 * (abs_zx GE 1 AND abs_zy GE 1 AND c2 GE 0) + $
        5 * (abs_zx GE 1 AND abs_zy GE 1 AND c2 LT 0) + $
        1 * (abs_zx LT 1 AND abs_yx LT 1 AND c0 GE 0) + $
        3 * (abs_zx LT 1 AND abs_yx LT 1 AND c0 LT 0) + $
        2 * (abs_zy LT 1 AND abs_yx GE 1 AND c1 GE 0) + $
        4 * (abs_zy LT 1 AND abs_yx GE 1 AND c1 LT 0)
			; determine face number

nface_0 = fix(nface EQ 0)
nface_1 = fix(nface EQ 1)
nface_2 = fix(nface EQ 2)
nface_3 = fix(nface EQ 3)
nface_4 = fix(nface EQ 4)
nface_5 = fix(nface EQ 5)
			; set logical vectors

eta = (c2/c0) * (nface_1 - nface_3) + $
      (c2/c1) * (nface_2 - nface_4) - $
      (c0/c2) * (nface_0 + nface_5)
			; calculate eta

xi  = (c1/c0) * (nface_1 + nface_3) - $
      (c0/c1) * (nface_2 + nface_4) + $
      (c1/c2) * (nface_0 - nface_5)
			; calculate xi

INCUBE,XI,ETA,X,Y
X=(X+1.)/2.
Y=(Y+1.)/2.

RETURN
END
;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.


