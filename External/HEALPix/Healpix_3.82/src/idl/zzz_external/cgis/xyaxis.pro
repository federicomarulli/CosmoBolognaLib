pro XYAXIS, NFACE,XI,ETA,C
;+                                                                  
;  NAME:
;    xyaxis
;
;  PURPOSE:
;    CONVERTS FACE NUMBER NFACE (0-5) AND XI, ETA (-1. - +1.)
;    INTO A UNIT VECTOR C
;
;  CALLING SEQUENCE:
;    xyaxis, nface,xi,eta,c
;
;  INPUT:
;    nface - sky-cube face number
;    xi, eta - tangent place coordinates
;
;  OUTPUT:
;    c - unit vector
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;    A.C. Raugh		8/90		Original FORTRAN version
;    J.M. Gales		1/92		IDL version
;    J.M. Gales		4/92		Eliminated CASE statement
;					Improved norm calculation
;-
;

	n = n_elements(nface)

	nface_0 = fix(nface EQ 0)
	nface_1 = fix(nface EQ 1)
	nface_2 = fix(nface EQ 2)
	nface_3 = fix(nface EQ 3)
	nface_4 = fix(nface EQ 4)
	nface_5 = fix(nface EQ 5)

row0 = eta * (nface_5 - nface_0) + $
        xi * (nface_4 - nface_2) + $
	     (nface_1 - nface_3)

row1 =  xi * (nface_0 + nface_1 - nface_3 + nface_5) + $
             (nface_2 - nface_4)

row2 = eta * (nface_1 + nface_2 + nface_3 + nface_4) + $
             (nface_0 - nface_5)

norm = sqrt(1 + xi*xi + eta*eta)

row0 = row0 / norm
row1 = row1 / norm
row2 = row2 / norm

c = [[row0],[row1],[row2]]

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


