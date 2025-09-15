pro fwdcube, X,Y,XI,ETA
;+                                                                  
;  NAME:
;    fwd_cube
;
;  PURPOSE:                                   
;     BASED ON POLYNOMIAL FIT FOUND USING FCFIT.FOR
;     TAKEN FROM FORWARD_CUBE.FOR
;
;  CALLING SEQUENCE:
;    fwd_cube, x,y,ix,eta
;
;  INPUT:
;    X,Y IN RANGE -1 TO +1 ARE DATABASE CO-ORDINATES
;
;  OUTPUT:
;    XI, ETA IN RANGE -1 TO +1 ARE TANGENT PLANE CO-ORDINATES
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;    J.M. Gales		1/92		IDL version converted from FORTRAN
;-
;

	p = fltarr(29)		; p(0) not used

	p(1) =  -0.27292696
	p(2) =  -0.07629969
	p(3) =  -0.02819452
	p(4) =  -0.22797056
	p(5) =  -0.01471565
	p(6) =   0.27058160
	p(7) =   0.54852384
	p(8) =   0.48051509
	p(9) =  -0.56800938
	p(10) = -0.60441560
	p(11) = -0.62930065
	p(12) = -1.74114454
	p(13) =  0.30803317
	p(14) =  1.50880086
	p(15) =  0.93412077
	p(16) =  0.25795794
	p(17) =  1.71547508
	p(18) =  0.98938102
	p(19) = -0.93678576
	p(20) = -1.41601920
	p(21) = -0.63915306 
	p(22) =  0.02584375
	p(23) = -0.53022337
	p(24) = -0.83180469
	p(25) =  0.08693841 
	p(26) =  0.33887446
	p(27) =  0.52032238
	p(28) =  0.14381585

	XX=X*X
	YY=Y*Y

	XI=X*(1.+(1.-XX)*( $
	P(1)+XX*(P(2)+XX*(P(4)+XX*(P(7)+XX*(P(11)+XX*(P(16)+XX*P(22)))))) + $
	YY*( P(3)+XX*(P(5)+XX*(P(8)+XX*(P(12)+XX*(P(17)+XX*P(23))))) + $
	YY*( P(6)+XX*(P(9)+XX*(P(13)+XX*(P(18)+XX*P(24)))) + $
	YY*( P(10)+XX*(P(14)+XX*(P(19)+XX*P(25))) + $
	YY*( P(15)+XX*(P(20)+XX*P(26)) + $
	YY*( P(21)+XX*P(27) + YY*P(28))))) )))

	ETA=Y*(1.+(1.-YY)*( $
	P(1)+YY*(P(2)+YY*(P(4)+YY*(P(7)+YY*(P(11)+YY*(P(16)+YY*P(22)))))) + $
	XX*( P(3)+YY*(P(5)+YY*(P(8)+YY*(P(12)+YY*(P(17)+YY*P(23))))) + $
	XX*( P(6)+YY*(P(9)+YY*(P(13)+YY*(P(18)+YY*P(24)))) + $
	XX*( P(10)+YY*(P(14)+YY*(P(19)+YY*P(25))) + $
	XX*( P(15)+YY*(P(20)+YY*P(26)) + $
	XX*( P(21)+YY*P(27) + XX*P(28))))) )))

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


