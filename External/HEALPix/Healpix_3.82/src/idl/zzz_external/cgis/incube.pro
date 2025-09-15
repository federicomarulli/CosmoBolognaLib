;******************************************************************************
;
;******************************************************************************
PRO INCUBE,ALPHA,BETA,X,Y
;C
;C	SEE CSC "EXTENDED STUDY ... NOTE THAT THE TEXT HAS TYPOS.  I HAVE
;C	TRIED TO COPY THE FORTRAN LISTINGS
;C
;
;       Copied, converted to IDL
;         from CSDR$SOURCE:[upx]upx_incube.for
;
	GSTAR=1.37484847732
	G=-0.13161671474
	M=0.004869491981
	W1=-0.159596235474
	C00=0.141189631152
	C10=0.0809701286525
	C01=-0.281528535557
	C11=0.15384112876
	C20=-0.178251207466
	C02=0.106959469314
	D0=0.0759196200467
	D1=-0.0217762490699
	R0=0.577350269
	AA=ALPHA^2
	BB=BETA^2
	A4=AA^2
	B4=BB^2
	ONMAA=1.-AA
	ONMBB=1.-BB
	X=ALPHA*(GSTAR+AA*(1.-GSTAR)+ONMAA*(BB*(G+(M-G)*AA  $
	  +ONMBB*(C00+C10*AA+C01*BB+C11*AA*BB+C20*A4+C02*B4)) $
	  +AA*(W1-ONMAA*(D0+D1*AA))))
	Y=BETA*(GSTAR+BB*(1.-GSTAR)+ONMBB*(AA*(G+(M-G)*BB  $
	  +ONMAA*(C00+C10*BB+C01*AA+C11*BB*AA+C20*B4+C02*A4))  $
	  +BB*(W1-ONMBB*(D0+D1*BB))))
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


