;******************************************************************************
PRO loadsky
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
;*NAME:
;
;     	loadsky
;
;*PURPOSE:
;
;      	Calculate transformation matrices between ecliptic, galactic
;	and equatorial coordinate systems.
;
;*CALLING SEQUENCE:
;
;	load_sky_conv
;
;*PARAMETERS:
;
;    INPUT:
;
; 	None
;
;*COMMON BLOCKS:
;
;     	sky_conv_com
;
;*INTERACTIVE INPUT:
;
;	None.
;
;*SUBROUTINES CALLED:
;
;     	None.
;
;*FILES USED:
;
;     	None.
;
;*RESTRICTIONS:
;
;       None.
;
;*NOTES:
;
;       Transformation Matrices taken from USA program source code in
;           CSDR$SOURCE:[USALIB]
;
;       IDL performs matrix multiplication from Column to Row
;
;
;*MODIFICATION HISTORY:
;
;       Created 8/90 by Weihsueh Chiu
;	Modified for coord_conv 1/92 by J.M. Gales
;
;-------------------------------------------------------------------------------

COMMON sky_conv_com, e2g,e2q,g2e,g2q,q2e,q2g,load_flag

	eps = 23.452294 - 0.0130125 - 1.63889E-6 + 5.02778E-7
	eps = eps * 3.141592653589793 / 180.

	e2g =      [[-0.054882486, -0.993821033, -0.096476249],  $
                    [ 0.494116468, -0.110993846,  0.862281440],  $
                    [-0.867661702, -0.000346354,  0.497154957]]  

			; ecliptic to galactic


	e2q =      [[1.,     0.    ,      0.         ], $
	            [0., cos( eps ), -1. * sin( eps )], $
	            [0., sin( eps ),    cos( eps )   ]]

			; ecliptic to equatorial



	g2e = invert(e2g)
			; galactic to ecliptic

	g2q = g2e # e2q
			; galactic to equatorial



	q2e = invert(e2q)
			; equatorial to ecliptic

	q2g = q2e # e2g
			; equatorial to galactic



	load_flag = -1
			; set load flag

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


