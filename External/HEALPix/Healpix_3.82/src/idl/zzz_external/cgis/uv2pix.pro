FUNCTION uv2pix, vector,resolution
;+                                                                  
;  NAME:
;    uv2pix
;
;  PURPOSE:                                   
;        Routine returns pixel number given unit vector pointing to center 
;        of pixel resolution of the cube. 
;
;  CALLING SEQUENCE:
;    pixel = uv2pix(vector, resolution)
;
;  INPUT:
;    vector - Unit vector to center of pixel 
;    resolution - Quad-cube resolution
;
;  OUTPUT:
;    pixel -  Pixel number 
;
;  SUBROUTINES CALLED:
;    axisxy
;
;  REVISION HISTORY
;    A.C. Raugh		9/90		Original FORTRAN version
;    J.M. Gales		1/92		IDL version
;-
;

	two = 2l
	n_vec = n_elements(vector)/3

	res1 = resolution - 1
	num_pix_side = long(2) ^ res1
	num_pix_face = num_pix_side ^ 2

		axisxy, vector,face,x,y
			; face - face # (0-5)
			; x,y - coordinates on face (0-1)

		i = fix(x * num_pix_side) < (2^res1 - 1)
		j = fix(y * num_pix_side) < (2^res1 - 1)

		pixel = fij2pix([[face],[i],[j]],resolution)

RETURN, pixel
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


