FUNCTION pix2uv, pixel,resolution
;+                                                                  
;  NAME:
;    pixel_2_uvec
;
;  PURPOSE:                                   
;        Routine returns unit vector pointing to center of pixel given pixel
;        number and resolution of the cube. 
;
;  CALLING SEQUENCE:
;    vector = pix2uv(pixel, resolution)
;
;  INPUT:
;    pixel -  Pixel number
;    resolution - Quad-cube resolution
;
;  OUTPUT:
;    vector - Unit vector to center of pixel
;
;  SUBROUTINES CALLED:
;    fwd_CUBE
;    XYAXIS
;
;  REVISION HISTORY
;    A.C. Raugh		9/90		Original FORTRAN version
;    J.M. Gales		1/92		IDL version
;-
;
;
;	  Variables:
;
;      Integer*4      FACE             ! Face number
;      Real           X,Y              ! Face coordinates
;      Real*4         SCALE              
;      Integer*4      PIXELS_PER_FACE  ! Number of pixels on a single face
;      Integer*4      FPIX             ! Pixel number within the face
;      Real           XI,ETA           ! 'Database' coordinates
;                                      !   overflow for large pixel numbers

	n_pixel = n_elements(pixel)

	bit_mask = long(1)

	scale = 2^(resolution-1) / 2.0

	out = pix2fij(pixel,resolution)

	x = (out(*,1) - scale + 0.5) / scale
	y = (out(*,2) - scale + 0.5) / scale

	fwdcube, x,y,xi,eta

	xyaxis, out(*,0),xi,eta,vector

	RETURN, vector
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


