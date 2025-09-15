FUNCTION fij2pix, fij,res

; This function takes an n by 3 element vector containing the face, 
; column, and row number (the latter two within the face) of a pixel
; and converts it into an n-element pixel array for a given resolution.
;
; SPR 10476  Add Documentation  J.M. Gales  01/21/93
;
n = n_elements(fij)/3
	; get number of pixels

pixel = lonarr(n)
pixel_1 = lonarr(n)
	; generate output pixel and intermediate pixel arrays

i = fij(*,1)
j = fij(*,2)
	; get input column and row numbers

num_pix_face = long(4) ^ (res-1)
	; calculate the number of pixels in a face

pow_2 = lonarr(16)
FOR k=0,15 DO pow_2(k) = 2l^k
	; generate an array containing powers of 2 (bit masks)

FOR bit=0,res-2 DO BEGIN

	pixel_1 = pixel_1 OR ishft((pow_2(bit) AND i),bit)
	pixel_1 = pixel_1 OR ishft((pow_2(bit) AND j),bit+1)
	; if col bit set then set corresponding even bit in pixel_l
	; if row bit set then set corresponding odd bit in pixel_l

ENDFOR

pixel = fij(*,0)*num_pix_face + pixel_1
	; add face number offset

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


