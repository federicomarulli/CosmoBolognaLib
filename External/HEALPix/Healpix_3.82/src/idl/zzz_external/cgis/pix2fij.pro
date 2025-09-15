FUNCTION pix2fij, pixel,res

; This function takes an n-element pixel array and generates an n by 3
; element array containing the corresponding face, column, and row 
; number (the latter two within the face).
;
; SPR 10476  Add Documentation  J.M. Gales  01/21/93
;
n = n_elements(pixel)
	; get number of pixels
output = intarr(n,3)
	; generate output array

res1 = res - 1
num_pix_face = long(4) ^ res1
	; calculate number of pixels per face

face = pixel / num_pix_face
fpix = pixel - num_pix_face*face
	; calculate face number and pixel number within face

output(*,0) = face
	; store face number in output array

pow_2 = LONARR(16)
FOR i=0,15 DO pow_2(i) = 2l^i
	; generate power of 2 vector (bit mask)

i = intarr(n)
j = intarr(n)
	; generate column and row arrays

FOR bit=0,res-1 DO BEGIN

	i = i OR (pow_2(bit) * (1l AND fpix))
	; set col bit if corresponding even bit set in fpix
	; Note: LSB is numbered 0

	fpix = ishft(fpix,-1)

	j = j OR (pow_2(bit) * (1l AND fpix))
	; set row bit if corresponding odd bit set in fpix

	fpix = ishft(fpix,-1)

ENDFOR

output(*,1) = i
output(*,2) = j
	; store col and row in output arrays

RETURN, output
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


