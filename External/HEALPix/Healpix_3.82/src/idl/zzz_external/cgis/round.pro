function round, indata
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     ROUND rounds a set of data to nearest absolute integers.
;
;DESCRIPTION:  
;     Rounds data to the nearest absolute integer.
;
;CALLING SEQUENCE:  
;     RESULT = round(indata)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     OUTDATA       O   (Fix)      rounded value
;     INDATA        I              value to round
;
;WARNINGS:
;     None
;
;EXAMPLE:
;
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;     Written by BA Franz;  SPR 11132
;
;.TITLE
; Routine ROUND
;-
;
; Check on input parameters
;
on_error, 2

outdata = fix(indata)
frac = indata - outdata
iminus = where(frac le -0.5)
iplus  = where(frac ge  0.5)

if (iminus(0) ne -1) then outdata(iminus)=outdata(iminus)-1
if (iplus(0)  ne -1) then outdata(iplus)=outdata(iplus)+1

return,outdata
end
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


