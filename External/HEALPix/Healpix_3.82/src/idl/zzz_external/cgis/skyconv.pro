function skyconv, in_uvec, inco=in_coord, outco=out_coord 
;+                                                                  
;  NAME:
;    skyconv
;
;  PURPOSE:                                   
;    Conversion routine between coordinate systems (ecl, gal, equ)
;
;  CALLING SEQUENCE:
;    out_uvec = skyconv(in_uvec,inco='in_coord',outco='out_coord')
;
;  INPUT:
;    in_uvec - array of input unit vectors
;    inco -  input coordinate system
;    outco - output coordinate system
;
;  OUTPUT:
;    out_uvec - array of output unit vectors
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;    J.M. Gales
;    Jan 92
;
; SPR 9741  Change PRINT statements to MESSAGE statements in
;           error situations.
; 2 JUN 1992
;
; SER 9997   Return input arrays directly for "no-op" transformations
; 17-SEP-92  JMG
;
;-
;

COMMON sky_conv_com, e2g,e2q,g2e,g2q,q2e,q2g,load_flag

n = n_elements(in_uvec) / 3

in_coord = strupcase(strmid(in_coord,0,1))
out_coord = strupcase(strmid(out_coord,0,1))


IF (in_coord eq ' ' AND out_coord eq ' ') THEN RETURN,in_uvec

case in_coord of

'E' :	begin

	case out_coord of

	'E' :	RETURN,in_uvec
	'G' :	out_uvec = in_uvec # e2g
	'Q' :	out_uvec = in_uvec # e2q

	else :	begin
		str = 'Unknown or improper output coordinate system: '
		str = str + '"' + out_coord + '"'
		MESSAGE,str,/cont
		RETURN,in_uvec
		end

	endcase

	end

'G' :	begin

	case out_coord of

	'G' :	RETURN,in_uvec
	'Q' :	out_uvec = in_uvec # g2q
	'E' :	out_uvec = in_uvec # g2e

	else :	begin
		str = 'Unknown or improper output coordinate system: '
		str = str + '"' + out_coord + '"'
		MESSAGE,str,/cont
		RETURN,in_uvec
		end

	endcase

	end



'Q' :	begin

	case out_coord of

	'Q' :	RETURN,in_uvec
	'E' :	out_uvec = in_uvec # q2e
	'G' :	out_uvec = in_uvec # q2g

	else :	begin
		str = 'Unknown or improper output coordinate system: '
		str = str + '"' + out_coord + '"'
		MESSAGE,str,/cont
		RETURN,in_uvec
		end

	endcase

	end


else :	begin 
	str = 'Unknown or improper input coordinate system: '
	str = str + '"' + in_coord + '"'
	MESSAGE,str,/cont
	RETURN,in_uvec
	end

endcase

exit:

RETURN, out_uvec
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


