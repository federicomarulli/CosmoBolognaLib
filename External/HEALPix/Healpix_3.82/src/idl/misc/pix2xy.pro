; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2013  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.sourceforge.net
;
; -----------------------------------------------------------------------------
pro pix2xy,pixel,x_out,y_out,res=res,data=data,raster=raster, $
	   face=face,sixpack=sixpack,bad_pixval=bad_pixval
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    PIX2XY creates a sky cube or face from a pixel/data list.
;
;DESCRIPTION:
;    This IDL procedure creates a raster image (sky cube or face)
;    given a pixel list and data array (provided, for example, by 
;    READ_SKYMAP).  The data array can be either a vector or 
;    two-dimensional array.  In the latter case, the data for each 
;    raster image can be stored in either the columns or rows.  The
;    procedure also returns the x and y raster coordinates of each 
;    pixel.  The rasterized images can be viewed with the IDL TV or 
;    TVSCL commands.  Only right oriented, ecliptic coordinate images
;    are built.
;
;CALLING SEQUENCE:
;    pix2xy,pixel,x_out,y_out,res=res,[data=data],[raster=raster], $
;           [bad_pixval=bad_pixval],[/face],[/sixpack]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    pixel       I     int/long arr        Input pixel list
;    [x_out]     O     int arr             x-coordinate array
;    [y_out]     O     int arr             y-coordinate array
;    res         I     int                 cube resolution
;    [data]      I     int/flt/dbl arr     input data array
;    [bad_pixval] I    flt                 Bad pixel value
;    [raster]    O     int/flt/dbl arr     rasterized output array
;    [face]      I     qualifier           Face output qualifier
;    [sixpack]   I     qualifier           Sixpack output qualifier
;
;WARNINGS:
;
;EXAMPLE: 
;
; To create x and y coordinate arrays from a FIRAS (res 6) pixel list
; without returning a rasterized image (since no data is provided):
;
; pix2xy,firas_pixel,x_out,y_out,res=6
;
;
; To create a rasterized DIRBE (res 9) sixpack from a pixel and data 
; list without returning the coordinate (x,y) arrays:
;
; pix2xy,dirbe_pixel,data=dirbe_data,res=9,raster=image,/sixpack
;
;
; To create a rasterized DIRBE face from a pixel and data list:
;
; pix2xy,dirbe_pixel,data=dirbe_data,res=9,raster=image,/face
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Checks data array to determine if vector array or 2D array
;    is provided and then whether each raster image is stored in the
;    columns or rows.
;    Calls RASTER (IDL procedure or c module if available) to 
;    perform rasterization
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: pix2fij, rastr
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 92
;
; SPR 9741  Change PRINT statements to MESSAGE statements in
;           error situations.
;
; SPR 9749  Fix output array allocation statement for single face
;           complex data.
;
; SPR 9833  Allows user to designate sixpacked output using command
;           line qualifier.
;
; 24 JUL 1992
;
; SER 9982     Implement CRASTER (C module rasterizer)
; 14-SEP-1992  JMG
;
; SPR 10059    Too large pixel number (for given resolution)
;              causes program to exit IDL.  Check max pix #.
; 29-SEP-1992  JMG
;
; SPR 10212    Call RASTER to generate rasterized image.
;              RASTER is implemented as both IDL routine and C module.
; 13-NOV-1992  JMG
;
; SPR 10212    Change call from 'RASTER' to 'RASTR'
; 02-MAR-1993  JMG
;
; SPR 10729    REFORM pixel list if column vector
;              Add BAD_PIXVAL keyword
; 22-MAR-1993  JMG
;
; 12-oct-2001  EH : replace switch by switch1 to remain compatible
;              with idl 5.4
;
;.TITLE
;Routine PIX2XY
;-
;

pix_sz = SIZE(pixel)
IF (pix_sz(0) NE 1) THEN BEGIN
	pixel = REFORM(pixel)
	pix_sz = SIZE(pixel)
ENDIF

switch1 = 0

IF (N_ELEMENTS(res) EQ 0) THEN BEGIN
	MESSAGE,'Map resolution must be specified',/CONT
	GOTO, exit
ENDIF

IF (MAX(pixel) GT 6*LONG(4)^(res-1)) THEN BEGIN
	MESSAGE,'Maximum pixel number too large for resolution',/CONT
	GOTO, exit
ENDIF


; Determine size and "orientation" of data array
; ----------------------------------------------
IF (n_elements(data) NE 0) THEN BEGIN

	data_info = SIZE(data)


	CASE data_info(0) OF

	1: 	BEGIN

		IF (pix_sz(1) NE data_info(1)) THEN BEGIN
			str = 'Pixel and Data Arrays of Incompatible Size'
			MESSAGE,str,/cont
			GOTO, exit
		ENDIF

		END


	2:	BEGIN

		IF (pix_sz(1) EQ data_info(2)) THEN switch1 = 1

		IF ((pix_sz(1) NE data_info(1)) AND $
		    (pix_sz(1) NE data_info(2))) THEN BEGIN
			str = 'Pixel and Data Arrays of Incompatible Size'
			MESSAGE,str,/cont
			GOTO, exit
		ENDIF

		END

	ELSE:	BEGIN

		MESSAGE,'Data array must be vector or 2D array',/CONT
		GOTO, exit

		END

	ENDCASE

ENDIF

IF (N_ELEMENTS(data) EQ 0) THEN data = [-1]
IF (KEYWORD_SET(face) EQ 0) THEN face = -1
IF (NOT KEYWORD_SET(sixpack)) THEN sixpack = 0
IF (KEYWORD_SET(bad_pixval) EQ 0) THEN bad_pixval = 0.0
		; set up flag values for RASTR

pix_long = 0
IF (pix_sz(pix_sz(0)+1) EQ 2) THEN BEGIN
	pixel = LONG(PIXEL)
	pix_long = 1
ENDIF
		; if pix is int array then temporary transform
		; to long int for compatability with c module
		; RASTR

IF (switch1 EQ 1) THEN data = TRANSPOSE(data)

rastr,pixel,FIX(res),face,sixpack,data,raster,x_out,y_out,FLOAT(bad_pixval)
		; call rasterization routine

IF (pix_long EQ 1) THEN pixel = FIX(pixel)
		; switch back to short int if necessary

IF (switch1 EQ 1) THEN data = TRANSPOSE(data)


exit:

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


