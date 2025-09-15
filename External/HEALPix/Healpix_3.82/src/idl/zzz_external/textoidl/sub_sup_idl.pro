;
;+
; NAME:
;       SUB_SUP_IDL
; PURPOSE:
;       Return the proper IDL font positioning command for TeX
;       sub/superscripts. 
; CATEGORY:
; CALLING SEQUENCE:
;       fnt = sub_sup_idl( strn )
; INPUTS:
;       strn -- Either '^' or '_', the TeX super/subscript       in
;               characters
; KEYWORD PARAMETERS:
;       /FORCE_UD -- Set this to use !U/!D instead of !E/!I for
;                    sub/superscripts .
; OUTPUTS:
;       fnt -- Either '!U' or !E' for superscripts,              out
;              or '!D' or '!I' for subscripts.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
; EXAMPLE:
; LIBRARY FUNCTIONS CALLED:
;        
; MODIFICATION HISTORY:
;       $Id: sub_sup_idl.pro,v 1.1 1996/01/31 18:47:37 mcraig Exp $
;       $Log: sub_sup_idl.pro,v $
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: Rel_2_1_2 $
;-
FUNCTION Sub_sup_idl, token,  FORCE_UD = force_ud

IF keyword_set(force_ud) THEN BEGIN 
    IF (token EQ '^') THEN return, '!U' 
    IF (token EQ '_') THEN return, '!D'
    return, ''
ENDIF ELSE BEGIN
    IF (token EQ '^') THEN return, '!E' 
    IF (token EQ '_') THEN return, '!I'
    return, ''
ENDELSE

END


