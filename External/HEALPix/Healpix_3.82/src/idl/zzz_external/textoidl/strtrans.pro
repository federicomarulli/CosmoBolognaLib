;
;+
; NAME:
;       STRTRANS
; PURPOSE:
;       Translate all occurences of one substring to another.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       new = strtrans(oldstr,from,to,ned)
; INPUTS:
;       oldstr -- string on which to operate.              in
;                 May be an array.
;       from   -- substrings to be translated. May be      in
;                 an array.
;       to     -- what strings in from should be           in
;                 translated to. May be an array.
; KEYWORD PARAMETERS:
;       /HELP  -- Set this to print useful message and 
;                 exit.
; OUTPUTS:
;       new    -- Translated string. Array if oldstr is    out          
;                 an array.
;       ned    -- number of substitutions performed in     out
;                 oldstr.  Array if oldstr is an array.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       - Any of old, from, and to can be arrays.  
;       - from and to must have the same number of elements.
; EXAMPLE:
;       inp='Many*bad!chars+in_here'
;       from=['*','!','+','_']
;       to  =[' ',' ',' ',' ']
;       out = strtrans(inp,from,to,ned)
;       Will produce out='Many bad chars in here', and set ned to 4.
; MODIFICATION HISTORY:
;       $Id: strtrans.pro,v 1.7 2004/06/15 17:25:54 mcraig Exp $
;       $Log: strtrans.pro,v $
;       Revision 1.7  2004/06/15 17:25:54  mcraig
;       Fixed bug in regular expression, changed array notation to square brackets
;
;       Revision 1.6  2004/01/11 01:49:00  mcraig
;       Changed format of one array to newer [] style to avoidf conflict with function name in astro library.
;
;       Revision 1.5  2001/11/23 21:14:35  mcraig
;       Added keywords /EXTRACT, /PRESERVE_NULL, /REGEX to call to
;       strsplit. This comes very close to reproducing the behavior of the
;       obsolete routine str_sep.
;
;       Revision 1.4  2001/11/21 19:13:23  mcraig
;       Changed str_sep to strsplit. The former is now considered obsolete by RSI.
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Sped up significantly by using str_sep to handle the translation.  No longer
;       relies on routines fromother user libraries.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: Rel_2_1_2 $
;
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
;
FUNCTION strtrans, InputString, from, to, ned,  $
                   HELP=Help

; Bomb out to caller if error.
    On_error, 2

; Offer help if we don't have at least InputString, from, and to, or
; if the user asks for it.
    IF (n_params() LT 3) OR keyword_set(help) THEN BEGIN
        offset = '   '
        print, offset+'Translate all occurences of one substring to another.'
        print, offset+'new = strtrans(oldstr,from,to,ned)'
        print, offset+'Inputs:'
        print, offset+offset+'oldstr -- string on which to operate.              in'
        print, offset+offset+'          May be an array.'
        print, offset+offset+'from   -- substrings to be translated. May be      in'
        print, offset+offset+'          an array.'
        print, offset+offset+'to     -- what strings in from should be           in'
        print, offset+offset+'          translated to. May be an array.'
        print, offset+'Outputs:'
        print, offset+offset+'new    -- Translated string. Array if oldstr is    out'
        print, offset+offset+'          an array.'
        print, offset+offset+'ned    -- number of substitutions performed in     out'
        print, offset+offset+'          oldstr.  Array if oldstr is an array.'
        print, offset+'Notes:'
        print, offset+offset+'- Any of old, from, and to can be arrays. ' 
        print, offset+offset+'- from and to must have the same number of elements.'
        return, -1
    ENDIF 
    
    strn = InputString

;  Check that From/To have same number of elements.  RETURN if they don't.
    NFrom = n_elements(from)
    NTo = n_elements(to)
    IF (NFrom EQ 0) OR (NTo EQ 0) THEN return, strn
    IF NFrom NE NTo THEN BEGIN
        print,'Error: Number of elements in from/to unequal'
        return,-1
    ENDIF

;  Make sure there are no null strings in From.  RETURN if there are.   
    FromLen = strlen(From)
    IF (total(FromLen EQ 0) GT 0) THEN BEGIN
        print, 'Error: elements of From must have nonzero length.'
        return, -1
    ENDIF 

    NStrings = n_elements(strn)
    ned = lonarr(NStrings)
    tmpned = 0L

; Say strn='a#b#c', from='#' and to='@'.  Then the approach here is to
; first split strn at all occurances of '#', then recombine the pieces
; with '@' inserted instead.  Do this for all elements of strn, and
; all elements of from.
    FOR i = 0L, NStrings-1 DO BEGIN
        ned[i] = 0L
        FOR j=0L, NFrom-1 DO BEGIN
            SepStr = strsplit(strn[i], from[j], $
                              /EXTRACT, /REGEX, /PRESERVE_NULL)
            NSubs = n_elements(SepStr) - 1
            strn[i] = SepStr[0]
            FOR k=1L, NSubs DO strn[i] = strn[i] + To[j] + SepStr[k]
            ned[i] =  ned[i] + NSubs
        ENDFOR 
    ENDFOR

    return, strn
END 
