;
;+
; NAME:
;       MATCHDELIM
; PURPOSE:
;        Match open/close delimiters in a string.
; CATEGORY:
;        text/strings
; CALLING SEQUENCE:
;        position = matchdelim( strn, [openpos])
; INPUTS:
;        strn        -- a string containing an open                 in
;                       delimiter (e.g. '{') in which you 
;                       want to find the matching closing  
;                       delimiter (e.g. '}')
; KEYWORD PARAMETERS:
;        OPEN_DELIM  -- A single character containing the opening   in
;                       delimiter (e.g. '(').  Default is '{'
;        CLOSE_DELIM -- A single character containing the closing   in
;                       delimiter (e.g. ')').  Default is '}'
; OUTPUTS:
;        position -- returns the position in strn of the            out
;                    closing delimiter, -1 if no closing found.
;        openpos  -- Set to a named variable to receive the         out
;                    position of the first opening delimiter.
;                    Optional.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;        - Any pair of (nonidentical) characters can be used as
;          delimiters. 
; EXAMPLE:
;        matchdelim('{one{two}}three') returns 9, the character just
;        before 'three'.  
; MODIFICATION HISTORY:
;       $Id: matchdelim.pro,v 1.3 1996/06/14 20:00:27 mcraig Exp $
;       $Log: matchdelim.pro,v $
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Removed restriction that open delim must be first char.  Added argument
;       to allow for return of position of open delim.
;
;       Revision 1.1  1996/01/31 18:41:06  mcraig
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
FUNCTION Matchdelim, InString, OpenPos, $
                     OPEN_DELIM=OpenDelim, $
                     CLOSE_DELIM=CloseDelim, $
                     HELP=Help

; Return to caller if error.
    On_error, 2

    IF (n_params() LT 1) OR keyword_set(Help) THEN BEGIN
        offset = '   '
        print, offset+'Match open/close delimiters in a string.'
        print, offset+'position = matchdelim( strn, [openpos])'
        print, offset+'Inputs:'
        print, offset+offset+'strn        -- a string containing an open                 in'
        print, offset+offset+"               delimiter (e.g. '{') in which you "
        print, offset+offset+'               want to find the matching closing  '
        print, offset+offset+"               delimiter (e.g. '}')"
        print, offset+'Keywords:'
        print, offset+offset+'OPEN_DELIM  -- A single character containing the opening   in'
        print, offset+offset+"               delimiter (e.g. '(').  Default is '{'"
        print, offset+offset+'CLOSE_DELIM -- A single character containing the closing   in'
        print, offset+offset+"               delimiter (e.g. ')').  Default is '}'"
        print, offset+'Outputs:'
        print, offset+offset+'position -- returns the position in strn of the            out'
        print, offset+offset+'            closing delimiter, -1 if no closing found.'
        print, offset+offset+'openpos  -- Set to a named variable to receive the         out'
        print, offset+offset+'            position of the first opening delimiter.'
        print, offset+offset+'            Optional.'
        print, offset+'Example:'
        print, offset+offset+"matchdelim('a{one{two}}three') returns 10, the character just"
        print, offset+offset+"  before 'three'.  "
        print, offset+offset+$
          "a=matchdelim('aaa[bbb(ccc)]ddd[eee]',f,OP='[',CL=']')"
        print, offset+offset+"  returns a=12 (just before ddd), f=3 "+$
          "(just before bbb)."  
        return, -1
    ENDIF 

; Set default delimiters.
    IF n_elements(OpenDelim) EQ 0 THEN OpenDelim =  '{'
    IF n_elements(CloseDelim) EQ 0 THEN CloseDelim =  '}'

; Make sure InString has more than 1 character.
    length = strlen(InString) 
    IF (length LE 1) THEN return,-1

; Return if no open delimiter
    OpenPos = strpos( InString, OpenDelim )
    IF (OpenPos EQ -1) THEN BEGIN 
        print, 'Error: No opening delimiter'
        return, -1
    ENDIF 
    
; Convert strings to array of integers to speed processing.
    OpenDelim = fix((byte(OpenDelim))(0))
    CloseDelim = fix((byte(CloseDelim))(0))
    TmpStr = fix(byte(strmid( InString, OpenPos, length)))
; Leave the -1* in here.  This forces conversion from BYTE to INTEGER,
; necessary because there are no negative BYTEs.
    TmpStr = (TmpStr EQ OpenDelim) $
              -1*(TmpStr EQ CloseDelim)
    length = n_elements(TmpStr) 

; Initialize count of number of delimiters.  We've found one, the
; first opener.
    BraceCnt = 1
    i=0
    WHILE (BraceCnt GT 0) AND (i LT length-1) DO BEGIN 
        i = i+1
        BraceCnt = BraceCnt + TmpStr(i)
    ENDWHILE 
    
    i = i + OpenPos
    IF (BraceCnt GT 0) THEN i = -1
    return, i
END

	



