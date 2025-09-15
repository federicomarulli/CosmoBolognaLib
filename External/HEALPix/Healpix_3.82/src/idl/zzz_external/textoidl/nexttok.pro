;
;+
; NAME:
;       NEXTTOK
; PURPOSE:
;       Find the next occurance of any of a set of characters in a
;       string and return the character which occurs next.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       tok = nexttok( strn, tokens )
; INPUTS:
;       strn   -- string to be searched for sub/superscripts    in
;       tokens -- string containing characters to be found.     in
; KEYWORD PARAMETERS:
;       POSITION -- Set to a named variable to get position     out
;                   of next token, or -1 if none found.
;       /HELP    -- Print useful message and exit.
; OUTPUTS:
;       tok    -- Contains the character among tokens which     out
;                 occurs next in strn, or null '' if none found.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
; EXAMPLE:
;       nexttok( 'x^2 + N_j^3', '^_', position=pos ) returns '^' and sets
;       pos to 1.
; MODIFICATION HISTORY:
;       $Id: nexttok.pro,v 1.4 2004/06/15 17:25:54 mcraig Exp $
;       $Log: nexttok.pro,v $
;       Revision 1.4  2004/06/15 17:25:54  mcraig
;       Fixed bug in regular expression, changed array notation to square brackets
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Generalized so that the next occurence of any of a set of characters will
;       be returned.
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
FUNCTION nexttok, strn, tokens, $
                  POSITION=position, $
                  HELP=Help

;  Return to caller on error.
    On_error, 2

;  Help those in need of it.
    IF (n_params() NE 2) OR keyword_set(Help) THEN BEGIN 
        offset = '   '
        print, offset+'Find the next occurance of any of a set of characters in a'
        print, offset+'string and return the character which occurs next.'
; CALLING SEQUENCE:
        print, offset+'tok = nexttok( strn, tokens )'
; INPUTS:
        print, offset+'Inputs:'
        print, offset+offset+'strn   -- string to be searched for sub/superscripts    in'
        print, offset+offset+'tokens -- string containing characters to be found.     in'
; KEYWORD PARAMETERS:
        print, offset+'Keywords:'
        print, offset+offset+'POSITION -- Set to a named variable to get position     out'
        print, offset+offset+'            of next token, or -1 if none found.'
        print, offset+offset+'/HELP    -- Print useful message and exit.'
; OUTPUTS:
        print, offset+'Outputs:'
        print, offset+offset+'tok   -- Contains the character among tokens which      out'
        print, offset+offset+"         occurs next in strn, or null '' if none found."
; EXAMPLE:
        print, offset+'Example:'
        print, offset+offset+"nexttok( 'x^2 + N_j^3', '^_', position=pos ) returns '^' and sets"
        print, offset+offset+'pos to 1.'
        return, ''
    ENDIF 

    TmpStr = byte(strn)
    TmpTok = byte(tokens)
    NumToks = n_elements(TmpTok) 

    MatchIdx = 0L
    Matches = 0L
    FOR j=0, NumToks-1 DO BEGIN 
        TmpMatch = where(TmpStr EQ TmpTok(j),  TmpCnt)
        IF (TmpCnt GT 0) THEN BEGIN
            MatchIdx = [MatchIdx, Replicate(j, TmpCnt)]
            Matches = [Matches, TmpMatch]
        ENDIF 
    ENDFOR 

    IF n_elements(MatchIdx) EQ 1 THEN BEGIN 
        Position = -1
        return, ''
    ENDIF 

    MatchIdx = MatchIdx[1:*]
    Matches = Matches[1:*]

    SortInd = sort(Matches)

    Position = Matches[SortInd[0]]

    Tok = string(TmpTok[MatchIdx[SortInd[0]]])
    
    return, Tok
END



