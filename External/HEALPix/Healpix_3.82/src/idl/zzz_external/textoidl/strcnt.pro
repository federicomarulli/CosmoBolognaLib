;
;+
; NAME:
;       STRCNT
; PURPOSE:
;       Count number of occurrences of a substring in a string.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       num = strcnt(strn, substring, [pos])
; INPUTS:
;       string    -- The string in which to count occurences.     in
;       substring -- The substring to count occurrences of.       in
;       pos       -- the position at which to begin the search.   [in]
;                    If not supplied, start at beginning of
;                    string.
; KEYWORD PARAMETERS:
;       /HELP     -- Print useful message and return.
; OUTPUTS:
;       num       -- Number of occurances of substring in string. out
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       Overlapping occurances are not counted separately.  For
;       example, counting occurances of 'bb' in 'blah bbb' returns one
;       occurance. 
; EXAMPLE:
; MODIFICATION HISTORY:
;       $Id: strcnt.pro,v 1.3 1996/06/14 20:00:27 mcraig Exp $
;       $Log: strcnt.pro,v $
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Added fast processing using BYTE arrays if we are counting occurences of
;       a single character.  Added error handling.
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
FUNCTION Strcnt, strn, substrn, startpos, $
                 HELP=Help

; Return to caller if error.
    On_error, 2

; Help user, if needed.
    IF (n_params() LT 2) OR keyword_set(Help) THEN BEGIN
        offset = '   '
        print, offset+'Count number of occurrences of a substring in a string.'
        print, offset+'num = strcnt(strn, substring, [pos])'
        print, offset+'Inputs:'
        print,offset+offset+'string    -- The string in which to count occurences.     in'
        print,offset+offset+'substring -- The substring to count occurrences of.       in'
        print,offset+offset+'pos       -- the position at which to begin the search.   [in]'
        print,offset+offset+'             If not supplied, start at beginning of'
        print,offset+offset+'             string.'
        print, offset+'Keywords:'
        print,offset+offset+'/HELP     -- Print useful message and return.'
        print, offset+'Outputs:'
        print,offset+offset+'num       -- Number of occurances of substring in string. out'
        return, -1
    ENDIF 

    IF n_params() EQ 2 THEN startpos = 0

;  return if we weren't really given a substring to search for. . .
    IF strlen(substrn) EQ 0 THEN BEGIN 
        print, "Error: Can't count occurances of null string."
        return, -1
    ENDIF 

; . . .or if we were told to start at the end of the string.
    tmpstrn = strmid(strn, startpos, strlen(strn))
    IF strlen(tmpstrn) EQ 0 THEN return, 0

; If looking for occurences of single character, process using BYTE
; array.
    IF strlen(substrn) EQ 1 THEN BEGIN
        tmpstrn = byte(TmpStrn)
        count = n_elements(where(TmpStrn EQ (byte(substrn))(0))) 
    ENDIF ELSE BEGIN 
        count = 0L
        pos = rstrpos(tmpstrn, substrn)
        WHILE pos GE 0 DO BEGIN
            count = count + 1
            pos = rstrpos(tmpstrn, substrn, pos)
        ENDWHILE
    ENDELSE 

    return, count
END 








