;  NOTE to future maintainers:
;   Make sure sub_sup_idl stays before translate_sub_super.  At least
;   for now, when IDL encounters a function and automatically compiles
;   it, it only compiles the functions in the file up to the named
;   function.  So even if sub_sup_idl was declared with
;   FORWARD_FUNCTION in translate_sub_super, it would not properly
;   compile. 
;
;+
; SPECIAL NOTE: 
;       The file translate_sub_super.pro contains two functions,
;       translate_sub_super, and sub_sup_idl.  The former is the
;       generic routine for processing TeX sub/superscripts, the
;       latter is used only by translate_sub_super and has no general
;       utility.  Hence it lives here.  You will see documentation for
;       translate_sub_super second if you use DOC_LIBRARY.
;-
;
;
;+
; NAME:
;       SUB_SUP_IDL
; PURPOSE:
;       Return the proper IDL font positioning command for TeX
;       sub/superscripts. 
; CATEGORY:
;       TeXtoIDL
; CALLING SEQUENCE:
;       fnt = sub_sup_idl( strn )
; INPUTS:
;       strn      -- Either '^' or '_', the TeX super/subscript       in
;                    characters
; KEYWORD PARAMETERS:
;       /FORCE_UD -- Set this to use !U/!D instead of !E/!I for
;                    sub/superscripts .
;       /HELP     -- Set to print useful message and exit.
; OUTPUTS:
;       fnt       -- Either '!U' or !E' for superscripts,             out
;                    or '!D' or '!I' for subscripts.
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       Used only by translate_sub_super.  Should be kept in same
;       file. 
; EXAMPLE:
; MODIFICATION HISTORY:
;       $Id: translate_sub_super.pro,v 1.5 2000/06/14 19:09:22 mcraig Exp $
;       $Log: translate_sub_super.pro,v $
;       Revision 1.5  2000/06/14 19:09:22  mcraig
;       Changed name of strtok str_token to avoid conflict in IDL 5.3.
;
;       Revision 1.4  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.3  1996/05/09 00:22:17  mcraig
;       Changed some function calls to reflect changes in those functions, moved
;       some code out of the main loop that didn't need to be there, added
;       documentation.
;
;       Revision 1.1  1996/01/31 18:47:37  mcraig
;       Initial revision
;
; RELEASE:
;       $Name: Rel_2_1_2 $
; COPYRIGHT:
;  Copyright (C) 1996 The Regents of the University of California, All
;  Rights Reserved.  Written by Matthew W. Craig.
;  See the file COPYRIGHT for restrictions on distrubting this code.
;  This code comes with absolutely NO warranty; see DISCLAIMER for details.
;-
FUNCTION Sub_sup_idl, token,  FORCE_UD = force_ud

; provide help if needed.
    IF (n_params() NE 1) OR keyword_set(Help) THEN BEGIN
        offset = '   '
        print, offset+'Return the proper IDL font positioning command for TeX'
        print, offset+'sub/superscripts. '
        print, offset+'fnt = sub_sup_idl( strn )'
        print, offset+'Inputs:'
        print, offset+offset+"strn      -- Either '^' or '_', the TeX super/subscript       in"
        print, offset+offset+'             characters'
        print, offset+'Keywords:'
        print, offset+offset+'/FORCE_UD -- Set this to use !U/!D instead of !E/!I for'
        print, offset+offset+'             sub/superscripts .'
        print, offset+offset+'/HELP     -- Set to print useful message and exit.'
        print, offset+'Outputs:'
        print, offset+offset+"fnt       -- Either '!U' or !E' for superscripts,             out"
        print, offset+offset+"             or '!D' or '!I' for subscripts."
        return, -1
    ENDIF 

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

;
;+
; NAME:
;       TRANSLATE_SUB_SUPER
; PURPOSE:
;       Translate TeX sub/superscripts to IDL sub/superscripts.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       new = translate_sub_super( old )
; INPUTS:
;       old       -- string to be translated from TeX to IDL.   in
; KEYWORD PARAMETERS:
;       /RECURSED -- set if this function is being called 
;                    recursively.                  
;       /HELP     -- Set to print useful message and exit.
; OUTPUTS:
;       new       -- string old converted from TeX to IDL       out
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       - For best results, when both a sub and superscript are used,
;         place the shorter of the two first (e.g. 'N^{a}_{bbbb}' is
;         better than 'N_{bbbb}^{a}').
;       - Single character sub/super scripts do not need to be
;         protected by braces.
;       - Sub/superscripts may be nested (e.g. 'N^{N_1^N}').
; EXAMPLE:
;       out = translate_sub_super( 'N^2_{big}' )
;       Then out='N!U2!N!Dbig!N' which looks like it should on the
;       display. 
; LIBRARY FUNCTIONS CALLED:
;       str_token      -- Text/string (mcraig)
;       sub_sup_idl -- contained in this file
; MODIFICATION HISTORY:
;       $Id: translate_sub_super.pro,v 1.5 2000/06/14 19:09:22 mcraig Exp $
;       $Log: translate_sub_super.pro,v $
;       Revision 1.5  2000/06/14 19:09:22  mcraig
;       Changed name of strtok str_token to avoid conflict in IDL 5.3.
;
;       Revision 1.4  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.3  1996/05/09 00:22:17  mcraig
;       Changed some function calls to reflect changes in those functions, moved
;       some code out of the main loop that didn't need to be there, added
;       documentation.
;
;       Revision 1.2  1996/02/08 18:54:20  mcraig
;       Changed default sub/superscript size to be !D/!U rather than !I/!E to
;       improve readability of plat annotations.
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
FUNCTION Translate_sub_super, InputString, $
                              RECURSED=recursed, $
                              HELP=Help

; Return to caller if error.
    On_error, 2

; Offer help if needed and/or desired
    IF (n_params() NE 1) OR keyword_set(help) THEN BEGIN
        offset = '   '
        print, offset+'Translate TeX sub/superscripts to IDL sub/superscripts.'
        print, offset+'new = translate_sub_super( old )'
        print, offset+'Inputs:'
        print, offset+offset+'old       -- string to be translated from TeX to IDL.   in'
        print, offset+'Keywords:'
        print, offset+offset+'/RECURSED -- set if this function is being called '
        print, offset+offset+'             recursively.                  '
        print, offset+offset+'/HELP     -- Set to print useful message and exit.'
        print, offset+'Outputs:'
        print, offset+offset+'new       -- string old converted from TeX to IDL       out'
        print, offset+'Notes:'
        print, offset+offset+'- For best results, when both a sub and superscript are used,'
        print, offset+offset+"  place the shorter of the two first (e.g. 'N^{a}_{bbbb}' is"
        print, offset+offset+"  better than 'N_{bbbb}^{a}')."
        print, offset+offset+'- Single character sub/super scripts do not need to be'
        print, offset+offset+'  protected by braces.'
        print, offset+offset+"- Sub/superscripts may be nested (e.g. 'N^{N_1^N}')."
        return, -1
    ENDIF 

;  To allow for nested scripts, use !E/!I instead of !U/!D for scripts
;  when called recursively.
    IF (NOT keyword_set(recursed)) THEN $
      ud = 1 $
    ELSE $
      ud = 0

;  Return to the normal level after making sub/superscript unless we
;  are recursed, which indicates we are processing a nested script.
    IF keyword_set(recursed) THEN fontRestore = '' ELSE fontRestore = '!N'

;  Initialize vars for processing scripts.
    SpcByte = (byte(' '))(0)    ;We need the BYTE value for a space below.
    strn = InputString
    pos = 0
    StorePos = ''
    RecallPos = ''
    OldToken =  ''
    LenLastScript = 0

; Grab next sub/superscript.  Token will be either '^' or '_'.
; RETURN if no scripts.
    Token = nexttok(strn,  '^_', pos = pos)
    if pos EQ -1 then return, InputString ;nothing to process.

    FntChange =  sub_sup_idl(Token)

; Our approach will be to grab the input string up to the next '^' or
; '_', then process the script we've found.
    NewString=str_token(strn,Token)

    WHILE  strlen(strn) GT  0 DO  BEGIN
;  Grab first char of sub/superscript.
        Script = strmid(strn, 0, 1)
        EndOfScript = 0         ;Position of end of this script.
        IF (Script EQ '{') THEN BEGIN   ; Scripts of more than 1 char.
            EndOfScript = matchdelim(strn)      
            Script = translate_sub_super(strmid(strn, 1, EndOfScript-1), $
                                         /recursed )
        ENDIF 
;     Grab rest of string _after_ the end of the script.        
        strn = strmid(strn, EndOfScript+1, $
                      strlen(strn)-EndOfScript-1)

;     Find the next script and prepare for processing it.
        FntChange = sub_sup_idl(Token, FORCE_UD = ud)
        OldToken = Token
        Token = nexttok(strn, '^_', POS = pos)

;     If the input is 'n^2_j', we want the '2' to be directly above
;     the 'j', rather than having the 'j' below and to the right of
;     the 2.  In other words, we want the first below, not the second.
;              2               2
;             N               N
;              J                J
;     To accomplish this, we need to save the position at which we
;     begin writing the 2 with a !S, and restore that position with a
;     !R after writing the 2.  The first section in the IF block below
;     handles the 'J' above, the thing after the first script.  We
;     don't care if there is another script following.  We also padd
;     the second script with spaces if it is shorter than the first to
;     make sure that whatever comes out after the scripts starts in
;     the proper place.  The worry is that without the spaces, the
;     input 'N^{looong}_{s} + 1' will end up with the + starting right
;     the 's' ends.
        IF (StorePos EQ '!S') THEN BEGIN
            StorePos = ''
            RecallPos = ''
;     calculate the difference in length between this script and the 
;     previous stacked one, removing font change commands (crudely by
;     guessing that the number of characters this takes is twice the
;     number of exclamation points).  The  + 1 below is a kludge.  I
;     don't know why, but I need one extra space.
            NumSpaces = LenLastScript - (strlen(script) - 2*strcnt(Script,'!'))
            NumSpaces = (NumSpaces + 1) > 0
            IF NumSpaces GT 0 THEN $
              Script = Script + string( replicate(SpcByte, NumSpaces) )
        ENDIF ELSE BEGIN
            IF (Token NE OldToken) AND (pos EQ 0) THEN BEGIN
;             The next script immediately folows this one.  Arrange to
;             save the position of the current script so that both begin
;             with the same horizontal position.
                StorePos = '!S'
                RecallPos = '!R'
                LenLastScript = strlen(Script) - 2*strcnt(Script,'!')
            ENDIF
        ENDELSE  

;  Continue building the IDL string, adding on our just processed script.
        NewString = NewString + StorePos + FntChange + Script + RecallPos $
          + FontRestore

        IF ( pos NE -1 ) THEN BEGIN     ; more left to process     
            NewString = NewString $
              + str_token(strn, Token)   
        ENDIF ELSE BEGIN                ; we are done
            NewString = NewString + strn
            strn = ''
        ENDELSE
    ENDWHILE 
    
    return, NewString
END




