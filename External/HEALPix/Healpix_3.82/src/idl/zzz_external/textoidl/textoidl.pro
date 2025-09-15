;
;+
; NAME:
;       TEXTOIDL
; PURPOSE:
;       Convert a valid TeX string to a valid IDL string for plot labels.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       new = textoidl(old)
; INPUTS:
;       old            -- TeX string to be converted.  Will not be     in
;                         modified.  old may be a string array.
; KEYWORD PARAMETERS:
;       FONT           -- Set to 0 to use hardware font, -1 to use 
;                         vector.  Note that the only hardware font 
;                         supported is PostScript.
;       /TEX_SEQUENCES -- return the available TeX sequences
;       /HELP          -- print out info on use of the function
;                         and exit.
; OUTPUTS:
;       new            -- IDL string corresponding to old.             out
; COMMON BLOCKS:
; SIDE EFFECTS:
; NOTES:
;       - Use the procedure SHOWTEX to get a list of the available TeX
;         control sequences.  
;       - The only hardware font for which translation is available is
;         PostScript. 
;       - The only device for which hardware font'
;         translation is available is PostScript.'
;       - The FONT keyword overrides the font selected'
;         by !p.font'
; EXAMPLE:
;       out = TeXtoIDL('\Gamma^2 + 5N_{ed}')
;       The string out may be used in XYOUTS or other IDL text
;       display routines.  It will be an uppercase Gamma, with an
;       exponent of 2, then a plus sign, then an N with the subscript
;       ed.
; MODIFICATION HISTORY:
;       $Id: textoidl.pro,v 1.7 2004/06/15 17:25:54 mcraig Exp $
;       $Log: textoidl.pro,v $
;       Revision 1.7  2004/06/15 17:25:54  mcraig
;       Fixed bug in regular expression, changed array notation to square brackets
;
;       Revision 1.6  2004/01/11 01:49:00  mcraig
;       Changed format of one array to newer [] style to avoidf conflict with function name in astro library.
;
;       Revision 1.5  2001/11/23 21:10:55  mcraig
;       Added backslash '\' to tex sequences in translation table to protect them during regexp search in strsplit.
;
;       Revision 1.4  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.3  1996/05/09 00:22:17  mcraig
;       Added error handling, cleaned up documentation.
;
;       Revision 1.2  1996/02/08 18:52:50  mcraig
;       Added ability to use hardware fonts for PostScript device.
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
FUNCTION Textoidl, InputString, $
                   FONT=fnt, $
                   HELP=hlp, $
                   TEX_SEQUENCES=tex_seq

;  Return to caller if there is an error.
    On_error, 2
;  We begin by deciding on the font.  PostScript = 0 means use vector.
    PostScript = 0
    IF n_elements(fnt) EQ 0 THEN BEGIN     ; get font from !p.font
        IF !p.font NE -1 THEN BEGIN        ; User wants hardware font.
            PostScript=1
        ENDIF
    ENDIF ELSE BEGIN                       ; get font from FONT keyword
        IF fnt NE -1 THEN PostScript = 1
    ENDELSE

;  Bomb out if user wants non-PostScript hardware font.
    IF (PostScript EQ 1) AND (!d.name NE 'PS') THEN BEGIN   
                                              ; Device isn't postscript 
                                              ; and user wants hardware
                                              ; font.  Not good.
        print,'Warning: No translation for device: ',!d.name
        return,InputString               
    ENDIF 
    
    IF keyword_set (tex_seq) THEN BEGIN
        table=textable()
        return,table[0,*]
    ENDIF 

    IF keyword_set(hlp) OR (n_params() EQ 0) THEN BEGIN
        print, '   Convert a TeX string to an IDL string'
        print, '   new = TeXtoIDL(old)'
        print, '     old = TeX string to translate.                 in'
        print, '     new = resulting IDL string.                    out'
        print, '   Keywords:'
        print, '      FONT       set to -1 to translate for vector fonts '
        print, '                 (DEFAULT) .  Set to 0 to translate for'
        print, '                 hardware font.'
        print, '      /TEX_SEQUENCES -- return the available TeX sequences'
        print, '      /HELP      print this message and exit.'
        print, '   NOTES:  '
        print, '      - Use SHOWTEX to obtain a list of the available'
        print, '        TeX control sequences.'
        print, '      - old may be a string array.  If so, new is too.'
        print, '      - The only device for which hardware font'
        print, '        translation is available is PostScript.'
        print, '      - The FONT keyword overrides the font selected'
        print, '        by !p.font'
        return, -1
    ENDIF
    
; PostScript has been set to 1 if PostScript fonts are desired.
    strn = InputString
    table = textable(POSTSCRIPT=PostScript)
    
;   Greek sub/superscripts need to be protected by putting braces
;   around them if they are unbraced.  This will have the result the
;   it will be difficult to use \ as a sub/superscript.  Get over it.

;   V2.11 Must include the '\' in from of translation table TeX
;   sequences to ensure that strsplit properly treats the '\' in the
;   TeX sequence. Since strsplit is doing a regexp replace, and '\' is
;   special in regexps, need to escape it.
    strn =  strtrans(strn, '\^'+'\'+table[0, *], '^{'+table[0, *]+'}')
    strn =  strtrans(strn, '\_'+'\'+table[0, *], '_{'+table[0, *]+'}')

;  First we translate Greek letters and the like.  This makes guessing
;  alignment of sub/superscripts easier, as all special characters will then
;  be one character long.

;   V2.11 Must include the '\' in from of translation table TeX
;   sequences to ensure that strsplit properly treats the '\' in the
;   TeX sequence. Since strsplit is doing a regexp replace, and '\' is
;   special in regexps, need to escape it.
    strn = strtrans(strn, '\'+table[0, *], table[1, *])


    FOR i = 0L, n_elements(strn)-1 DO BEGIN 
        strn[i] = translate_sub_super(strn[i]) ; Take care of sub/superscripts
    ENDFOR 

    return,strn
END 
