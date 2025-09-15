;
;+
; NAME:
;       SHOWTEX
; PURPOSE:
;       Display TeX sequence translation table on current graphics device.
; CATEGORY:
;       text/strings
; CALLING SEQUENCE:
;       showtex
; INPUTS:
; KEYWORD PARAMETERS:
;       /HELP -- print out info on use of the function
;                and exit.
;       FONT  -- Set to 0 to use hardware font, -1 to use vector.
;                Note that the only hardware font supported is
;                Postscript
; OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
;       Plot is created.
; NOTES:
;       Hardware fonts are supported only for device PS (PostScript)
; EXAMPLE:
; MODIFICATION HISTORY:
;       $Id: showtex.pro,v 1.4 2004/06/15 17:25:54 mcraig Exp $
;       $Log: showtex.pro,v $
;       Revision 1.4  2004/06/15 17:25:54  mcraig
;       Fixed bug in regular expression, changed array notation to square brackets
;
;       Revision 1.3  1996/06/14 20:00:27  mcraig
;       Updated Copyright info.
;
;       Revision 1.2  1996/05/09 00:22:17  mcraig
;       Added error handling and updated built in help.
;
;       Revision 1.1  1996/02/08 18:55:12  mcraig
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
PRO Showtex, FONT=fnt, HELP=help

; Return to caller on error.
    On_error, 2

; Print help if needed.
    IF keyword_set(help) THEN BEGIN
        print, '    Display TeX sequence translation table on current graphics device.'
        print, '    showtex'
        print, '    Keywords:'
        print, '       /HELP       print this message and return'
        print, '       FONT        set to 0 to use hardware fonts for current device,'
        print, '                   -1 to use vector fonts (DEFAULT)'
        print, '    NOTES:  - The only hardware font supported is PostScript.'
        print, '            - The FONT keyword overrides the font selected in !p.font'
        return
    ENDIF
    
;  We begin by deciding on the font.  PostScript = 0 means use vector.
    PostScript = 0
    PlotTitle = 'Vector Fonts'
    IF n_elements(fnt) EQ 0 THEN BEGIN ; get font from !p.font
        IF !P.font NE -1 THEN BEGIN ; User wants hardware font.
            PostScript = 1
            PlotTitle = 'PostScript Fonts'
        ENDIF 
    ENDIF ELSE BEGIN            ; get font from FONT keyword
        IF fnt NE -1 THEN BEGIN
            PostScript = 1
            PlotTitle = 'PostScript Fonts'
        ENDIF
    ENDELSE
    
;  Bomb out if user wants hardware font for non-PostScript device.
    IF (PostScript EQ 1) AND (strupcase(!D.name) NE 'PS') THEN BEGIN   
                                              ; Device isn't postscript 
                                              ; and user wants hardware
                                              ; font.  Not good.
        print, 'Warning: No translation for device: ', !D.name
        return
    ENDIF
   
; Set !P.font to value indicated by FONT keyword, saving surrent
; setting to reset at end.
    OldPFont = !p.font
    !p.font = PostScript - 1

    erase
    seq = textoidl(/tex)
    DisplayString = seq + '  ' + textoidl(seq)

    nseq = n_elements(seq)
    nrows = nseq/5 + 1          ; Five sequences per row.
    dx = .9/5.
    dy = .9/nrows
    y=.95
    xyouts,.5,y,PlotTitle,align=.5,/norm,size=2.5
    count=0
    FOR i = 1L, nrows DO BEGIN
        y= y - dy
        x = .1
        FOR j = 1, 5 DO BEGIN
            IF (count LT nseq ) THEN xyouts, x, y, DisplayString[count], align = .5, /norm
            count = count+1
            x = x + dx
        ENDFOR
    ENDFOR

; Restore old !P.font.
    !p.font = OldPFont
END
