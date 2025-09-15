function htmlcolor2rgb, hex, help=help
;+
; NAME:
;      htmlcolor2rgb
;
; PURPOSE:
;      converts HTML like hexadecimal color coding into RGB triplet
;      (see http://www.w3schools.com/html/html_colors.asp)
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;    rgb = htmlcolor2rgb(hex [, HELP=])
;
; INPUTS:
;    hex = 7 character string starting with '#' and containing a 6 digit hexadecimal number
;      (not case-sensitive: eg '#09aCdF')
;
; KEYWORD PARAMETERS:
;     HELP=  if set, this documentation header is printed and the routine exits
;
; OUTPUTS:
;    rgb = a 3-element vector, where each element is an integer in [0,255]
;
; PROCEDURE:
;
;
; EXAMPLE:
;      print,htmlcolor2rgb('#ffFffF'),htmlcolor2rgb('#102030')
;  will return
;    255 255 255
;     16  32  48
;
; MODIFICATION HISTORY:
;   2013-09-05: version 1.0, Eric Hivon
;
;-

routine = 'htmlcolor2rgb'
syntax  = 'rgb = '+routine+'(hex [, HELP=])'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if n_params() eq 0 then begin
    print,syntax
    return,-1
endif

hx = strtrim(hex[0],2) ; for future loop
if (size(hex,/tname) ne 'STRING' || strlen(hx) ne 7 || strmid(hx,0,1) ne '#') then begin
    print, syntax
    print, 'hex ='+hx
    message,'hex must be a 7 character string starting with #'
    return,-1
endif

hx = strupcase(strmid(hx, 1)) ; drop leading #

; converts HEX into long integer, adapted from strnumber.pro
On_IOerror, ioerror                 ;Go to ioerror if conversion error occurs
val = 0L
reads, hx, val, Format="(Z)"

; splits 24-bit into 3 bytes
r = ishft(val,-16) and 255b ; left most bits  -> Red
g = ishft(val, -8) and 255b ; middle bits     -> Green
b =       val      and 255b ; right most bits -> Blue

; returns 3-element vector
rgb = byte([r, g, b])
return, rgb

; handles IO errors
ioerror: message,'invalid input string'

end
