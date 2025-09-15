pro message_patch, text, level=level, informational=informational
;+
;  message_patch, text [,LEVEL=, INFORMATIONAL= ]
;
; wrapper around MESSAGE that emulates the LEVEL= keyword not supported by GDL
;
; Eric Hivon, 2012-03-02
;-

mylevel = keyword_set(level) ? level : 0
if mylevel lt 0 then mylevel -=1 ; make sure this routine does not appear in back tracing

if is_gdl() then begin
    help,call=trace
    nt = n_elements(trace)
    ii = (mylevel lt 0) ? abs(mylevel) : nt - mylevel - 1
    ii <= (nt - 1) 
    ii >= 0
    line = trace[ii]
    name = (strsplit(line, /extract))[0]
    message, name+': '+text, informational=informational, /noname
endif else begin

    message, text, level=mylevel, informational=informational
endelse


return
end

