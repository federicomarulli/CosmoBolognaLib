pro help_st, var, longname=longname, old=old
;+
; NAME:
;    HELP_ST
;
; PURPOSE:
;    provides HELP-like information on sub-structures of structure
;
; CATEGORY:
;
; CALLING SEQUENCE:
;    HELP_ST, Var
;
; INPUTS:
;    Var: IDL variable, of any kind
;
; SIDE EFFECTS:
;    prints out information
;
; RESTRICTIONS:
;    
;
; PROCEDURE:
;    If Var is an IDL structure: does a recursive HELP,/STRUCTURES on Var and
;     each of its substructure
;    If Var is NOT a structure: does the equivalent of HELP, Var
;
;
; EXAMPLE:
;    init_healpix  ; make sure !healpix is defined
;    help_st, !healpix
;
;    a=0
;    help_st, a+1
;
; MODIFICATION HISTORY:
;    2017-05-05, EH, IAP: list_structures
;    2017-05-10: became help_st in Healpix
;-

routine = 'help_st'
narg   = n_params()
syntax = 'SYNTAX:   '+routine+', var'
if (narg ne 1) then begin
    print,syntax
    return
endif
root='  .'
        
if size(var,/tname) eq 'STRUCT' then begin
    ; deal with structure
    if ~keyword_set(old) then begin
        ; only on first entry
        help, var, output=output,/structure
        print, output[0]
    endif

    ; recursive call to HELP
    names = tag_names(var)
    for i=0, n_tags(var)-1 do begin
        current = defined(longname) ? longname+'.'+names[i] : root+names[i]
        help_st,var.(i),longname=current, old=1
    endfor

endif else begin
    ; deal with non-structure
    if defined(var) then _myjunk_ = var
    help, _myjunk_, output=output
    if keyword_set(old) then begin
        ; describe non-structure tag
        arg = string(longname,form='(a-30)') ; left justified in 30 blanks
        print,strtrans(output,['_MYJUNK_','='],[arg,''])
    endif else begin
        ; treat non-structure entry argument
        help, /recall, output=history
        command = strupcase(strcompress(strmid(history[1],1),/remove_all) )
        arg     = strtrans(command, strupcase(routine)+',', '')
        if (strlen(arg) gt 30) then arg='<Expression>'
        print,strtrans(output,'_MYJUNK_',arg)
    endelse
endelse


return
end

