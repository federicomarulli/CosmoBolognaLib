;+
; handle_hash
; 
; wrapper for HASH related functions only supported in IDL >= 8.0 
; and in GDL and FL
;-
function handle_hash, h,create=create, tostructure=tostructure
if keyword_set(create) then begin
    return,hash()
endif
if keyword_set(tostructure) then begin
    return,h.tostruct()
endif
end
