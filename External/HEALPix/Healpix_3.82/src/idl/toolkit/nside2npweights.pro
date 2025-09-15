;===============================================
function nside2npweights, nside, error=error, help=help
;+
; npweights = nside2npweights(nside, error=)
;  returns the number of pixel-based weights (in compact non-redundant form) for a given Nside.
; Given the symmetries of the Healpix lay out, the number of non-redundant weights
; is  nf = ((3*Nside+1)*(Nside+1))/4  for any Nside
;        ~ Npix/16
; 
; Nside      1, 2,  4,  8,  16,  32,   64,   128,   256,    512,   1024,    2048,     4096,     8192
; Npweights  2, 5, 16, 56, 208, 800, 3136, 12416, 49408, 197120, 787456, 3147776, 12587008, 50339840
;
; version 1.0, 2017-11-09, adapted from n_fullweights in src/cxx/Heapix_cxx/weight_utils.cc
; version 1.0, 2018-05-22, much simplified
;-

routine = 'nside2npweights'
syntax = 'npweights = '+routine+'(nside, ERROR=, HELP=)'

if keyword_set(help) then begin
    error = 0
    doc_library,routine
    return,-1
endif

error = 1
if n_params() ne 1 then begin
    print,syntax
    return,-1
endif

; test validity of Nside
npix = nside2npix(nside, error=error)
if error then begin
    message,/info,'Unvalid Nside ('+strtrim(nside,2)+') provided'
    message,'Aborting.'
    return,-1
endif

; actual calculation
one = max(nside) gt 8192 ? 1LL : 1L
lnside = round(nside)
nf = (  ( 3*lnside + one ) * ( lnside + one )  ) / (4 * one)

error = 0
return, nf
end

