pro uniq2nest, puniq, nside, pnest, help=help
;+
; NAME:
;      Uniq2Nest 
;
; PURPOSE:
;      Transforms Unique HEALPix pixel ID number(s) into Nested pixel number(s) and their Nside(s)
;
; CATEGORY:
;
; CALLING SEQUENCE:
;       Uniq2Nest, Puniq, Nside, Pnest, HELP=
;
; INPUTS:
;       Puniq: Unique HEALPix pixel ID number(s)
;            in [4*Nside^2, 16*Nside^2-1]
;
; KEYWORD PARAMETERS:
;       HELP: prints this documentation header and exits
;
; OUTPUTS:
;       Nside: HEALPix resolution parameter (power of 2)
;       Pnest: HEALPix NESTED pixel number (in [0, 12*Nside^2-1])
;
; SIDE EFFECTS:
;       None
;
; RESTRICTIONS:
;
; PROCEDURE:
;   Nside = 2 ^ Int(  Log_2(Puniq/4) / 2 )
;   Pnest = Puniq - 4 * Nside^2
;   tweaked to avoid round-off errors
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;   Aug 2005: original Unique Identifier equation
;   Sep 2015: first validated code
;
;-

routine = 'uniq2nest'
syntax = routine+', Puniq, Nside, Pnest, HELP='

if (keyword_set(help)) then begin
    doc_library,routine
    return
endif

if (n_params() ne 3) then begin
    print,syntax
    return
endif

if (min(puniq) lt 4) then begin
    print,syntax
    print,'Puniq must be 4 or larger'
    print, 'found range: ',min(puniq),max(puniq)
endif

tmp1= cheap_isqrt(puniq/4LL); use cheap_isqrt instead of sqrt to avoid round-off errors
tmp = alog(tmp1*1d0)/alog(2d0) ; must be done in double precision
nside = 2L ^ long( tmp   )
pnest = puniq - (4LL * nside) * nside


return
end

