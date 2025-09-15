pro nest2uniq, nside, pnest, puniq, help=help
;+
; NAME:
;      Nest2Uniq 
;
; PURPOSE:
;      Transforms Nside(s) and Nested pixel number(s) into  Unique HEALPix pixel ID number(s)
;
; CATEGORY:
;
; CALLING SEQUENCE:
;       Nest2Uniq, Nside, Pnest, Puniq, HELP=
;
; INPUTS:
;       Nside: HEALPix resolution parameter (power of 2)
;           Scalar or array of integer values
;
;       Pnest: HEALPix NESTED pixel number (in [0, 12*Nside^2-1])
;           if Nside is a scalar, Pnest can be scalar or array
;           if Nside is an array, Pnest must have the same number of elements
;
; KEYWORD PARAMETERS:
;       HELP: prints this documentation header and exits
;
; OUTPUTS:
;       Puniq: Unique HEALPix pixel ID number(s), 
;            in [4*Nside^2, 16*Nside^2-1]
;          scalar or array, same size as Pnest
;
; SIDE EFFECTS:
;       None
;
; RESTRICTIONS:
;
; PROCEDURE:
;   Puniq = Pnest + 4 * Nside^2
;
; EXAMPLE:
;
;
; MODIFICATION HISTORY:
;   Aug 2005: original Unique Identifier equation
;   Sep 2015: first validated code
;
;-

routine='nest2uniq'
syntax = routine+', nside, pnest, puniq, HELP='

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if (n_params() ne 3) then begin
    print,syntax
    return
endif

n1 = n_elements(nside)
n2 = n_elements(pnest)

if (n1 ne n2 and n1 ne 1) then begin
    print,syntax
    print,'Nside must be a scalar, or Nside and Pnest must have same length'
    return
endif

ns = long(nside)

; test that Nside is a power of 2
defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix
bad = where( (ns and (ns-1L)) ne 0, nbad)
if (min(nside) le 0 || max(nside) gt max(!healpix.nside) || nbad gt 0) then begin
    print,syntax
    print,'Invalid Nside value(s) (must be a power of 2  >0)'
endif
; npix = lon64arr(n1)
; for i=0L, n1-1 do begin
;     npix[i] = nside2npix(nside[i])
; endfor
; if (min(npix) lt 0) then begin
;     print,syntax
;     print,'Invalid Nside value(s) (must be a power of 2)'
;     return
; endif

offset = (4LL * ns) * ns
if (n1 eq 1) then offset = offset[0]
puniq = offset + pnest

return
end

