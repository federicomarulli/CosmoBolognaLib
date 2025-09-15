pro degrade_ring, nside_in, ir_in, ip_in, nside_out, ir_out, ip_out, help=help
;+
; NAME:
;          Degrade_ring
;
;
; PURPOSE:
;          Downgrade RING based pixel indices
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;          degrade_ring, Nside_in, Ir_in, Ip_in, Nside_out, Ir_out, Ip_out [, HELP=]
;
;
; INPUTS:
;         Nside_in: (scalar integer)
;                input Healpix resolution parameter (a power of 2)
;
;         Ir_in: (scalar or vector integer)
;                index of ring in wich the pixel lies [1,4*Nside_in-1]
;
;         Ip_in: (integer, same size as Ir_in)
;                index of pixel within the ring in [0,4*Nside_in-1]
;
;         Nside_out: (scalar integer)
;                output Healpix resolution parameter (a power of 2)
;                Must be smaller than (or equal to) Nside_in
;
;
; KEYWORD PARAMETERS:
;         Help: if set, this header is printed out and the routine exits
;
;
; OUTPUTS:
;         Ir_out: (integer, same size as Ir_in)
;                index of ring in wich the pixel lies [1,4*Nside_out-1]
;
;         Ip_out: (integer, same size as Ir_in)
;                index of pixel within the ring in [0,4*Nside_out-1]
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;          2011-04-01: EH/IAP, version 1.0
;
;-


routine = 'degrade_ring'
syntax = routine+', Nside_in, Ir_in, Ip_in, Nside_out, Ir_out, Ip_out [,HELP=]'
syntax2 = 'Ir_in and Ip_in arrays must have same size. Current sizes: '
syntax3 = 'Nside_in and Nside_out must be scalar integer'
syntax4 = 'For extended help, type          '+routine+', /help'

if keyword_set(help) then begin
    doc_library,routine
    return
end

; number of arguments
if n_params() lt 6 then begin
    print,syntax
    print,syntax4
    return
endif

; consistency of input arguments
np1  = n_elements(ir_in)
np2  = n_elements(ip_in)
if (np1 ne np2) then begin
    print,syntax
    print,syntax2+ string([np1, np2],form='(i10,1x,i10)')
    print,syntax4
    return
endif
if (n_elements(nside_in) gt 1 ||  n_elements(nside_out) gt 2) then begin
    print,syntax
    print,syntax3
    print,syntax4
    return
endif
junk = nside2npix(nside_in,  error=err1)
junk = nside2npix(nside_out, error=err2)
if (err1 || err2) then begin
    print,syntax
    print,'Nside_in or Nside_out are not power of 2: '+string(nside_in, nside_out,form='(i10,1x,i10)')
    return
endif

;    ---------------------------------------

lnside_in  = long(nside_in[0])
lnside_out = long(nside_out[0])
ratio = lnside_in / lnside_out

if (ratio eq 1) then begin
    ir_out = ir_in
    ip_out = ip_in
    return
endif

one = (nside_in le 8192) ? 1L : 1LL

if (ratio lt 1) then begin
    message,'Nside_out must be smaller than Nside_in'
    return
endif

nli4 = 4L * lnside_in
nlo4 = 4L * lnside_out

; put all pixels on a single virtual Healpix face at input resolution
; -------------------------------------------------------------------
; ir_in is unchanged, ip_in is changed in polar caps
ip = ip_in * one ; general case

; deal with polar caps
k_n = where(ir_in lt   lnside_in, nk_n) ; North cap
if (nk_n gt 0) then begin
    nir  = ir_in[k_n]
    iphi = ip[k_n]
    ip[k_n] = (iphi / nir) * lnside_in + (iphi mod nir) + (lnside_in - nir)/2
    nir = 0 & iphi = 0
endif

k_s = where(ir_in gt 3L*lnside_in, nk_s) ; South cap
if (nk_s gt 0) then begin
    nir  = nli4-ir_in[k_s]
    iphi = ip[k_s]
    ip[k_s] = (iphi / nir) * lnside_in + (iphi mod nir) + (lnside_in - nir)/2
    nir = 0 & iphi = 0
endif

shift = ~(ir_in and 1) or (ir_in lt lnside_in) or (ir_in gt 3*lnside_in) ; 1 in caps, (0,1) alternating in equat. region
if (lnside_in eq 1) then shift = ~shift

; compute (X=North-East, Y=South-East) at output resolution
; ---------------------------------------------------------
scale = 2 * ratio
up  = nli4      & dn  = up / scale ; offset used to avoid dealing with Euclidian division of <0 numbers
up2 = lnside_in & dn2 = up2 / scale ; offset to give correct results for Nside_out = 1
x = (2L * ip - ir_in + shift + up + up2) / scale  - dn - dn2
y = (2L * ip + ir_in + shift      + up2) / scale       - dn2

; ring based indices at new resolution
; ------------------------------------
; general case
ir_out =  y - x
shift_out = (ir_out and 1) and (ir_out ge lnside_out) and (ir_out le 3*lnside_out) ; 0 in caps, (1,0) alternating in equat. region
if (lnside_out eq 1) then shift_out = ~shift_out
ip_out = (y + x +shift_out) / 2
ip_out = ip_out mod nlo4

; deal with polar caps
if (nk_n gt 0) then begin
    nir  = ir_out[k_n]
    face = ip_out[k_n] /lnside_out
    ip_out[k_n] += face*(nir - lnside_out) - (lnside_out - nir)/2
endif
if (nk_s gt 0) then begin
    nir  = nlo4 - ir_out[k_s]
    face = ip_out[k_s] /lnside_out
    ip_out[k_s] += face*(nir - lnside_out) - (lnside_out - nir)/2
endif

return
end

