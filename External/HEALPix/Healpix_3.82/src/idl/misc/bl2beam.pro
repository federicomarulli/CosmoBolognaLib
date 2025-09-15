function bl2beam, bl, theta, arcmin=arcmin, degrees=degrees, HELP=help, radians=radians
;+
; NAME:
;      bl2beam
;
; PURPOSE:
;   Computes a circular beam profile b(theta) from its transfer (or window) function b(l)
;
; CATEGORY:
;
; CALLING SEQUENCE:
;    beam = bl2beam(bl, theta [, ARCMIN=, DEGREES=, HELP=, RADIANS=])
;
; INPUTS:
;     bl: b(l) window function of beam
;     theta: angle
;
; KEYWORD PARAMETERS:
;     ARCMIN:  if set, theta is in arc-minutes
;     DEGREES: if set, theta is in Degrees
;     HELP:    if set, print out this doc header and leaves
;     RADIANS: if set, theta in is Radians
;
; OUTPUTS:
;     beam: (circular) beam profile value at radius theta
;
; RESTRICTIONS:
;     None
;
; PROCEDURE:
;     Since beam(vec) = \sum_lm b_lm Y_lm(vec)
;   then, for a circular beam
;   beam(theta) = \sum_l  B(l) P_l(theta) (2l+1)/(4 Pi): beam profile
;   where P_l is Legendre Polynomial
;   and B(l) = b_l0 / \sqrt{2l+1/(4 Pi)} : beam window (or transfer) function
;
;    BEAM2BL performs the inverse transform
;
; EXAMPLE:
;   bl = gaussbeam(15.d0, 4000, 1) ; window function of a 15arcmin-FWHM gaussian beam
;   theta = dindgen(3000)/100.    ; angle in arcmin
;   plot, theta,  bl2beam(bl, theta, /arcmin) ; its beam profile = a 15arcmin-FWHM gaussian
;
;
; MODIFICATION HISTORY:
;   EH, 2003-03-28: w2beam
;   EH, 2011-05-09: added Equations in doc header
;   EH, 2011-12-07: turned into bl2beam, 
;      added ARCMIN, DEGREES, HELP, RADIANS
;      added to Healpix release
;
;-

routine = 'bl2beam'
syntax = 'beam = '+routine+'(window, theta [, ARCMIN=, DEGREES=, HELP=, RADIANS=] )'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if n_params() lt 2 then begin
    print, syntax
    return,-1
endif

deg2rad = !DPi / 180.d0
do_am  = keyword_set(arcmin)
do_deg = keyword_set(degrees)
do_rad = keyword_set(radians)
if ((do_am + do_deg + do_rad) ne 1) then begin
    print, syntax
    message, 'Choose either ARCMIN, DEGREES or RADIANS'
endif
;;;if (((do_am + do_deg + do_rad) eq 0) then do_rad = 1
if do_deg then theta_rad = theta * deg2rad
if do_am  then theta_rad = theta * deg2rad / 60.d0
if do_rad then theta_rad = theta * 1.d0

lmax = n_elements(bl)-1
nx = n_elements(theta_rad)
x = cos(theta_rad)

; do P_l recursion and sum over l
p0 = replicate(1.d0,nx)
p1 = x

beam = bl[0]*p0 + bl[1]*p1*3.d0

for l=2L,lmax do begin
    fl = double(l)
    p2 = x * p1 * (2*fl-1.d0)/fl - p0 * (fl-1.d0)/fl
    p0 = p1
    p1 = p2
    beam += bl[l] * p2 * (2.d0*fl+1.d0)
endfor

beam /=  (4.d0*!dpi)

return,beam
end

