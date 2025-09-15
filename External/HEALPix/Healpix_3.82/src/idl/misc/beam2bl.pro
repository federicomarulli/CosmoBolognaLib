function myintegrate, f, x

res = int_tabulated(x, f, /double)

return, res
end

function beam2bl, beam, theta, lmax, arcmin=arcmin, degrees=degrees, HELP=help, radians=radians
;+
; NAME:
;   beam2bl
;
; PURPOSE:
;   computes a transfer (or window) function from its beam profile b(theta)
;
; CATEGORY:
;
; CALLING SEQUENCE:
;   bl = beam2bl(beam, theta, lmax [, ARCMIN=, DEGREES=, HELP=, RADIANS=]) 
;
; INPUTS:
;   beam:  Vector:              circular beam profile in theta
;   theta: Vector of same size: radius
;   lmax: Scalar integer, maximum multipole of bl
;
; KEYWORD PARAMETERS:
;     ARCMIN:  if set, theta is in arc-minutes
;     DEGREES: if set, theta is in Degrees
;     HELP:    if set, print out this doc header and leaves
;     RADIANS: if set, theta in is Radians
;
;
; OUTPUTS:
;   bl: beam window function B(l)
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;   Since b_lm = \int vec B(vec) Y_lm^*(vec)
; then, for a circular beam
; b_l0 = \int  B(theta) P_l(theta) sin(theta) dtheta 2Pi \sqrt{(2l+1)/(4 Pi)}
; where P_l is Legendre Polynomial
; and B(l) = b_l0 / \sqrt{2l+1/(4 Pi)} : beam window (or transfer) function
;          = \int  B(theta) P_l(theta) sin(theta) dtheta 2Pi
;
; calls INT_TABULATED to perform integration
;
; BL2BEAM  performs the inverse transform
;
; EXAMPLE:
;   bl = gaussbeam(15.d0, 4000, 1) ; window function of a 15arcmin-FWHM gaussian beam
;   theta = dindgen(4000)/100.     ; angle in arcmin
;   beam = bl2beam(bl, theta, /arcmin) ; its beam profile = a 15arcmin-FWHM gaussian
;   bl1 = beam2bl(beam, theta, 4000, /arcmin) ; back to window function
;   plot, bl1-bl ; absolute discrepancy between starting B(l) and the one 
;                ; reconstructed from real space beam profile
;                ; 
;
; MODIFICATION HISTORY:
;
; EH, 2011-09-09: beam2w
; EH, 2011-12-08: turned into beam2bl, 
;      added ARCMIN, DEGREES, HELP, RADIANS
;      added to Healpix release
;
;-

routine = 'beam2bl'
syntax = 'bl = '+routine+'(beam, theta, lmax [, ARCMIN=, DEGREES=, HELP=, RADIANS=] )'

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if n_params() lt 3 then begin
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

nx = n_elements(theta_rad)
nb = n_elements(beam)
if (nb ne nx) then begin
    message, 'Beam and Theta must have same size'
endif
x = cos(theta_rad)
st = sin(theta_rad)
window = dblarr(lmax+1)

; do P_l recursion
p0 = replicate(1.d0,nx)
p1 = x

; window[0] = myintegrate(beam*p0,x)
; window[1] = myintegrate(beam*p1,x)
window[0] = myintegrate(beam*p0*st,theta_rad)
window[1] = myintegrate(beam*p1*st,theta_rad)

for l=2L,lmax do begin
    fl = double(l)
    p2 = x * p1 * (2*fl-1.d0)/fl - p0 * (fl-1.d0)/fl
    ; window[l] = myintegrate(beam*p2,x)
    window[l] = myintegrate(beam*p2*st,theta_rad)
    p0 = p1
    p1 = p2
endfor

window *= (2.d0*!dpi)

return,window
end

