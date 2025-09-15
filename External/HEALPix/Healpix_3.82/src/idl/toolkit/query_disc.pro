; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2013  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.sourceforge.net
;


;===============================================================
pro pixels_per_ring, nside, npr, kshift, npc
;-
; nside: input
; npr: number of pixel per ring, vector of length 4*nside+1
; kshift: shift in phi (0 or 1), vector of length 4*nside+1
; npc: number of pixels between North pole and current ring, vector of length 4*nside+1
;   with ring=0 -> North Pole, ring=4*nside -> South Pole
;-
nl3 = 3L*nside
nl4 = 4L*nside
iring  = (nside le 8192) ? lindgen(nl4+1) : l64indgen(nl4+1)
npr = iring < (nl4 - iring)
npr <= nside
npr *= 4
kshift = ~(iring and 1)
if (nside eq 1) then kshift = ~kshift
kshift[0:nside] = 1
kshift[nl3:nl4] = 1
if (arg_present(npc)) then begin
    npc = total(npr, /cumul, /integer)
endif

return
end

;===============================================================
pro pixels_on_edge, nside_in, irlist, phi0, dphi, ringphi, ngr

pixels_per_ring, nside_in, npr, kshift ;, npc
nl4 = 4L*nside_in
badvalue = -9999 ; NEGATIVE value given to phi index on empty rings


count = histogram(dphi/(!dpi-1.0d-8)>(-1)<1, rev=rev, min=-1, max=1)
n_empty = count[0] & n_normal = count[1] & n_full = count[2] & n_non_empty = n_normal+n_full
if (n_empty gt 0)  then empty  = rev[rev[0]:rev[1]-1]
if (n_normal gt 0) then normal = rev[rev[1]:rev[2]-1]
if (n_full gt 0)   then full   = rev[rev[2]:rev[3]-1]
if (n_non_empty gt 0) then non_empty = rev[rev[1]:rev[3]-1]

nrings = n_elements(irlist)
ip_low = lonarr(nrings) ; pixel range on each ring
ip_hi = ip_low

if (n_empty gt 0) then begin
    ip_low[empty] = badvalue
    ip_hi[empty]  = badvalue
endif

twopi = 2.d0*!DPI
if (n_full gt 0) then begin
    nr = npr[irlist[full]]
    ip_low[full] = 0
    ip_hi[full] = nr-1
endif

if (n_normal gt 0) then begin
    ir = irlist[normal]
    nr = npr[ir]
    shift = (kshift[ir]) * .5d0
    phi_low = phi0[normal] - dphi[normal] ; in [-Pi,2Pi]
    phi_hi  = phi0[normal] + dphi[normal] ; in [0,3Pi]
    iphi_low = CEIL  (nr * phi_low / twopi - shift)
    iphi_hi  = FLOOR (nr * phi_hi  / twopi - shift)
    miss = where(iphi_low gt iphi_hi, nmiss) ; phi range is >0 but does not contain any pixel center
    iphi_low = modulo (iphi_low, nr)     ; in {0,nr-1}
    iphi_hi  = modulo (iphi_hi,  nr)
    if (nmiss gt 0) then begin
        iphi_low[miss] = badvalue
        iphi_hi[miss]  = badvalue
    endif
    ip_low[normal] = iphi_low
    ip_hi [normal] = iphi_hi
endif 


good = where(ip_low ge 0, ngr)
if (ngr eq 0) then begin
    ringphi = [-1L, -1L, -1L]
endif else begin
    ringphi = lonarr(ngr, 3)
    ringphi[0,0] = irlist[good]
    ringphi[0,1] = ip_low[good]
    ringphi[0,2] = ip_hi [good]
endelse


; if (n_non_empty eq 0) then begin
;     ngr = 0
;     ringphi = [-1L,-1L,-1L]
; endif else begin
;     good = non_empty[sort(non_empty)] ; sort non-empty rings in North to South order
;     ngr  = n_elements(good)
;     ringphi = lonarr(ngr, 3)

;     ringphi[0,0] = irlist[good]
;     ringphi[0,1] = ip_low[good]
;     ringphi[0,2] = ip_hi [good]
; endelse


return

end


;===============================================================

pro discedge2fulldisc, nside, ringphi, ngr, pir, npir

npir = 0
if (ngr eq 0) then begin
    pir[0] = -1
    return
endif

pixels_per_ring, nside, npr, kshift, npc
nl4 = 4L*nside
ramp = lindgen(nl4)

nj = n_elements(ringphi[*,0])
for j=0L, nj-1 do begin
    ir = ringphi[j,0]             ; current ring, in [1, nl4-1]
    nr = npr[ir]                ; number of pixels available on this ring
    my_low = ringphi[j,1]
    if (my_low ge 0) then begin
        my_hi  = ringphi[j,2]
        np = my_hi - my_low     ; in [-nr+1, nr-1]
        np = (np ge 0) ? np + 1 : nr + np + 1 ; deal with periodic BC
        np <= nr
        ip = modulo(my_low + ramp[0:np-1], nr)
        pir[npir] = npc[ir-1] + ip ; fill final list
        npir += np
    endif
endfor

if npir eq 0 then pir[0] = -1


return
end

;===============================================================
pro find_pixel_bounds, nside, nsboost, iring, iphi, phiw, phie, np=np

common pixbounds, ns_old, nsb_old, npr, kshift, f
if undefined(ns_old)  then ns_old  = -1
if undefined(nsb_old) then nsb_old = -1
if (nside ne ns_old || nsboost ne nsb_old) then begin
    ;;;;;print,'Init variables',ns_old,nsb_old
    pixels_per_ring, nside, npr, kshift
    f = (dindgen(2*nsboost+1) - nsboost)/nsboost ; in [-ns, ns]/(ns)
    ns_old = nside
    nsb_old = nsboost
endif

halfpi = !dpi * .5d0
nq = npr[iring]/4 ; number of pixels on current ring in [0,Pi/2] (quadrant)
transition = (iring eq nside || iring eq nside*3)
np = npr[iring] ; number of pixel on  current ring

if (nq eq nside || transition) then begin ; equatorial region (and transition rings)
    k0 = kshift[iring]*.5d0

    ff = (1.d0-abs(f))*0.5d0    ; triangle of height 1/2
    c0 = halfpi * (iphi + k0) / nq
    f1 = halfpi * ff / nq
    phiw = c0 - f1
    phie = c0 + f1
    if (transition) then begin ; store for future use
        phiw_t = phiw
        phie_t = phie
    endif
endif
if (nq lt nside || transition) then begin ; polar regions and transition rings
    ip = iphi mod nq ; in [0,nq-1]
    quad = iphi / nq ; quadrant in [0,3]
    curve = (iring le nside*2) ? halfpi / (nq + f) : halfpi / (nq - f); swap sign for South pole
    phiw1 = curve * ip     < halfpi
    phie1 = curve *(ip+1)  < halfpi
    phiw2 = curve *(nq-ip-1) < halfpi
    phie2 = curve *(nq-ip)  < halfpi
    
    phiw = phiw1 > (halfpi - phie2)
    phie = phie1 < (halfpi - phiw2)
    phiw += quad * halfpi
    phie += quad * halfpi
endif

if (transition) then begin
;     if (iring eq nside) then begin ; transition in N hemisphere
;         phiw = [phiw[0:nsboost],phiw_t[nsboost+1:2*nsboost]]
;         phie = [phie[0:nsboost],phie_t[nsboost+1:2*nsboost]]
;     endif else begin ; transition in S hemisphere
;         phiw = [phiw_t[0:nsboost],phiw[nsboost+1:2*nsboost]]
;         phie = [phie_t[0:nsboost],phie[nsboost+1:2*nsboost]]
;     endelse
    if (iring eq nside) then begin ; transition in N hemisphere
        phiw[nsboost+1] = phiw_t[nsboost+1:2*nsboost]
        phie[nsboost+1] = phie_t[nsboost+1:2*nsboost]
    endif else begin ; transition in S hemisphere
        phiw[0] = phiw_t[0:nsboost]
        phie[0] = phie_t[0:nsboost]
    endelse
endif

return
end
;===============================================================
pro check_edge_pixels, nside, wnside, izlist, phi0list, dphilist, ringphi, ngr

nsboost = wnside/nside
if (nsboost le 1) then return

; print, ringphi
;;;sub0 = lindgen(2*nsboost-1) - (nsboost-1) - izlist[0]
sub0 = lindgen(2*nsboost+1) - nsboost - izlist[0]
for i=0, ngr-1 do begin ; loop on low-res rings
    subrings = ringphi[i,0]* nsboost + sub0
    for k=-1,1,2 do begin ; West and East side of disc
        nshrink = 0
        kk = (k+3)/2 ; 1 or 2
rep:
;        if (ringphi[i,kk] ge 0 && ringphi[i,1] le ringphi[i,2]) then begin
        if (ringphi[i,kk] ge 0) then begin ; corrected 2020-08-28
            find_pixel_bounds, nside, nsboost, ringphi[i,0], ringphi[i,kk], phiw, phie,np=np
            phic = (phie+phiw)*.5d0 ; in [0,2Pi[
            dphi = (phie-phiw)*.5d0
            ;print,nside,izlist[0],max(izlist),k,n_elements(dphilist),i,minmax(subrings)
            dd = abs(phi0list[subrings]-phic) ; distance from disc center to pixel border sample
            dd <= (2.d0 * !dpi - dd) ; in [0,Pi]
            touching = max(dd le (dphilist[subrings]+dphi)) ; 0:out or 1:in
            if (~touching) then begin;  debugged 2020-08-28
                ;ringphi[i,kk] -= k ; move edge pixel inwards
                ringphi[i,kk] = (ringphi[i,kk] + np - k) MOD np ; move edge pixel inwards in [0,np-1]
                nshrink += 1
                if (nshrink le 2) then begin
                    goto, rep   ; repeat with next big pixel inward
                endif else begin
                    ringphi[i,kk] = -999
                endelse
            endif
        endif
    endfor
endfor
; print,' ***'
; print, ringphi

; remove empty rings
; diff = ringphi[*,2] - ringphi[*,1]
; bug correction 2020-08-28
;valid = where(diff ne -2 and diff ne -1 and ringphi[*,1] ge 0 and ringphi[*,2] ge 0, ngr)
valid = where(ringphi[*,1] ge 0 and ringphi[*,2] ge 0, ngr)
if (ngr gt 0) then begin
    ringphi = ringphi[valid,*]
endif else begin
    ringphi = -1L
endelse

return
end


;===============================================================


pro query_disc, nside, vector0, radius_in, listpix, nlist, deg = deg, inclusive=inclusive, help=help, nested=nested, walltime=walltime, boost=boost
; -----------------------------------------------------------------------------
;+
; NAME:
;  query_disc
;
; PURPOSE:
;  finds the Healpix pixels that lie within a radius (radius_in) from the
;  vector (vector_0) in the ring scheme for the resolution Nside
;
; CATEGORY:
;  Healpix pixel toolkit
;
; CALLING SEQUENCE:
;     QUERY_DISC, Nside, Vector0, Radius_In, Listpix, Nlist, [Deg=, Nest=, Help=,
;     Inclusive=, Walltime=]
;
; INPUTS:
;     Nside : scalar integer : Healpix resolution (power of 2)
;     Vector0 : 3-element vector (float or double) : position of the
;          disc center on the sphere (north pole in [0.,0.,1.])
;          the norm of Vector0 does not have to be one, what is
;          consider is the intersection of the sphere with the line of
;          direction Vector0
;     Radius_in : radius of the circle (in radians, unless Deg is set)
;
; KEYWORD PARAMETERS:
;     Deg : if set, the disc radius is in degrees instead of radians.
;
;     /Help: if set, prints this documentation header and exits
;
;     /INCLUSIVE :0 by default, only the pixels whose center 
;                       lie in the triangle are listed on output
;                  if set to 1, all pixels overlapping the triangle are output
;
;     /NESTED  :0 by default, the output list is in RING scheme
;                  if set to 1, the output list is in NESTED scheme
;
; OUTPUTS:
;     Listpix : list of pixels found within a radius Radius_in from
;          Vector0,
;     Nlist = number of elements of Listpix
;       if no pixel is found (too small a circle) Nlist is 0 and
;       listpix is -1 .
;
; OPTIONAL OUTPUTS:
;     WALLTIME=  contains on output the wall clock time used by the routine [s]   
;
; COMMON BLOCKS:
;     none.
;
; RESTRICTIONS:
;     On output Listpix is (supposed to be) ranked in pixel index
;     rather than ranked in distance to the disc center.
;
; SIDE EFFECTS:
;     calls : ring_num, discphirange_at_z, ring2nest, modulo
;          and   pixels_in_phi_range2 (included in this file)
;
; PROCEDURE:
;     uses the particular layout of the pixels in parallel rings,
;     to find rapidly the selected pixels.
;     Is much faster than a simple search among all the pixels.
;
; EXAMPLE:
;     query_disc, 256L, [.5,.5,0.], 10., listpix, nlist, /Deg
;       outputs in listpix the Healpix pixel numbers of the pixels
;       with 10 deg of the point on the sphere 
;       having the direction [.5,.5,0.]
;
; MODIFICATION HISTORY:
;     1998??     EH, TAC, 1st f90 version
;     1999-??    EH, Caltech, traduction in IDL
;     1999-12-07 : correction of a bug in the build in subroutine in_ring
;     2002-09    : added inclusive keyword, added nest keyword, renamed
;     query_disc
;     2008-03-30: fixed bug appearing when disc centered on either pole
;                 use 'work' array instead of expanding output array piece by piece
;     2009-04-08: actually returns -1 if nlist = 0
;     2011-09-22: improved performance (fewer false positive) in /INCLUSIVE mode
;          by testing border of edge pixels (sampled at larger Nside)
;     2013-01-14: avoid crashes is some configurations with empty Listpix list,
;            systematically returns Listpix=[-1], Nlist=0 in case of problem
;     2014-07-07: bug correction in discedge2fulldisc (reported by F. Atrio-Barandella)
;     2020-08-28: bug correction in check_edge_pixels (reported by R. de Belsunce)
;-

tstart = systime(1)
code = 'query_disc'
syntax= 'SYNTAX = '+code+', Nside, Vector0, Radius_In, Listpix, [Nlist, DEG=, HELP=, INCLUSIVE=, NESTED=, WALLTIME=]'
nlist = 0 & listpix = [-1]

if keyword_set(help) then begin
    doc_library,code
    return
endif

if (n_params() lt 4 || n_params() gt 5) then begin
    print,syntax
    return
endif

prompt = strupcase(code)+'> '
if (n_elements(vector0) ne 3) then begin
    print,prompt+'vector0 should be a 3 element vector'
    return
endif

if (radius_in lt 0.) then begin
    print,prompt+'radius should be > 0 : ',radius_in
    return
endif

if (keyword_set(deg) && radius_in gt 180.) then begin
    print,prompt+'radius (deg) too big :',radius_in
    return
endif

if (~keyword_set(deg) && radius_in gt !pi) then begin
    print,prompt+'radius (radian) too big :',radius_in
    return
endif

do_inclusive = keyword_set(inclusive)

npix = nside2npix(nside,err=errpix)
if (errpix) then begin
    print,prompt+'invalid Nside:',nside
    return
endif

defsysv,'!healpix',exists=yes
if (~yes) then init_healpix

;radius is in radian
radius = keyword_set(deg) ? radius_in*!DPI/180.d0 : radius_in

halfpi = !DPI*.5d0
twopi  = !DPI*2.d0
lnside = long(nside)

radius_eff = do_inclusive ? fudge_query_radius(lnside, radius) : radius
cosang     = cos(radius_eff)

; create work array used to store valid pixels
; (faster than expanding final array)
worksize = long64( (1.d0 - cosang)/2.d0 * 1.2 * npix  + 4*lnside + 50)
worksize = worksize < npix

if (lnside gt 8192) then begin
    work = lon64arr(worksize)
endif else begin
    work = lonarr(worksize)
endelse
;time1 = systime(1)

; circle center
norm_vect0 = sqrt(total(vector0^2))
z0 = vector0[2]/norm_vect0

; find upper and lower z
rlat0 = asin(z0) ; lat in RAD
rlat1 = rlat0 + radius_eff
rlat2 = rlat0 - radius_eff
zmax = (rlat1 ge halfpi)  ? 1.d0  : sin(rlat1)
zmin = (rlat2 le -halfpi) ? -1.d0 : sin(rlat2)

; fills list of rings and dphi
irmin = ring_num(lnside, zmax, shift=+1)
irmax = ring_num(lnside, zmin, shift=-1)
nz = (irmax-irmin+1) > 1  ; bug correction 2013-01
izlist = irmin + lindgen(nz)    ; list of rings
zlist  = ring2z(lnside, izlist) ; list of z
dphilist = discphirange_at_z (vector0, radius_eff, zlist, phi0=phi0) ; phi range in each ring
phi0list = replicate(phi0, nz)
; identify edge pixel at nominal resolution
pixels_on_edge, lnside, izlist, phi0list, dphilist, ringphi, ngr
if do_inclusive then begin
    boost_def  = 2L^ ceil( alog(1.1d0/(radius*lnside))/alog(2.d0) ) ; make pixel smaller than disc (2020-08-31)
    boost_def >= 16 ; minimal boost
    nsboost    = defined(boost) ? boost : boost_def ; user defined value of boost
    wnside = (lnside * nsboost) <  max(!healpix.nside)
    radius2 = fudge_query_radius(wnside, radius, /quadratic)

    irmin = ring_num(wnside, zmax, shift=+1)
    irmax = ring_num(wnside, zmin, shift=-1)
    nz = (irmax-irmin+1) > 1  ; bug correction 2013-01
    izlist = irmin + lindgen(nz) ; list of active rings
    zlist  = ring2z(wnside, izlist) ; list of z
    dphilist = discphirange_at_z (vector0, radius2, zlist, phi0=phi0) ; phi range in each ring
    phi0list = replicate((phi0 + twopi) mod twopi, nz) ; in [0,2Pi] (2020-08-28)

;   check boundary of edge pixels
;   ringphi and ngr computed in non-inclusive configuration (with lnside) will be updated
    check_edge_pixels, lnside, wnside, izlist, phi0list, dphilist, ringphi, ngr
endif

;time2 = systime(1)
nlist = 0LL
discedge2fulldisc, lnside, ringphi, ngr, work, nlist

if (nlist gt 0) then begin
    if keyword_set(nested) then ring2nest, lnside, work[0:nlist-1], listpix else listpix=work[0:nlist-1]
endif else begin
    listpix = -1L
endelse

tend = systime(1)
walltime = tend - tstart
;;print,'total, init, find, fill:',walltime*1.e6,(time1-tstart)*1.e6,(time2-time1)*1.e6,(tend-time2)*1.e6
return
end

