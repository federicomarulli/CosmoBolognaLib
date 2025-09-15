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
; -----------------------------------------------------------------------------
pro ud_grade_cut4_old, pixel, signal, n_obs, serror, nside_in = nside_in, $
          nside_out = nside_out, order_in = order_in, order_out = order_out, pessimistic = pessimistic
;+
;
; UD_GRADE_CUT4_OLD,    $
;     Pixel, Signal [, N_Obs, Serror],   $
;     Nside_in=, [Nside_out=,] Order_in=, [Order_out=], [Pessimistic=]
;
;
; up_gradation, de_gradation for the fields
;  Pixel (pixel index), 
;  Signal (pixel temperature = average of several observations), 
;  N_obs (number of observations in pixel), 
;  Serror (error on temperature = rms of average)
;
;
;  degradation :  K small pixels (i) -> 1 big pixel (j)
;    N_Obs :    n_i  -> nf_j = sum_i (n_i)
;    Signal:    t_i  -> tf_j = sum_i (n_i t_i)  / sum_i (n_i)
;    Serror:    
;     e_i  -> ef_j = sqrt[ sum_i(n_i(t_i-tf_j)^2) + sum_i (n_i s_i)^2] / sum_i n_i
;    Pixel : at least one exposed small pixel -> 1 exposed big pixel    
;
;  upgradation : 1 big pixel(j) -> K small *identical* pixels (i)
;    N_obs :    nf_j   -> n_i = nf_j / K (approximated to closest integer)
;    Signal:    tf_j   -> t_i = tf_j
;    Serror:    ef_j   -> e_i = ef_j * sqrt(K)
;
;
;   EH, 2000-11-23
;    June 2003,  EH, replaced STOPs by MESSAGEs
;
;    2009-03-25: fixed bug in prograding, allow large Nside
;
;-
tstart = systime(1)
routine = 'UD_GRADE_CUT4_OLD'
if (undefined(nside_in) or undefined(order_in)) then begin
    print,routine+': Pixel, Signal [, N_Obs, Serror], Nside_in=, [Nside_out=,] Order_in=, [Order_out=]'
    print,'******** No processing ********'
    return
endif
if undefined(nside_out) then nside_out = nside_in
if undefined(order_out) then order_out = order_in
npix_in = nside2npix(nside_in)
npix_out = nside2npix(nside_out)
if (npix_in le 0 or npix_out le 0) then begin
    print,routine+': Invalid Nside_in or Nside_out :',nside_in,nside_out
    message,'Abort'
endif

ord_in  = strmid(order_in,0,4)
ord_out = strmid(order_out,0,4)

if ((ord_in eq ord_out) and (nside_in eq nside_out)) then begin
; ------------
; nothing
; ------------
    return
endif


ok_obs = 0
if defined(n_obs) then begin
    if (max(n_obs) gt 0) then ok_obs = 1
endif
if (nside_out lt nside_in) then begin 
; ------------
; degradation
; ------------
    ratio = (nside_in/nside_out)
    ratio2 = ratio * ratio

    if (ord_in eq 'RING') then begin
        ring2nest, nside_in, pixel, pixel
    endif

    minhit = 1                                         ; at least one good sub-pixel
    if (keyword_set(pessimistic)) then minhit = ratio2 ; all sub pixel should be good

    hitin = bytarr(ratio2,npix_out)
    hitin[pixel] = 1B
    hitout = total(hitin,1)
    if (ok_obs) then begin
        obsin = lonarr(ratio2,npix_out)
        obsin[pixel] = n_obs
        obsout = total(obsin,1) ; total number of hits per large pixel
    endif else begin
        obsin = hitin
        obsout = hitout         ; number of exposed subpixels per large pixel
    endelse
    
    mapin = fltarr(ratio2,npix_out)
    mapin[pixel] = signal
    mapout = total(obsin * mapin,1)
    
    pixelout  = where(obsout gt 0 and hitout ge minhit, obs_npix)
    signalout = mapout[pixelout]/obsout[pixelout] ; average temperature per large pixel
    
    if defined(serror)  then begin
        mapav = fltarr(ratio2,npix_out)
        errin = fltarr(ratio2,npix_out)
        ;;;mapav[pixel] = signalout[pixel/ratio2]
        mapav[0:ratio2-1,pixelout] = replicate(1.0,ratio2)#signalout
        errin[pixel] = serror
        count = total(obsin gt 0, 1)>1.e-6 ; number of non-empty small pixels in big pixel
        errout = sqrt(   total( obsin*errin*errin,1,/double)/count  + total(obsin*(mapin-mapav)^2,1,/double)     )
        errin = 0
        mapav = 0
    endif
    
    mapin = 0
    obsin = 0
    
    pixel = temporary(pixelout)
    signal = temporary(signalout)
    if defined(n_obs)  then n_obs = obsout[pixel]
    if defined(serror) then serror = float(errout[pixel])/obsout[pixel]
    nside = nside_out
    npix = npix_out
    
    if (ord_out eq 'RING') then begin
        nest2ring, nside, pixel, pixr
        
        ind = sort(pixr)
        pixel = pixr[ind]
        signal = signal[ind]
        if defined(n_obs) then n_obs = n_obs[ind]
        if defined(serror) then serror = serror[ind]
        ind = 0
    endif
endif else begin
; ------------
; upgradation
; ------------
    ratio = (nside_out/nside_in)
    ratio2 = ratio * ratio

    ; go to NEST
    if (ord_in eq 'RING') then begin
        ring2nest, nside_in, pixel, pixel
    endif

    npix = n_elements(pixel)
    npix_out = long64(npix) * ratio2

    increm = lindgen(ratio2)
    if ((nside_out le 8192) && (size(/type,pixel) ne 14)) then begin
        pixel_out = lonarr(ratio2,npix) ; corrected 2009-03-25
        lratio2   = ratio2
    endif else begin
        pixel_out = lon64arr(ratio2,npix)
        lratio2   = long64(ratio2)
    endelse
    for i=0LL, npix-1LL do begin
        pixel_out[*,i] = pixel[i] * lratio2 + increm
    endfor
    pixel  = reform(pixel_out, npix_out)
    pixel_out = 0
    increm = 0

    subpix = replicate(1.,ratio2)
    signal = signal ## (subpix)       
    signal = reform(signal, npix_out, /overwrite)

    if defined(n_obs) then begin
        n_obs  = n_obs  ## (subpix/ratio2)
        n_obs  = reform(round(n_obs) , npix_out, /overwrite)
    endif
    if defined(serror) then begin
        serror = serror ## (subpix*ratio) 
        serror = reform(serror, npix_out, /overwrite)
    endif
    npix = npix_out
    nside = nside_out

    ; go to output order
    if (ord_out eq 'RING') then begin
        nest2ring, nside, pixel, pixel
    endif

    ; sort 
    if (ord_in eq 'RING' or ord_out eq 'RING') then begin
        ind = sort(pixel)
        pixel = pixel[ind]
        signal = signal[ind]
        if defined(n_obs) then n_obs = n_obs[ind]
        if defined(serror) then serror = serror[ind]
    endif

endelse
print,'Time [s]:',systime(1)-tstart
    
return
end
