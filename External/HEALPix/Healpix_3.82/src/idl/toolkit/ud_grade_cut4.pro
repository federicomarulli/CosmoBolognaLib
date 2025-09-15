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
pro ud_grade_cut4, pixel, signal, n_obs, serror, nside_in = nside_in, $
          nside_out = nside_out, order_in = order_in, order_out = order_out, pessimistic = pessimistic
;+
;
; UD_GRADE_CUT4,    $
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
;     e_i  -> ef_j = sqrt[ sum_i(n_i(t_i-tf_j)^2) + sum_i (n_i (e_i)^2) / sum_i (1)] / sum_i (n_i)
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
;    2009-03-27: fixed bug in prograding, allows large Nside
;                fixed bug in Serror degradation
;                degrading now works for large Nside
;-
tstart = systime(1)
routine = 'UD_GRADE_CUT4'
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

ord_in  = strupcase(strmid(order_in, 0,4))
ord_out = strupcase(strmid(order_out,0,4))
check = ['RING','NEST']
junk = where(check eq ord_in,  n_in)
junk = where(check eq ord_out, n_out)
if (n_in  eq 0) then message,'Invalid Input  ordering '+order_in
if (n_out eq 0) then message,'Invalid Output ordering '+order_out


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
do_obs = defined(n_obs)
do_err = defined(serror)

if (nside_out lt nside_in) then begin 
; ------------
; degradation
; ------------
    ratio = long(nside_in/nside_out)
    ratio2 = ratio * ratio

    if (ord_in eq 'RING') then begin
        ring2nest, nside_in, pixel, pixel
    endif

    minhit = 1                                         ; at least one good sub-pixel
    if (keyword_set(pessimistic)) then minhit = ratio2 ; all sub pixel should be good

    ; subdivide maps in smaller areas (containing an integer number of big output pixels)
    ;nside_wrk = (nside_in / 64) > 1 < nside_out
    nside_wrk = 16
    nsub = nside2npix(nside_wrk)
    subsize = npix_out / nsub

    ; get an upper limit on the size of ouput map
    ; and create arrays to collect results
    h1 = histogram(pixel/ratio2, min=0, max=npix_out-1, binsize=1)
    npo = long64(total(h1 ge minhit))
    l64 = (nside_out gt 8192)
    pixel_cumul  = make_array(npo, l64=l64, long=1-l64)
    signal_cumul = make_array(npo, /float)
    if (do_obs) then begin
        l64_obs = ((max(n_obs) * double(ratio2)) gt 2.e9)
        obs_cumul = make_array(npo, l64=l64_obs, long=1-l64_obs)
    endif
    if (do_err) then err_cumul = make_array(npo, /float)

    ; identify small input pixels belonging to each area
    hist = histogram(pixel/ratio2, min=0, max=npix_out - 1, binsize= subsize, locations=locations, rev=rev,/l64)

    iter = 0
    obs_npix = 0LL
    for isub=0L, nsub-1 do begin
        ; find out if output submap contains some input pixels
        ngood = hist[isub]
        if (ngood gt 0) then begin
            good_in = rev[rev[isub]:rev[isub+1]-1]
            lowsub = locations[isub]

            pix0 = lowsub * ratio2
            local = pixel[good_in] - pix0
            hitin = bytarr(ratio2, subsize)
            hitin[local] = 1B
            hitout = total(hitin,1)
            if (ok_obs) then begin
                obsin = lonarr(ratio2, subsize)
                obsin[local] = n_obs[good_in]
                obsout = round(total(obsin,1,/double),l64=l64_obs) ; total number of hits per large pixel
            endif else begin
                obsin = hitin
                obsout = hitout ; number of exposed subpixels per large pixel
            endelse
            
            mapin = fltarr(ratio2, subsize)
            mapin[local] = signal[good_in]
            mapout = total(obsin * mapin,1,/double)
            
            pixelout  = where(obsout gt 0 and hitout ge minhit, obs_npix1)
            if (obs_npix1 gt 0) then begin
                signalout = mapout[pixelout]/obsout[pixelout] ; average temperature per large pixel
            
                if do_err  then begin
                    mapav = fltarr(ratio2, subsize)
                    errin = fltarr(ratio2, subsize)
                    mapav[0:ratio2-1,pixelout] = replicate(1.0,ratio2)#signalout
                    errin[local] = serror[good_in]
                    count = total(obsin gt 0, 1)>1.e-6 ; number of non-empty small pixels in big pixel
                    errout = sqrt(   total( obsin*errin*errin,1,/double)/count  + total(obsin*(mapin-mapav)^2,1,/double)     )
                endif

                pixel_cumul [obs_npix:obs_npix+obs_npix1-1] = pixelout +lowsub
                signal_cumul[obs_npix:obs_npix+obs_npix1-1] = float(signalout)
                if (do_obs) then obs_cumul[obs_npix:obs_npix+obs_npix1-1] = obsout[pixelout]
                if (do_err) then err_cumul[obs_npix:obs_npix+obs_npix1-1] = float(errout[pixelout])
            endif
            obs_npix = obs_npix + obs_npix1
        endif
    endfor
    mapin = 0
    obsin = 0
    pixelout = 0 & signalout = 0 & obsout = 0 & errout = 0
    if (obs_npix eq 0) then begin
        message,/info,'Output map is empty'
        return
    endif
    pixel  = pixel_cumul[0:obs_npix-1]   & pixel_cumul  = 0
    signal = signal_cumul[0:obs_npix-1]  & signal_cumul = 0
    if do_obs then n_obs  = obs_cumul[0:obs_npix-1]
    if do_err then serror = err_cumul[0:obs_npix-1]/obs_cumul[0:obs_npix-1]
    obs_cumul = 0
    err_cumul = 0

    nside = nside_out
    npix = npix_out
    
    if (ord_out eq 'RING') then begin
        nest2ring, nside, pixel, pixr
        
        ind = sort(pixr)
        pixel = pixr[ind]
        signal = signal[ind]
        if do_obs then n_obs = n_obs[ind]
        if do_err then serror = serror[ind]
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

    if do_obs then begin
        n_obs  = n_obs  ## (subpix/ratio2)
        n_obs  = reform(round(n_obs) , npix_out, /overwrite)
    endif
    if do_err then begin
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
        if do_obs then n_obs = n_obs[ind]
        if do_err then serror = serror[ind]
    endif

endelse
;print,'Time [s]:',systime(1)-tstart
    
return
end
