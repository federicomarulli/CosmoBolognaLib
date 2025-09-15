;===============================================
function unfold_weights, arg1, dim,  $
  directory = dir_usr, $
  help      = help,    $
  pixel     = pixel,   $
  ring      = ring,    $
  scheme    = scheme,  $
  silent    = silent
;+
; NAME:
;       UNFOLD_WEIGHTS
;
; PURPOSE:
;      Unfold the weights to be applied to a HEALPix map in order to improve the quadrature.
;      The input weights are either ring based (the same for all pixels on the same iso-latitude ring),
;      and packed in a FITS file with 2*Nside values centered on 0,
;      OR pixel based (different for each pixel, up to some symmetries),
;      and packed in a FITS file with ~ 0.75 * Nside^2 values centered on 0.
;      After unpacking, the ouput is a full sky map with Npix = 12*Nside^2 values centered on 1.
;
; CATEGORY:
;
; CALLING SEQUENCE:
;   2 alternative sequences are available
;    Weight_map = unfold_weights (Nside, [PIXEL=|RING=|SCHEME=] [, Dim, DIRECTORY=, HELP=, SILENT=])'
;    Weight_map = unfold_weights (File                          [, Dim,             HELP=, SILENT=])'
;
; INPUTS:
;    File:  input weight file to be read (character string)
;    Nside: HEALPix resolution parameter (integer number, power of 2),
;       the routine determines the file to be read using Nside, DIRECTORY, 
;       and RING, PIXEL or SCHEME
;
; OPTIONAL INPUTS:
;    Dim:   1 or 2, dimension of output, default=1
;        if Dim=1, the output is a  vector of size   Npix
;        if Dim=2, the output is an arrray of shape (Npix,3) with 3 identical columns
;
; KEYWORD PARAMETERS:
;    DIRECTORY: directory in which to look for the weight file (default: !healpix.path.data)
;
;    HELP:  if set, this documentation header is printed out
;
;    PIXEL: if set, the code will look for the pixel-based weight file corresponding the the Nside provided,
;         in the default or provided directory 
;    RING:  if set, the code will look for the ring-based  weight file corresponding the the Nside provided,
;         in the default or provided directory 
;    SCHEME: can be either 'PIXEL' or 'RING', setting the type of weight file the code will look for.
;
;    SILENT: if set, the routine works silently
;
;
; OUTPUTS:
;    Weight_map: full sky weight map
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;     reads a weight FITS file, and expand its content into a full sky map
;
; EXAMPLE:
;     mollview,unfold_weights(16,/ring), title='Ring-based weights', /silent
;     mollview,unfold_weights(16,/pixel),title='Pixel-based weights',/silent
;   will show respectively the ring-based and pixel-based weights 
;   located in !healpix.path.data and applicable to HEALPix maps with Nside=16.
;
; MODIFICATION HISTORY:
;    2017-12-22: created
;    2018-01-08: documented
;
;-

routine = 'unfold_weights'
syntax  = ['Weight_map = '+routine+' (nside, [PIXEL=|RING=|SCHEME=] [, dim, DIRECTORY=, HELP=])',$
           'Weight_map = '+routine+' (file                          [, dim,             HELP=])']

if keyword_set(help) then begin
    doc_library,routine
    return,-1
endif

if n_params() eq 0 then begin
    print,syntax
    print,'with default directory= '+!healpix.path.data
    return,-1
endif

be_silent = keyword_set(silent)

; ------------------------------------
; look for relevant file, and read it
; ------------------------------------

if (size(arg1,/tname) eq 'STRING') then begin    ; file given by user

    file = arg1
    ; read file
    if (~file_test(file)) then begin
        print,'Could not find the file '+file
        message,'Abort'
    endif
    read_fits_map, file, w8, h, xh, nside=nside, /silent
    nw   = n_elements(w8[*,0])
    do_ring = (nw eq 2*nside)
    do_pix  = (nw eq nside2npweights(nside))

    if ~(do_ring || do_pix) then begin
        print,'Input file '+file+' is neither a ring-based nor pixel-based weight file'
        message,'Abort'
    endif

endif else begin    ; guess file name from Nside

    nside = arg1
    uscheme = defined(scheme) ? strupcase(strmid(scheme,0,1)) : ''
    do_ring = keyword_set(ring)  || uscheme eq 'R'
    do_pix  = keyword_set(pixel) || uscheme eq 'P'
    if (do_ring eq do_pix) then begin
        print,'Choose one and only one among '
        print,"- Ring  (with ring=1  or scheme='RING')  "
        print,"- Pixel (with pixel=1 or scheme='PIXEL') "
        message,'Abort'
    endif
    
    
    defsysv, '!healpix', exists = exists
    if (exists ne 1) then init_healpix
    
    ns = long(nside)
    
    ; check nside
    err = 0
    npix = nside2npix(ns,err=err)
    
    ; test, complain and exit
    if (err ne 0) then begin
        print,routine +': invalid Nside ',nside
        return,-1
    endif
    
    if (ns le 1) then begin
        print,routine +print,'Nside out of range ',nside
        return,-1
    endif
    
    
    ; form file name
    if (do_ring) then begin
        snside = string(ns,form='(i5.5)')
        file1 = 'weight_ring_n'+snside+'.fits'
    endif else begin
        snside = string(ns,form='(i5.5)')
        file1 = 'weight_pixel_n'+snside+'.fits'
    endelse
    
    ; look for file
    defdir = defined(dir_usr) ? dir_usr : !healpix.path.data
    files  = defdir+['','/']+file1
    found  = 0
    for i=0,n_elements(files)-1 do begin
        file = files[i]
        found = found || file_test(file)
    endfor
    if (~found) then begin
        print,'Could not find the file '+file1
        print,' in the directory '+defdir
        if ~defined(dir_usr) then begin
            print,'Are you sure the Healpix system variable is correcly defined'
            print,'current value HEALPIX = '+getenv('HEALPIX')
        endif
        message,'Abort'
    endif
    
    ; read file
    if ~be_silent then begin
        print,'Reading weights from '+file
    endif
    read_fits_map, file, w8, h, xh, nside=nside1, /silent
    nw   = n_elements(w8[*,0])

    if (nside1 ne nside) then begin
        print,'Expected Nside = ',nside
        print,'Found ',nside1,' in '+file
        message,'Abort'
    endif

endelse

; -------------------------
; do the actual unfolding
; -------------------------

if (do_ring) then begin

    if (nw ne 2*nside) then begin
        print,'Expected ',strtrim(2*nside,2),' ring weights.'
        print,'Read ',strtrim(nw,2),' instead.'
        message,'Aborting.'
    endif

    twod = (defined(dim) && dim ge 2)
; fill ring map
    npix = nside2npix(nside)
    map  = (twod) ? dblarr(npix,3) : dblarr(npix)
    start = 0LL
    n4 = 4L * nside
    for ir=1L,n4-1 do begin     ; loop on all rings
        it = ir < (n4 - ir)     ; ring index, counting from closest pole
        np = 4L* (it < nside)   ; pixels/ring
        if (twod) then begin
            for k=0,2 do map[start:start+np-1,k] = w8[it-1,k]
        endif else begin
            map[start:start+np-1] = w8[it-1,0]
        endelse
        start += np
    endfor

endif else begin

    nf = nside2npweights(nside)
    if (nw ne nf) then begin
        message,/info,' expected '+strtrim(nf,2)+' weights for Nside='+strtrim(nside,2)+$
                ' and got '+strtrim(nw,2)
        message,'Aborting.'
        return,-1
    endif

    twod = (defined(dim) && dim ge 2)
    npix = nside2npix(nside)
    pnorth  = npix*0            ; position in expanded weights
    vpix    = npix*0            ; position in compress list

    map  = (twod) ? dblarr(npix,3) : dblarr(npix)
    for ring=0L, 2*nside-1 do begin ; loop on Northern rings
        qpix    = (ring+1) < nside ; 1/4 of number of pixels per ring
        shifted = ring lt (nside-1L) || (ring+nside) mod 2
        odd     = qpix mod 2
        qp4     = 4L*qpix       ; number of pixels per ring
                                ; fill the weight map
        if (qpix le nside) then begin ; reuse the same ones in the tropical region
            ramp = lindgen(qp4)
            j4   = ramp mod qpix
        endif
        rpix = j4 < (qpix - shifted - j4) ; seesaw
        if (twod) then begin
            for k=0,2 do map[pnorth,k] = w8[vpix + rpix]
        endif else begin
            map[pnorth] = w8[vpix + rpix]
        endelse
        if (ring lt 2*nside-1) then begin ; south part (except equator)
            psouth     = npix - pnorth - qp4
            if (twod) then begin
                for k=0,2 do map[psouth,k] = map[pnorth:pnorth+qp4-1,k]
            endif else begin
                map[psouth] = map[pnorth:pnorth+qp4-1]
            endelse
        endif
                                ; locations on next ring
        pnorth  += qp4
        wpix     = (qpix+1L)/2  + 1 - (odd || shifted)
        vpix    += wpix
    endfor
    
endelse

map += 1.d0

return, map
end



