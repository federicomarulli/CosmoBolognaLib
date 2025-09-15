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
;+
; NAME:
;       Median_Filter
;
; PURPOSE:
;       Produces the median filtered map of an Healpix map
;
; CATEGORY:
;
; CALLING SEQUENCE:
;       median_filter, InputMap, Radius, MedianMap [, ORDERING=, /RING,
;       /NESTED, /FILL_HOLES, /DEGREES, /ARCMIN]
;
; INPUTS:
;       InputMap:   either an IDL array containing a Healpix map to filter, 
;       or a FITS file name containing the map to filter
;
;       Radius: radius of the disk on which the median is computed.
;            It is in Radians, unless /DEGREES or /ARCMIN are set
;      
; KEYWORD PARAMETERS:
;
;  ORDERING=  Healpix map ordering, should be either 'RING' or 'NESTED'
;  /RING      same as ordering='RING'
;  /NESTED    same as ordering='NESTED'
;
;  BAD_DATA=:  flag value of missing pixels (default=!healpix.bad_value).
;              if flagged pixels are enclosed in the disc used to compute the
;              median they will be ignored (and the disc radius remain unchanged).
;              if the central pixel is flagged itself, it will remain unchanged,
;              unless /FILL_HOLES is set, in which case it will be replaced by the
;              median of the disc
;
;  /FILL_HOLES       if set, missing pixels will take the value of the local median
;              (if it is defined),
;              otherwise leave bad pixel unchanged
;
;  /DEGREES    if set, the Radius is understood to be provided in Degrees
;
;  /ARCMIN     if set, the Radius is understood to be provided in Arcmin
;
;
; OUTPUTS:
;       MedianMap: either an IDL variable containing on output the filtered map,
;       or a FITS file name into which the filtered map will be written
;
; SIDE EFFECTS:
;      will write a FITS file if InputMap and MedianMap are FITS file names
;
; RESTRICTIONS:
;      InputMap and OutputMap should be of the same type (both arrays or both filenames)
;
; PROCEDURE:
;      loop over pixels and replace by median value of disc centered on that pixel
;
; EXAMPLE 1:
;      median_filter, 'map1.fits', 10, 'map1_mf.fits'
;
; EXAMPLE 2:
;      map = randomn(seed, nside2npix(128))
;      median_filter, map, 10, map_mf
;
; MODIFICATION HISTORY:
;      EH, IPAC, 2005-April
;      March 2006: M.A. Miville-Deschenes (IAS) and EH (IAP): bugs correction
;      July 2006: EH corrected definition of do_nest (pointed out by Kwangil Seon)
;      May 2017: replaced obsolete DATATYPE() with IDL's SIZE(/TNAME)
;      Frb 2020: use FINITE() to deal with NaN; deal with partial sky maps
;-

pro medfilt_sub, in_map, radius_radian, med_map, ring=ring, nested=nested, ordering=ordering, $
                 bad_data = bad_data, fill_holes=fill_holes


NaN    = !values.f_nan
dofill = keyword_set(fill_holes)
nmaps  = n_elements(in_map[0,*])
npix   = n_elements(in_map[*,0])
nside  = npix2nside(npix, error=error)
if (nside lt 0) then begin
    message,'Unvalid number of pixels in Healpix median filter: '+strtrim(npix,2)
endif
opix  = 4.d0*!dpi / npix
odisc = 4.d0*!dpi * sin(radius_radian/2.d0)^2
print,odisc / opix,' pixels per disc in average'
; count = lonarr(npix)

; replace bad pixels by NaN
if undefined(bad_data) then bad_data = !healpix.bad_value
bad = where(in_map eq bad_data, nbad)
if (nbad gt 0) then in_map[bad] = NaN

ord = keyword_set(ordering) ? strupcase(strmid(ordering,0,4)) : ''
do_nest = (keyword_set(nested) or ord eq 'NEST')
do_ring = (keyword_set(ring) or ord eq 'RING')

if (do_ring eq do_nest) then begin
    message,'choose map ordering: either RING or NESTED'
endif

; denoised map
med_map = in_map

; loop on pixel centers
for p=0L, npix-1 do begin

    ; find pixel location
    if (do_ring) then pix2vec_ring, nside, p, vector
    if (do_nest) then pix2vec_nest, nside, p, vector
    
    ; find disc centered on pixel
    query_disc, nside, vector, radius_radian, list, nlist, nested=do_nest
;     if ((p+1) mod (npix/12) eq 0) then print,p,npix, nlist
;     count[p] = nlist

    ; replace center pixel value by median value within disc
    ; Nan pixels are ignored in median calculation
    ; if even number of valid pixels, take average of 2 pixels closest to median
    if (nlist gt 0) then begin
        for j=0, nmaps-1 do begin
            if (finite(in_map[p, j]) or dofill) then begin
                med_map[p, j] = median(in_map[list,j], /even)
            endif
        endfor
    endif

endfor

; replace NaN pixels by original bad pixel flag
if (nbad gt 0 and finite(bad_data)) then begin
    in_map[bad] = bad_data
    bad2 = where(med_map eq bad_data or finite(med_map,/nan), nbad2)
    if (nbad2 gt 0) then med_map[bad2] = bad_data
endif

; print,moment(count)

return
end

;--------------------------------------------------------------------------------------------------

pro median_filter, datain, radius, medmap, ring=ring, nested=nested, ordering=ordering, fill_holes=fill_holes, degrees=degrees, arcmin=arcmin, bad_data=bad_data
code = 'MEDIAN_FILTER'

syntax = code+', InputMap, Radius, MedianMap [, ORDERING=, /RING, /NESTED, /FILL_HOLES, /DEGREES, /ARCMIN]'
if (n_params() ne 3) then begin
    print,syntax
    return
endif

if (size(/tname,datain) eq 'STRING' xor size(/tname,medmap) eq 'STRING') then begin
    message,'Inconsistent input and output data: should both be either FITS files or arrays'
endif

if (keyword_set(arcmin) && keyword_set(degrees)) then begin
    message,'Choose either DEGREES or ARCMIN'
endif

radius_radian = radius
if (keyword_set(degrees)) then radius_radian = radius * !DtoR
if (keyword_set(arcmin )) then radius_radian = radius * !DtoR / 60.
if (radius_radian lt 0 or radius_radian gt !pi) then begin
    message,'Invalid Radius provided: '+strtrim(radius)
endif
radius_deg = radius_radian * !Radeg

if size(/tname,datain) ne 'STRING' then begin
    ; data are memory array
    medfilt_sub, datain, radius_radian, medmap, ring=ring, nested=nested, ordering=ordering, $
      fill_holes=fill_holes, bad_data=bad_data
endif else begin
    ; data is read from FITS file and output is written on another FITS file
    file_in = datain
    npix = getsize_fits(file_in, type=type, nmaps=nmaps)
    case type of
        2: begin ; full sky
            read_fits_map, file_in, data_in, h, xh, /silent, nside=nside, ordering=ordering
        end
        3: begin ; cut sky
;             read_fits_cut4, file_in, pixel, signal, n_obs, serror, hdr=h, xhdr=xh, nside=nside, ordering=ordering
;             data_in = replicate(!healpix.bad_value, nside2npix(nside))
;             data_in[pixel] = signal
            names  = sxpar(headfits(file_in, ext=1),'TTYPE*')
            cut4 = array_equal(strupcase(strmid(names[0:3],0,1)),['P','S','N','S'])
            if (cut4) then begin
                read_fits_cut4, file_in, pixel, signal, n_obs, serror, hdr=h, xhdr=xh, nside=nside, ordering=ordering
                data_in = replicate(!healpix.bad_value, nside2npix(nside))
                data_in[pixel] = signal
            endif else begin
                read_fits_partial, file_in, pixel, tmp, hdr=h, xhdr=xh, nside=nside, ordering=ordering
                ncol = (size(tmp,/dim))[1]
                data_in = replicate(!healpix.bad_value, nside2npix(nside), ncol)
                help,pixel,tmp,data_in
                for i=0,ncol-1 do data_in[pixel, i] = tmp[*,i]
            endelse
        end
        else: begin
            message,'wrong file type for '+file_in
        end
    endcase
    ; do actual filtering
    medfilt_sub, data_in, radius_radian, data_out, ordering=ordering, fill_holes=fill_holes, bad_data=bad_data
    ; write output data
    file_out = medmap
    sxaddpar, xh,'COMMENT','----------------------------------'
    sxaddpar, xh,'HISTORY','median filtered map'
    sxaddpar, xh,'HISTORY','input file'+strtrim(file_in)
    sxaddpar, xh,'MFRADIUS',radius_deg,'[Deg] median filter radius'
    sxaddpar, xh,'COMMENT','----------------------------------'
    case type of
        2: begin
            case nmaps of
                1:    write_fits_map, file_out, data_out, xh
                else: write_tqu,      file_out, data_out, hdr= h, xhdr = xh
            endcase
        end
        3: begin
            pixel = where(data_out[*,0] ne !healpix.bad_value, obs_npix)
            print,obs_npix, minmax(data_out)
            if (cut4) then begin
                signal = data_out[pixel]
                write_fits_cut4, file_out, pixel, signal, hdr=h, xhdr=xh, nside=nside
            endif else begin
                tmp = data_out[pixel,0:ncol-1]
                write_fits_partial, file_out, pixel, tmp, hdr=h, xhdr=xh, nside=nside
            endelse
        end
    endcase
    

endelse

return
end
