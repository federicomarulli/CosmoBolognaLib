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
function sample_sparse_array, in_values, out_pix, id2, id3, in_pix=in_pix,  default=default, help=help
;+
; NAME:
;   sample_sparse_array
;
; PURPOSE:
;   returns value of array defined on a sparse set of indices, with a minimal
;   memory footprint:
;   out_values = sample_sparse_array(in_values, out_pix, i2, in_pix=in_pix, default=default)
;   is equivalent to
;     tmp_array = replicate(default, npix)  ; fill (big) array with default value
;     tmp_array[in_pix] = in_values[*,i2]   ; set subset of points to arbitrary values
;     out_values = tmp_array[out_pix]       ; read arbitrary set of points
;   without creating tmp_array which could be prohibitive
;
; CATEGORY:
;
; CALLING SEQUENCE:
;    out_values = sample_sparse_array(in_values, out_pix [, id2, id3, in_pix=, default=])
;
; INPUTS:
;    in_values : input array values
;    out_pix : index at which we want array values (integer >=0)
;
; OPTIONAL INPUTS:
;    id2: 2nd dimension index of in_values(default: 0)
;    id3: 3nd dimension index of in_values (default: 0)
;
; KEYWORD PARAMETERS:
;    in_pix: index at which input array is defined
;      - if in_pix is present, it must have the same size as in_values
;      - if in_pix is absent the input array is assumed to be defined for all
;      positive index values and the output is simply
;      out_values = in_values[out_pix]
;
;    default: default values of input array for indices outside of in_pix
;
;    help: prints out this help header
;
; OUTPUTS:
;    out_values : values of in_values where out_pix = in_pix, or default otherwise
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;           None
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;         finds common range of indices, and use histogram
;
; EXAMPLE:
;         p = [2, 30, 15]  &   x = p * 10.
;         print, sample_sparse_array(x, in=p, [15,14,16], default=-1)
;         will return
;              150.00   -1.00  -1.00
;
; MODIFICATION HISTORY:
;   2008-04-08: EH, IAP, first version
;-

routine = 'sample_sparse_array'
syntax = routine+', in_values, out_pix [, id2, id3, IN_PIX=,  DEFAULT=, HELP=]'

if (keyword_set(help)) then begin
    doc_library, routine
    return,-1
endif

nn = n_params()
if (nn lt 2) then begin
    print,syntax
    return,-1
endif

n_out = n_elements(out_pix)

iid2 = defined(id2)? id2 : 0
iid3 = defined(id3)? id3 : 0
;partial = (n_elements(in_pix) eq n_elements(in_values[*,iid2,iid3]))
partial = (n_elements(in_pix) eq (size(in_values))[1])

if (n_elements(in_pix) gt 0 && ~partial) then begin
    message,syntax,/info
    message,'In_values and In_pix should the same number of elements',/info
    return,-1
endif

mydefault = (keyword_set(default)) ? default : 0
mydefault = mydefault + 0*in_values[0] ; same type as input values

if ~partial then begin
    out_values = in_values[out_pix, iid2, iid3]
    return, out_values
endif else begin
    n_in = n_elements(in_pix)
    min_in = min(in_pix, max=max_in)
    min_out = min(out_pix, max=max_out)

    low = min_in > min_out ; highest of minima
    hi  = max_in < max_out ; lowest of maxima

    if (low gt hi) then begin
        ; no match
        out_values = replicate(mydefault, n_out)
        return, out_values
    endif else begin
        out_values = replicate(mydefault, n_out)
        np = hi - low + 1
        onegulp = long(1.e6)
        if (np le onegulp) then begin

            g_in  = where( in_pix ge low and in_pix le hi,  ng_in)
            g_out = where(out_pix ge low and out_pix le hi, ng_out)

            tmp = replicate(mydefault, np)

            if (ng_in gt 0) then tmp[in_pix[g_in]-low] = in_values[g_in, iid2, iid3]
            if (ng_out gt 0) then out_values[g_out] = tmp[out_pix[g_out] - low]
            return, out_values
        endif else begin
            nslices = (np + onegulp -1)/ onegulp
            bs = round(np*1.d0/nslices) > 1
            h_in  = histogram(in_pix,  min= low, max = hi, rev=r_in, binsize=bs, locations=low_in)
            h_out = histogram(out_pix, min= low, max = hi, rev=r_out, binsize=bs, locations=low_out)
            tmp = replicate(mydefault, bs+1)
            for i=0L, n_elements(h_in)-1 do begin
                if (h_in[i] gt 0) then begin
                    list1 = r_in[r_in[i]:r_in[i+1]-1]
                    tmp[in_pix[list1] - round(low_in[i]) ] = in_values[ list1, iid2, iid3 ]
                endif
                if (h_out[i] gt 0) then begin
                    list2 = r_out[r_out[i]:r_out[i+1]-1]
                    out_values[ list2 ] = tmp[out_pix[list2] - round(low_out[i]) ]
                endif
                tmp[*] = mydefault
            endfor
            return, out_values
        endelse
    endelse

endelse


return, out_values
end
