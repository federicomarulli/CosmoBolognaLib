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
pro oplot_line_with_label, u, v, $
  linelabel=linelabel, putlabel=putlabel, flush=flush, lines=line_type, _extra = oplot_kw, charsize=charsize, $
  latex=latex, ltxstc=ltxstc, color=color
; , labrot=labrot, fixminus=fixminus

; note on oplot: linestyle in _extra prevails silently over explicit linestyle

lensegment = n_elements(u)
minseglen = 20
chars = keyword_set(charsize) ? charsize : 0
do_label =  (keyword_set(putlabel) && size(/tname, linelabel) eq 'STRING' && lensegment gt minseglen && chars gt 0.)
if (do_label) then begin
    middle = lensegment* (keyword_set(flush) ? 0.1 : 0.45)
    ; find out room to leave in line around label
    chsize = chars*!d.y_ch_size/(!d.y_vsize*!y.s[1]) ; size of character in units of data coordinates
    labelsize = chsize*strlen(linelabel)
    blank = long(labelsize/sqrt((u[middle-1]-u[middle+1])^2+(v[middle-1]-v[middle+1])^2))> 0.6
    ; make sure that label is shorther than line to be labelled
    if (2*blank lt lensegment) then begin
                                ; break line under label
        if ((middle-blank) gt 0) then oplot, u[0:middle-blank], v[0:middle-blank], _extra = oplot_kw, lines=line_type, color=color
        if ((middle+blank) lt (lensegment-1)) then oplot, u[middle+blank:*], v[middle+blank:*], _extra = oplot_kw, lines=line_type, color=color
                                ; put label
        step = 3
        dv = v[middle+step]-v[middle-step]
        du = u[middle+step]-u[middle-step]
        ;angle = atan(dv/du) * !radeg ; angle in degree
        angle = atan(dv,du) * !radeg ; angle in degree
        if (angle lt -89.9)    then angle += 180.
        if (angle gt  89.9)    then angle -= 180.
        ; if keyword_set(labrot) then angle += 180.
        ; offset label position to be level with line
        xlab = u[middle] + 0.3*chsize * cos((angle-90.)*!dtor)
        ylab = v[middle] + 0.3*chsize * sin((angle-90.)*!dtor)
        ; if keyword_set(fixminus) then linelabel = StrJoin(StrSplit(linelabel, '-', /Regex, /Extract, /Preserve_Null), '!M-!X')
        xyouts_latex, xlab, ylab, linelabel, align=0.5,orientation=angle, noclip=0,charsize=chars, $
                      latex=latex, ltxstc=ltxstc
    endif else begin
        ; if label too big, drop it and only plot line
        oplot, u, v, _extra = oplot_kw, lines=line_type, color=color
    endelse
endif else begin
    ; no label, line only
    oplot, u, v, _extra = oplot_kw, lines=line_type, color=color
endelse

return
end

;=============================================

pro oplot_sphere, u, v,  $
                  line_type=line_type, _extra = oplot_kw, linelabel=linelabel, flush=flush, charsize=charsize, $
                  latex=latex, ltxstc=ltxstc, color=color

if undefined(line_type) then line_type = 0
if undefined(color)     then color = [!p.color, !p.background]

; find points where line should be interrupted (large step in u)
bad = where(abs(u-shift(u,1)) gt .1, nbad)

iw = index_word(tag_names(oplot_kw),'PSYM',err=errword)
if errword eq 0 then begin
    if (oplot_kw.psym gt 0) then nbad = 0
endif

if (nbad eq 0) then begin
    if (line_type lt 0) then begin
        oplot, u, v, _extra = oplot_kw, col=color[1] ; white background
        oplot, u, v, _extra = oplot_kw, col=color[0], lines=abs(line_type)
    endif else begin
        oplot_line_with_label, u, v, linelabel=linelabel, _extra = oplot_kw, lines=line_type, $
          putlabel=1,flush=flush,charsize=charsize, $
          latex=latex, ltxstc=ltxstc, color=color[0]
    endelse
endif else begin
;    bad = [0,bad,n_elements(u)-1]
    bad = [0,bad,n_elements(u)]
    already = 0
    for j=0,nbad do begin
        lensegment = bad[j+1] - bad[j]
        if (lensegment gt 1) then begin
            u1 = u[bad[j]:bad[j+1]-1]
            v1 = v[bad[j]:bad[j+1]-1]
            if (line_type lt 0) then begin
                oplot, u1, v1, _extra = oplot_kw, col=color[1] ; white background
                oplot, u1, v1, _extra = oplot_kw, col=color[0], lines=abs(line_type)
            endif else begin
                putlabel = (~already)
                oplot_line_with_label, u1, v1, linelabel=linelabel, _extra = oplot_kw, lines=line_type, $
                  putlabel=putlabel, flush=flush, charsize=charsize, $
                  latex=latex, ltxstc=ltxstc, color=color[0]
                already = 1
            endelse        
        endif
    endfor
endelse




return
end

