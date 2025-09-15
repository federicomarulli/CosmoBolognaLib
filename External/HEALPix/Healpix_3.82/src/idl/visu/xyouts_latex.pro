; ----------------------------------------------------------------------------

;+
; NAME:
;       XYOUTS_LATEX
;
; PURPOSE:
;       puts LaTeX string at arbitrary location
;
; CATEGORY:
;
; CALLING SEQUENCE:
;       XYOUTS_LATEX, X, Y, String, [ HELP=, LATEX=, LTXSTC=, +all XYOUTS keywords]
;
; INPUTS:
;       X, Y:  X and Y  position
;       String: character string to be put
;
; KEYWORD PARAMETERS:
;      LATEX: if set (to 1 or 2), the string should be interpreted as a LaTeX string
;         if set to 2 and in PS/PDF: genuine LaTeX
;         if set to 1 or in X/PNG/GIF/JPG:  TeXtoIDL emulation
;     
;      LTXSTC: structure containing information necessary for LaTeX handling by psfrag
;
; SIDE EFFECTS:
;      LTXSTC updated at each call (if LATEX is set)
;
; RESTRICTIONS:
;
; PROCEDURE:
;     invokes TeXtoIDL (for non-PS output) or 
;     prepares data for a subsequent call to LaTeX's psfrag (for PS or PDF output)
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;     2015-05-12: EFH, version 1.0
;-

; ----------------------------------------------------------------------------
function e2latex, instr, idl=idl
;
; turn C/fortran-like scientific notation into TeX/LaTeX or IDL power of 10
; eg:  5e3  -> 5$\ 10^{3}$
;           -> 5 10!u3!n
;

c_pos = ['e\+00', 'e\+01', 'e\+02', 'e\+03', 'e\+04', 'e\+05', 'e\+06', 'e\+07', 'e\+08', 'e\+09', $
         'e\+10' ,'e\+11', 'e\+12', 'e\+13', 'e\+14', 'e\+15']
t_pos = ['',     '10',   '10^{2}', '10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}','10^{9}',$
                   '10^{10}','10^{11}','10^{12}','10^{13}', '10^{14}', '10^{15}']
i_pos = ['',     '10',   '10!u2!n', '10!u3!n','10!u4!n','10!u5!n','10!u6!n','10!u7!n','10!u8!n','10!u9!n',$
                   '10!u10!n','10!u11!n','10!u12!n','10!u13!n', '10!u14!n', '10!u15!n']

c_neg = ['e-00', 'e-01', 'e-02', 'e-03', 'e-04', 'e-05',   'e-06',   'e-07',   'e-08',   'e-09', $
         'e-10' ,'e-11', 'e-12', 'e-13', 'e-14', 'e-15']
t_neg = ['',   '10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}','10^{-6}','10^{-7}','10^{-8}','10^{-9}',$
                   '10^{-10}','10^{-11}','10^{-12}','10^{-13}', '10^{-14}', '10^{-15}']
i_neg = ['',   '10!u-1!n','10!u-2!n','10!u-3!n','10!u-4!n','10!u-5!n','10!u-6!n','10!u-7!n','10!u-8!n','10!u-9!n',$
                   '10!u-10!n','10!u-11!n','10!u-12!n','10!u-13!n', '10!u-14!n', '10!u-15!n']

outstr = instr
if (keyword_set(idl)) then begin
    outstr = strtrans(outstr, c_pos, ' '+i_pos)
    outstr = strtrans(outstr, c_neg, ' '+i_neg)
endif else begin
    outstr = strtrans(outstr, c_pos, '$\ '+t_pos+'$')
    outstr = strtrans(outstr, c_neg, '$\ '+t_neg+'$')
endelse

return, outstr
end
; ----------------------------------------------------------------------------
function len_latex, instr
;
; estimate final length of LaTeX string
;
outstr = instr

; replace \command with single character
split = strsplit(outstr,'\\[A-Za-z]+',/regex, /extract, /preserve_null)
outstr = strjoin(split, 'x')

; remove $,{,},^,_
outstr = strtrans(outstr, ['\$','{','}','\^','_'], ['','','','',''])

; remove ' ' and '\' ?

;print,instr+' -> '+outstr
return, strlen(outstr)
end

; ----------------------------------------------------------------------------


pro xyouts_latex, x, y, str, $
                  alignment = alignment, $
                  charsize= charsize, $
                  color = color, $
                  font = font, $
                  help = help, $
                  latex=latex, $
                  ltxstc = ltxstc, $
                  width = width, $
                  _ref_extra = extra
;                   charthick = charthick, $
;                   text_axes = text_axes, $
;                   CLIP = clip, $
;                   DATA=data, DEVICE=device, $
;                   NOCLIP=noclip, NORMAL=normal, ORIENTATION=orientation, T3D=t3d, Z=z

routine = 'xyouts_latex'
syntax = routine+' x, y, string, HELP=, LATEX=, LTXSTC=, +all XYOUTS keywords'

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if n_params() ne 3 then begin
    print,syntax
    return
endif

do_latex = keyword_set(latex)
do_ps    = !d.name eq 'PS'

myfont = defined(font) ? font : !p.font
if (do_latex) then begin
    if (latex eq 2 && do_ps) then begin
                                ; for Postscript, create information to be used by LaTeX psfrag package
        if (size(/tname, ltxstc) ne 'STRUCT') then begin
            nmax = 100
            ltxstc={n:0, tag:replicate('void',nmax),  tex:replicate('void',nmax), $
                    scale:replicate(-1.,nmax), alignment:replicate(0.0,nmax)}
        endif 
        fudge_size = 1.0
        xxx = 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        tvlct,/get, rgb
        mycol = defined(color) ? color : !p.color
        rgb1 = strtrim(long(rgb[mycol,*]),2)
        scol = '\color[RGB]{'+rgb1[0]+','+rgb1[1]+','+rgb1[2]+'}'
        
        n = ltxstc.n
        tag1 = '#'+string(n,form='(i2.2)')+strmid(xxx,0, (len_latex(str)-3)>0)
        chars1 = defined(charsize) ? charsize*fudge_size : 1.*fudge_size
        ltxstc.tag[n]   = tag1
        ;ltxstc.tex[n]   = scol+e2latex(str)
        ltxstc.tex[n]   = scol+strtrans(e2latex(str),' ','{ }')
        ltxstc.scale[n] = chars1
        ltxstc.alignment[n] = defined(alignment) ? alignment : 0
        ltxstc.n        = n + 1
        myfont = 0  ; must be written in Device Font for psfrag to work
        mystr = tag1
    endif else begin
                                ; for non Postcript, use TeXtoIDL
        mystr = e2latex(textoidl(str),/idl)
        mystr = strtrans(mystr,['\$'],[''])
    endelse
endif else begin
    mystr = str
endelse

xyouts, x, y, mystr, $
        alignment = alignment, $
        charsize= charsize, $
        color = color, $
        font = myfont, $
        width = width, $
        _strict_extra = extra
;         charthick = charthick, $
;         text_axes = text_axes, $
;         CLIP = clip, $
;         DATA=data, DEVICE=device, $
;         NOCLIP=noclip, NORMAL=normal, ORIENTATION=orientation, T3D=t3d, Z=z

return
end

