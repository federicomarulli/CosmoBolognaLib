;+
; NAME:
;       HPX_LATEXIFY
;
;PURPOSE
;	to use PSfrag to do latex syntax correctly
;SYNTAX
;	hpx_latexify, filename, tags, tex, [scale, ALIGNMENT=, HEIGHT=, HELP=, OUTNAME=, WIDTH=, KEEP_TEMPORARY_FILES=]
;INPUTS
;	filename: name of eps file
;	tags: name of dummy placeholder
;	tex: the tex you want to replace it
;	scale: the scale you want the latex to be [default to 1]
;NOTES:
;	requires the following latex packages properly installed:
;		geometry, graphicx, psfrag, color, helvet
;	requires you to use font=0
;
;	Follows from discussion by Saurav
;		http://www.astrobetter.com/idl-psfrag/
;	Writes and deletes files called $IDL_TMPDIR/psfrag_*_temp*
;EXAMPLE
;
;   SET_PLOT, 'PS'
;   !P.FONT=0
;   DEVICE, FILENAME='figure.eps'
;   plot,findgen(10),xtitle='xtitle'
;   DEVICE, /CLOSE
;
;   hpx_latexify, 'figure.eps', 'xtitle', '$M_{\odot}$'
;	
; HISTORY:
; adapted from latexify.pro by R. da Silva, UCSC, 1-4-11
;  and from LS_latexify.pro by L. Spencer (2011)
;  version 1.0: 2015-05-12
;  2017-05-18: added ALIGNMENT keyword
;  2017-06-13: added KEEP_TEMPORARY_FILES keyword
;	
;-
; *****************************
pro check_commands, commands
nc = n_elements(commands)
default_paths = getenv('PATH')
for k=0,nc-1 do begin
    command = commands[k]
    command_path = strtrim(file_which(default_paths, command),2)
    if (command_path eq '') then begin
        message,/info,level=-1,''
        message,/info,level=-1,'WARNING: the command  '+command+'  was not found in PATH ('+default_paths+')'
        message,/info,level=-1,''
    endif
endfor
return
end
; *****************************
pro hpx_latexify, filename, tag, tex, scale, $
                  alignment = alignment, $
                  height = height, $
                  help = help, $
                  outname=outname, $
                  width = width, $
                  keep_temporary_files = keep_tmp_files

routine = 'hpx_latexify'
syntax = routine+', filename, tags, tex, [scale, ALIGNMENT=, HEIGHT=, HELP=, OUTNAME=, WIDTH=, KEEP_TEMPORARY_FILES=]'

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if n_params() lt 3 then begin
    print,syntax
    return
endif

if ~keyword_set(scale) then scale=replicate(1., n_elements(tag))
scale=strtrim(scale*1.,2)

; temporary files
tmpdir = getenv('IDL_TMPDIR')
stime  = string(long(systime(1)*1e3 mod 1.e8), form='(i8.8)')
tmprad = add_prefix(tmpdir,path_sep(),/suff)+'psfrag_'+stime+'_temp'
noname=0
if ~keyword_set(outname) then begin
    outname=tmprad + '2.ps'
    cpname =tmprad + 'BKUP.ps'
    noname=1
endif

; generate LaTeX file
openw, lun, tmprad + '.tex', /get_lun
printf, lun, '\documentclass{article}'
printf, lun, '\nonstopmode' ; will crash on error instead of prompting user
printf, lun, '\usepackage{geometry, graphicx, psfrag}'
printf, lun, '\usepackage{color}'
; printf, lun, '\usepackage{helvet}'
; printf, lun, '\renewcommand{\familydefault}{\sfdefault}'
printf, lun,'\pagestyle{empty}'
printf, lun,'\geometry{paperwidth='+strtrim(width*1.01,2)+'cm,'+$
                    'paperheight='+strtrim(height*1.01,2)+'cm,margin=0pt}'
printf, lun,'\begin{document}'
for i=0, n_elements(tag)-1 do begin
    pospos = '[cc][cc]' ; centers of tag and tex will match
    if (n_elements(alignment) ge (i+1)) then begin
        if (alignment[i] lt 0.1) then pospos = '[cl][cl]' ; left  edges of tag and tex will match
        if (alignment[i] gt 0.9) then pospos = '[cr][cr]' ; right edges of tag and tex will match
    endif
    printf, lun,'\psfrag{'+tag[i]+'}'+pospos+'['+scale[i]+']{'+tex[i]+'}'
    ;;;;;print,      '\psfrag{'+tag[i]+'}'+pospos+'['+scale[i]+']{'+tex[i]+'}'
endfor
printf, lun,'\begin{center}'
printf, lun,'\includegraphics[width='+$
    	strtrim(width,2)+'cm,height='+strtrim(height,2)+'cm]{'+filename+'}'
printf, lun,'\end{center}'
printf, lun,'\end{document}'
close, lun
free_lun, lun

; generate PS file
check_commands, ['cd','latex']
spawn, 'cd '+tmpdir+'; latex '               +tmprad +'.tex', ls0
if ~file_test(tmprad+'.dvi') then begin
    print, ls0,form='(a)'
    message,'Error in LaTeX. Aborting.'
endif else begin
    check_commands, ['cd','dvips']
    command = 'cd '+tmpdir+'; dvips -o '+outname+' '+tmprad +'.dvi'
    if (strupcase(!version.os_family) eq 'WINDOWS') then begin
        spawn,      command, ls1
    endif else begin
        spawn, /sh, command+' > '+tmprad+'_dvi.log  2>&1 '
    endelse
; spawn, 'dvips -o '+tmprad+'.ps '+tmprad+'.dvi', ls1
; spawn, 'ps2epsi  '+tmprad+'.ps '+tmprad+'.epsi, ls2
; spawn, "perl -ne 'print unless /^%%BeginPreview/../^%%EndPreview/' <"+$
;	 tmprad + '.epsi > '+outname, ls3
endelse

;garbage cleanup
tmp_files =  tmprad+['2.ps','.tex','.aux','.log','.dvi','_dvi.log']
;;;;;;; spawn,'cp '+filename+' '+tmprad+'_in.ps'
if noname then begin
    check_commands, ['cp','mv']
    spawn, 'cp -f '+filename+' '+cpname
    tmp_files = [tmp_files, cpname]
    spawn, 'mv -f '+outname+' '+filename
    outname=filename
endif

if keyword_set(keep_tmp_files) then begin
    print,'Files NOT deleted: ',tmp_files
endif else begin
    file_delete, tmp_files, allow_nonexistent=1
endelse


end
