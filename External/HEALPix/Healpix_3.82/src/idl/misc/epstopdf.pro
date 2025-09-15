;+
; NAME:
;       epstopdf
;
; PURPOSE:
;       PostScript to PDF conversion with a tight BoundingBox
;       Syntax is very close to that of cgPS2PDF
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;       epstopdf, ps_file [, pdf_file], [DELETE_PS=, HELP=, OPTIONS=, SHOWCMD=, SILENT=,
;       SUCCESS=]
;
; INPUTS:
;      ps_file:  input PostScript file
;
; OPTIONAL INPUTS:
;      pdf_file:  output PDF file
;       if absent, it will be generated from ps_file by replacing
;       (.ps or .eps) suffix with .pdf
;
;
; KEYWORD PARAMETERS:
;      DELETE_PS: if set, the PostScript file will be deleted
;
;      HELP: is set, this documentation header is printed out
;        and the routine leaves
;
;      OPTIONS: string containing options to be passed to epstopdf.pl
;
;      SHOWCMD: if set, the routine prints the epstopdf.pl debug info
;
;      SILENT: if set, the routine runs silently
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       SUCCESS: on output, will be 1 in case of success, 0 otherwise
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;      calls external script
;      creates and delete files
;
; RESTRICTIONS:
;
; PROCEDURE:
;     calls public domain epstopdf.pl script (which itself calls 'gs') to
;     convert a PS file into PDF.
;     If that fails, D.W. Fanning's (Coyote) cgps2pdf is called instead.
;
; EXAMPLE:
;     epstopdf, 'file.ps', 'file.pdf'
;
;
; MODIFICATION HISTORY:
;   2015-03-17:  E. Hivon, version 1.0
;
;-
pro epstopdf, ps_file, pdf_file, $
              delete_ps = delete_ps, $
              help = help, $
              options = options, $
              showcmd = showcmd, $
              silent = silent, $
              success = success


routine = 'epstopdf'
syntax = routine+', ps_file [,pdf_file], [delete_ps=, options=, help=, showcmd=, silent=, success=]'

if keyword_set(help) then begin
    doc_library,routine
    return
endif

np = n_params()
if (np ne 1 and np ne 2) then begin
    print,syntax
    return
endif

if ~file_test(ps_file[0]) then begin
    message,'Input PS file '+ps_file[0]+' not found.'
    return
endif

if (np eq 1) then begin
    root_name = cgRootName(ps_file[0], DIRECTORY=thisDir)
    pdf_file  = Filepath(ROOT_DIR=thisDir, root_name + ".pdf")
endif

if (pdf_file[0] eq ps_file[0]) then begin
    message,'PS and PDF files must have different names'
    return
endif

; first try epstopdf.pl which produces PDF files with BoundingBox
bin = 'epstopdf.pl'
bin_path = file_which(bin)
defsysv, '!healpix',          exists = exists
if ~exists then init_healpix
exists2 = 0
junk1 = where(tag_names(!healpix) eq 'PATH', exists1)
if (exists1 eq 1) then junk2 = where(tag_names(!healpix.path) eq 'SRC', exists2)
if (bin_path eq '' && exists2) then begin
    bin_path = file_search(!healpix.path.src, bin, count=count)
endif

success = 0
myoptions  = keyword_set(showcmd)             ? ' --debug '        : ' '
myoptions += size(/tname,options) eq 'STRING' ? ' '+options[0]+' ' : ' '
;myoptions += " --gsopts=-c <</Orientation 1>> setpagedevice " ;-c <</Orientation 1>> setpagedevice "
;myoptions += " --gsopts=--debug --debug " ;-c <</Orientation 1>> setpagedevice "
if (bin_path ne '') then begin
    command =  bin_path[0]+' '+ps_file[0]+myoptions+' --outfile='+pdf_file[0]
    file_delete, pdf_file[0], /quiet
    spawn, /sh, command, out
    if (file_test(pdf_file[0]) && keyword_set(delete_ps)) then file_delete,ps_file[0]
    success = file_test(pdf_file[0])
endif

if (success && ~keyword_set(silent)) then begin
    print,ps_file[0]+' -> '+pdf_file[0]
endif

; in case of failure, use cgps2dpf which produces PDF files
; at paper size (letter, a4, ...)
if (~success) then begin
    message,bin+' not available or not working! Trying cgps2pdf instead.',/info
    cgps2pdf, ps_file, pdf_file, delete_ps = delete_ps, silent = silent, success = success
endif


return
end

