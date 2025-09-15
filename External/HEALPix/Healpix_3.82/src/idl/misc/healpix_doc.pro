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
pro healpix_doc, help=help, epub=epub, html=html, pdf=pdf, whole=whole
;+
; NAME:
;         Healpix_doc
;
;
; PURPOSE:
;         displays epub, html or pdf IDL-Healpix documentation
;
;
; CALLING SEQUENCE:
;         healpix_doc, [EPUB= | HTML= | PDF=] [, HELP=,  WHOLE=]
;
;
; KEYWORD PARAMETERS:
;
;         HELP= shows this header and exits
;
;         EPUB=   if set, the EPUB documentation is shown using
;            the e-viewer found by the $HEALPIX/healpix_doc script
;
;
;         HTML=  if set, the HTML documentation is shown using
;            the web browser found by the $HEALPIX/healpix_doc script
;            If the browser is already in use, a new tab is open.
;
;
;         PDF=   if set, the PDF documentation is shown using
;            the pdf viewer found by the $HEALPIX/healpix_doc script
;
;            Either EPUB or HTML or PDF must be set.
;
;
;         WHOLE=  if set, the whole Healpix documentation is accessible,
;              not just the IDL related part.
;
; SIDE EFFECTS:
;      
;       will start (or open a new tab in) an epub reader, pdf viewer or html browser
;
;
; PROCEDURE:
;        
;      calls healpix_doc script
;
; EXAMPLE:
;
;       healpix_doc, /html
;
;
; MODIFICATION HISTORY:
;        created 2009-10-07 by EH.
;        2010-03-12: opens main.htm instead of index.htm when used with
;        /html,/whole
;        2012-03-01: replaces online_help with call to healpix_doc script
;        2013-03-28: recursive search of healpix_doc script
;        2020-06-09: added support for epub
;-

in_gdl = is_gdl()
in_fl  = is_fl()
in_idl = ~(in_gdl || in_fl)

routine = 'healpix_doc'
syntax = routine+', [HELP=, HTML=, PDF=, EPUB=, WHOLE=]'
do_whole = keyword_set(whole)
do_pdf = keyword_set(pdf)
do_html = keyword_set(html)
do_epub = keyword_set(epub)

if keyword_set(help) then begin
    doc_library,routine
    return
endif

if (n_params() ne 0 || (do_pdf+do_html+do_epub) ne 1) then begin
    print,syntax
    print, ' '
    print,'Choose either /epub, /html or /pdf'
    print, ' '
    print,'Type '+routine+',/help'
    print,'   for more help.'
    return
endif


init_healpix

found = where( tag_names(!healpix.path) eq 'DOC', nfound)
nfound2 = (nfound eq 1) ? N_tags(!healpix.path.doc) : 0

sep = path_sep()
if (nfound2 ge 2) then begin
    ; use !healpix.path.doc.(*) if they are defined
    path_pdf  = !healpix.path.doc.pdf
    path_html = !healpix.path.doc.html
    if (nfound2 eq 3) then begin
        path_epub = !healpix.path.doc.epub
    endif else begin
        path_epub  = !healpix.directory+sep+'doc'+sep+ 'epub'+sep
    endelse
endif else begin
    ; otherwise, build path with !healpix.directory
    path_main  = !healpix.directory+sep+'doc'+sep
    path_pdf   = path_main + 'pdf' +sep
    path_html  = path_main + 'html'+sep
    path_epub  = path_main + 'epub'+sep
endelse

if do_pdf then begin
    if (do_whole) then begin
        doc_path = path_pdf + 'pdf_index.pdf'
    endif else begin
        doc_path = path_pdf + 'idl.pdf'
    endelse
endif

if do_html then begin
    if (do_whole) then begin
;;        doc_path = path_html + 'index.htm'
        doc_path = path_html + 'main.htm'  ; corrected 2010-03-12
    endif else begin
        doc_path = path_html + 'idl.htm'
    endelse
endif

if do_epub then begin
    doc_path = path_epub + 'HEALPixDocumentation.epub'
    if (do_whole) then begin
        message,/info,'Warning: \Whole ignored for EPUB'
    endif
endif


if (~file_test(doc_path)) then begin
    print,'ERROR: documentation file ('+doc_path+') not found.'
    return
endif else begin
    ; cmd = !healpix.directory+'/'+'healpix_doc'
    ; recursively search subdirectories of !healpix.directory
    cmd = file_search(!healpix.directory,'healpix_doc',/FULLY_QUALIFY_PATH,count=count)
    if (count lt 1) then begin
        print,'ERROR: healpix_doc shell script was not found.'
        return
    endif 
    cmd = cmd[0]
    if do_pdf  then cmd += ' -p '
    if do_html then cmd += ' -h '
    if do_epub then cmd += ' -e '
    cmd += doc_path
    spawn, cmd, /sh

;     if (in_idl) then begin
;         online_help, /full_path, book=doc_path
;     endif else begin
;         bg = ' & '
;         if keyword_set(pdf) then begin
;             cmd = !dir+'/bin/online_help_pdf  '+doc_path+ bg
;         endif else begin
;             cmd = !dir+'/bin/online_help_html '+doc_path+ bg
;         endelse
;         spawn,cmd
;     endelse
endelse


return
end




