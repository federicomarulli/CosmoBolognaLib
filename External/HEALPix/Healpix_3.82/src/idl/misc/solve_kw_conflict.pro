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
pro solve_kw_conflict, name1, name2, k1=k1, k2=k2, kout=kout, set=set, present=present, defined=defnd
;
; v1.0: ?2008
; v1.1: July 2008: replace 'defined' variable with 'defnd' to avoid conflict
; with 'defined' function
;

uname1 = strupcase(name1)
uname2 = strupcase(name2)
do_set = keyword_set(set)
do_present = keyword_set(present)
do_defined = keyword_set(defnd)

if (do_set+do_present+do_defined) ne 1 then begin
    message,'Choose either Set or Present or Defined'
endif

if (do_set) then begin
    ; Are keywords both set ?
    ks1 = keyword_set(k1) 
    ks2 = keyword_set(k2)
    if (ks1 && ks2) then begin
        message_patch,level=-1,'Set either '+uname1+' or '+uname2+' but not both!'
    endif
endif

if (do_present) then begin
    ; Are keywords both present ?
    ks1 = arg_present(k1) 
    ks2 = arg_present(k2)
    if (ks1 && ks2) then begin
        message_patch,level=-1,uname1+' or '+uname2+' should be present, but not both!'
    endif
endif

if (do_defined) then begin
    ; Are keywords both defined ?
    ks1 = defined(k1) 
    ks2 = defined(k2)
    if (ks1 && ks2) then begin
        message_patch,level=-1,'Define either '+uname1+' or '+uname2+' but not both!'
    endif
endif

if (ks1) then kout = k1
if (ks2) then kout = k2

return
end
