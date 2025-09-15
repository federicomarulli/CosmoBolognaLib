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
function index_word, word_list, word, ERROR = error, VALUE = value
;+
; NAME:
;       INDEX_WORD
;
; PURPOSE:
;       returns the index of the first occurence of Word in the string
;       array word_list
;       (both strings are set to uppercase before comparison)
;
; CALLING SEQUENCE:
;       result = index_word(word_list, word[, ERROR=error, VALUE = value])
;
; INPUTS:
;      word_list = string array
;
;      word      = string looked for in Word_list   or integer
;
; OPTIONAL OUTPUTS:
;      error = named variable containing the error code
;            0 : no error
;            1 : invalid arguments
;            2 : Word out of range (if word is an integer)
;            3 : Word not found (it Word is a string)
;
;      value = named variable containing the value of the string(word)
;
; MODIFICATION HISTORY:
;      March 1999, Eric Hivon, Caltech
;      Sept  2000, bugs correction
;-


  if N_params() LT 2 then begin
      print,'error in index_word : INDEX_WORD, Word_list, Word, [ERR]'
      error = 1
      return, -1L
  endif

  if datatype(word_list,2) ne 7 then begin
      print,'error in index_word : Word_list should be a name'
      error = 1
      return,-1L
  endif

  n_words = N_ELEMENTS(word_list)
  if datatype(word,2) le 3 then begin  ; integer
      if (word[0] LT 0 or word[0] GT n_words-1) then begin
          error = 2
          return, -1L
      endif
      error = 0
      value = STRTRIM(word_list(word[0]),2)
      return,LONG(word)
  endif

  if datatype(word,2) eq 7 then begin  ; string
      nw = n_elements(word)
      for i=0,nw-1 do begin
          template = STRTRIM( STRUPCASE(word[i]), 2)
          len  = STRLEN(template)
          list = STRTRIM(STRUPCASE(word_list),2)
          list = STRMID(list,0,len)
          index = WHERE( list EQ template, ni)
          if (ni gt 0) then goto, found
      endfor
      found:
      if (ni LT 1) then begin
          error = 3
          return,-1
      endif else begin
          error = 0
          value = STRTRIM(word_list(index(0)),2)
          return,LONG(index(0))
      endelse
  endif


  error = 1
  return,-1   ; Word is not string nor an integer
end



