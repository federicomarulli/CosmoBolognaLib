pro hpx2dm, file_in, select_in, $
                     _extra = other_kwds, $
                     gif = gif, $ ; not allowed
                     help = help, $
                     preview = preview, $
                     pxsize = pxsize
;+
; NAME:
;      HPX2DM
;
; PURPOSE:
;    Turns a Healpix data set into DomeMaster compliant image
;   (azimuthal equidistant projection of the half-sphere)
;   in a PNG or lossless JPEG file.
;   See eg http://fulldome.ryanwyatt.net/fulldome_domemasterSpec_v05.pdf
;
; CATEGORY:
;
; CALLING SEQUENCE:
;   Hpx2dm, file [,select, COORD=, JPEG=, PNG=, PREVIEW=, PXSIZE=, ROT= ...
;   + all remaining AZEQVIEW's keywords (except CROP, HALF_SKY, PNG, RESO_ARCMIN, WINDOW)]
;
;
; INPUTS:
; 	File = 
;          by default,           name of a FITS file containing 
;               the healpix map in an extension or in the image field
;          or    name of a variable containing the healpix map
;          if Save is set   :    name of an IDL saveset file containing
;               the healpix map stored under the variable  data
;
; OPTIONAL INPUTS:
;       Select =  if the file read is a multi column BIN table, Select indicates
;                 which column is to be plotted (the default is to plot the
;                 first column containing some signal, as opposed to pixel index)
;               can be either a name : value given in TTYPEi of the FITS file
;                        NOT case sensitive and
;                        can be truncated, 
;                        (only letters, digits and underscore are valid)
;               or an integer        : number i of the column
;                            containing the data, starting with 1
;
; KEYWORD PARAMETERS:
;      COORD=
;      JPEG=
;      PNG=
;      PREVIEW=
;      PXSIZE=  number of pixels in each dimension. default=1024
;      ROT=
;
;
; OUTPUTS:
;       none
;
; OPTIONAL OUTPUTS:
;       none
;
; COMMON BLOCKS:
;       none
;
; SIDE EFFECTS:
;       creates a PNG or JPEG file and attempts to visualize it
;
; RESTRICTIONS:
;       none?
;
; PROCEDURE:
;       calls AzEqView
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;      2012-01-19: creation
;-

routine = 'HPX2DM'
syntax = [routine+', File [,select, COORD=, JPEG=, PNG=, PREVIEW=, PXSIZE=, ROT= ...,',$
'             + most of AZEQVIEW''s keywords (except CROP, HALF_SKY, PNG, RESO_ARCMIN, WINDOW)]']

if (keyword_set(help)) then begin
    doc_library,routine
    return
endif

if (n_params() eq 0) then begin
    print,syntax,format='(a)'
    print,'Type: '+routine+', /help'
    print,'   for more help.'
    return
endif

if keyword_set(gif) then begin
    print,syntax,format='(a)'
    message,' GIF output not allowed.'
endif
;-----------------

; forced : CROP, HALF_SKY, RESO_ARCMIN, WINDOW
mypxsize = defined(pxsize) ? pxsize : 1024
reso_arcmin = 180.d0 * 60.d0 / mypxsize
azeqview, file_in, select_in, $
          crop = 1, $
          half_sky = 1, $
          preview = preview, $
          pxsize = mypxsize, $
          reso_arcmin = reso_arcmin, $
          _strict_extra = other_kwds, $
          window = -1

return
end

