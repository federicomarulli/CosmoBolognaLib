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
pro mollview, file_in, select_in, $
              ASINH = asinh, $
              BAD_COLOR = bad_color, $
              BG_COLOR = bg_color, $
              CHARSIZE = charsize, $
              CHARTHICK = charthick, $
              COLT = colt, $
              COORD = coord, $
              CROP = crop, $
              CUSTOMIZE = customize, $
              DEFAULT_SETTINGS = default_settings, $
              EXECUTE = execute, $
              FACTOR = factor, $
              FITS = fits, $
              FG_COLOR=fg_color, $
              FLIP = flip, $
              GAL_CUT = gal_cut, $
              GIF = gif, $
              GLSIZE = glsize, $
              GRATICULE = graticule, $
              HELP = help, $
              HBOUND = hbound, $
              HIST_EQUAL = hist_equal, $
              HXSIZE = hxsize, $
              IGLSIZE = iglsize, $
              IGRATICULE=igraticule, $
              JPEG=jpeg, $
              LATEX=latex, $
              LOG = log, $
              MAP_OUT = map_out, $
              MAX = max_set, $
              MIN = min_set, $
              NESTED = nested_online, $
              NOBAR = nobar, $
              NOLABELS = nolabels, $
              NO_DIPOLE = no_dipole, $
              NO_MONOPOLE = no_monopole, $
              OFFSET = offset, $
              ONLINE = online, $
              OUTLINE = outline, $
              PDF = pdf, $
              PFONTS = pfonts, $
              PNG = png, $
              POLARIZATION = polarization, $
              PREVIEW = preview, $
              PS = ps, $
              PXSIZE = pxsize, $
              QUADCUBE = quadcube, $
              RETAIN = retain, $
              ROT = rot, $
              SAVE = save, $
              SILENT = silent, $
              SILHOUETTE = silhouette, $
              SUBTITLE = subtitle, $
              TITLEPLOT = titleplot, $
              TRANSPARENT = transparent, $
              TRUECOLORS = truecolors, $
              UNITS = units, $
              WINDOW = window, $
              XPOS = xpos, $
              YPOS = ypos

;+
; NAME:
;	MOLLVIEW, GNOMVIEW, CARTVIEW, ORTHVIEW
;
; PURPOSE:
; 	tools to view a Mollweide/gnomic/cartesian/orthographic projection of maps binned
;	in Healpix or COBE Quad-Cube pixelisation
;
; CALLING SEQUENCE:
; 	xxxxVIEW, File, [Select, ] $
;                       [ASINH=, $
;                       BAD_COLOR=, BG_COLOR=, $
;                       CHARSIZE=, COLT=, COORD=, CROP=, $
;                       CUSTOMIZE=, $
;                       DEFAULT_SETTINGS=, $
;                       EXECUTE=, $
;                       FACTOR=, FG_COLOR=, FITS=, FLIP=, $
;                       GAL_CUT=, GIF=, GLSIZE=, GRATICULE=, $
;                       HALF_SKY =, HBOUND =, HELP =, HIST_EQUAL=, HXSIZE=, $
;                       IGLSIZE=, IGRATICULE=, $
;                       JPEG=, $
;                       LATEX=, LOG=, $
;                       MAP_OUT=, MAX=, MIN=, $ 
;                       NESTED=, NOBAR=, NOLABELS=, NOPOSITION=, $
;                       NO_DIPOLE=, NO_MONOPOLE=, $
;                       OFFSET =, ONLINE=, OUTLINE=, $
;                       PFONTS=, PNG=, POLARIZATION=, PREVIEW=,$
;                       PS=, PXSIZE=, PYSIZE=, $
;                       QUADCUBE= , $
;                       RESO_ARCMIN= , ROT=, $
;                       SAVE=, SHADED=, SILENT=, SILHOUETTE=, STAGGER=, SUBTITLE=, $
;                       TITLEPLOT=, TRANSPARENT=, TRUECOLORS= $
;                       UNITS=, WINDOW=, XPOS=, YPOS=]
;                        
;  all the arguments and parameters are identical for all the
;  routines, excepted stated otherwise.
;
;
; INPUTS:
; 	File = 
;          by default,           name of a FITS file containing 
;               the healpix map in an extension or in the image field
;          if Online is set :    name of a variable containing
;               the healpix map
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
;               (see the Examples below)
;
; OPTIONAL INPUT KEYWORDS:
;
;       ASINH: if set, the color table is altered in to emulate the effect of replacing
;            the data by sinh^{-1}(data) in order to enhance the low contrast regions.
;            Can be used in conjonction with FACTOR and OFFSET, but can not be
;            used with /LOG nor /HIST_EQUAL
;
;       BAD_COLOR: color given to missing pixels (having
;            !healpix.bad_value=-1.6375e30 or NaN value on input).
;          The color can be provided as either:
;           - a single integer in [0,255], specifying the index to be used in
;             the color table chosen via COLT (in which the indexes 0, 1 and 2
;             are reserved for black, white and grey respectively),
;           - a 3 element vector, with each element in [0,255], specifying the
;             amount of RED, GREEN and BLUE
;           - a 7-character string, starting with '#', specifying the color in
;             HTML Hexadecimal fashion (eg, '#ff0000' for red).
;          Default: neutral grey (=2, =[175, 175, 175], ='#afafaf')
;          See also: BG_COLOR, FG_COLOR, TRANSPARENT
;
;       BG_COLOR: color given to pixels outside the sphere
;          See BAD_COLOR for expected format
;          Default: white (=1, =[255, 255, 255], ='#ffffff')
;          See also: BAD_COLOR, FG_COLOR, TRANSPARENT
;
;       CHARSIZE : overall multiplicative factor applied to the size of all
;               characters appearing on the plot
;                default = 1.0
;
;       CHARTHICK : character thickness (in TITLE, SUBTITLE and color bar labeling).  
;               Other characters thickness (such as graticule labels), can be 
;               controlled with !P.CHARTHICK.
;                default = 1
;
; 	COLT : color table index:
;              -Indexes [0,40] are reserved for standard IDL color tables, while
;               [41,255] are used for user defined color tables read from disc (created and
;               written to disc with MODIFYCT), if any.
;              -Indexes 1001 (or 'planck1', case insensitive) and 1002 (or 'planck2')
;               are reserved for Planck color tables 1 and 2 (see Planck_colors routine)
;              -If the index does not match any existing table, or if it is
;              above 255, the current
;               table (modifiable with TVLCT, XLOADCT, XPALETTE, ... 
;               or eg, J.Davenport's cubehelix.pro implementation of D. Green cubehelix
;               color scheme) is used instead.
;              -If not set, the color table will be 33 (Blue-Red).
;              -If colt<0, the IDL color table ABS(colt) is used, but the scale is
;              reversed (ie a red to blue scale becomes a blue to red scale).
;              (note: -0.1 can be used as negative 0)
;          See also: BAD_COLOR, BG_COLOR, FG_COLOR
;
;       COORD : vector with 1 or 2 elements describing the coordinate system of the map 
;                either 'C' or 'Q' : Celestial2000 = eQuatorial,
;                       'E'        : Ecliptic,
;                       'G'        : Galactic 
;               if coord = ['x','y'] the map is rotated from system 'x' to system 'y'
;               if coord = ['y'] the map is rotated to coordinate system 'y' (with the
;               original system assumed to be Galactic unless indicated otherwise in the file)
;                  see also : Rot
;
;       CROP : if set, the image produced (gif, jpeg, pdf, png, ps, X) only
;               contains the projected map and not the title, color bar, ...
;               (see also : GIF, JPEG, PDF, PNG, PS)
;
;       CUSTOMIZE: (input) IDL structure containing customization parameters
;          whose default values are listed in DEFAULT_SETTINGS
;           accept   aspos.x, aspos.y:               astronomical coordinates,
;                    cbar.dx, cbar.dy:               color bar X and Y sizes,
;                    cbar.spaces, cbar.ty:           management of text around color bar
;                    cbar.box:                       thickness of box around color bar
;                    cring.dx, cring.xll, cring.yll: color ring size and location
;                    pdf.debug:                      if set to 1, and SILENT not set, will help debugging PDF generation
;                    subtitle.x, subtitle.y, subtible.charsize: subtitle X,Y location and charsize
;                    title.x, title.y, title.charsize:        title X and Y location
;                    vscale.x, vscale.y:                X,Y location for polarization vector legend
;
;         eg   customize={cbar:{dx:1./3., dy:1./70.}, title:{x:0.5, y:0.95}, subtitle:{x:0.5, y:0.905}}
;
;       DEFAULT_SETTINGS: (output) IDL structure listing the default values
;          of the plot settings, like the color bar size, title and subtitle location
;          can be inspected with HELP_ST
;
;       EXECUTE: character string containing an IDL command to be executed in
;                the plotting window
;
;       FACTOR : multiplicative factor to be applied to the data (default = 1.0)
;               the data plotted is of the form FACTOR*(data + OFFSET)
;               see also : OFFSET, LOG
;
;       FG_COLOR: color of title and subtile characters, 
;          graticule lines and labels, units, outlines ...
;          See BAD_COLOR for expected format
;          Default: black (=0, =[0, 0, 0], ='#000000')
;          See also: BAD_COLOR, BG_COLOR
;
;       FITS : string containing the name of an output fits file with
;       the projected map in the primary image
;	      if set to 0 or not set : no .FITS done
;	      if set to 1            : output the plot in plot_XXX.fits
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;	      if set to a file name  : output the plot in that file 
;	   * For compatibility with standard FITS viewers (including STIFF), 
;        unobserved pixels, and pixels outside the sphere, take the value {\tt
;        NaN} (ie {\tt !values.f\_nan} in IDL).
;          * The resulting FITS file can be read in IDL with eg. map=readfits(filename). 
;          * In the case of orthographic projection, HALF_SKY must be set.
;
;       FLIP : if set, the longitude increases to the right, whereas by
;               default (astro convention) it increases towards the left
;
;       GAL_CUT: (positive float) specifies the symmetric galactic cut in degree
;             outside of which the the monopole and/or dipole fitting is done
;             (see also: NO_DIPOLE, NO_MONOPOLE)
;             ** mollview and orthview only **
;
;	GIF : string containing the name of a .GIF output
;	      if set to 0 or not set : no .GIF done
;	      if set to 1            : output the plot in plot_XXX.gif
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;	      if set to a file name  : output the plot in that file 
;             (see also : CROP, JPEG, PNG, PDF, PS and PREVIEW)
;
;       GLSIZE : character size of the graticule labels in units of CHARSIZE
;             can be a scalar (which applies to both parallel and meridian labels),
;             or a 2 element vector (interpreted as [meridian_label_size, parallel_label_size])
;             default: 0 (ie, no labeling of graticules)
;             (see also: CHARSIZE, GRATICULE)
;
; 	GRATICULE : if set, puts a graticule with delta_long = delta_lat = default
;         if graticule is set to a scalar > gmin delta_long = delta_lat = graticule
;         if set to [x,y] with x,y > gmin then delta_long = x and delta_let = y
;         ** cartview : default =  5, gmin =  0 **
;         ** gnomview : default =  5, gmin =  0 **
;         ** mollview : default = 45, gmin = 10 **
;         ** orthview : default = 45, gmin = 10 **
;         The graticule lines thickness is controlled with !P.THICK
;
;       HALF_SKY: if set, only shows only one half of the sky 
;          (centered on (0,0) or on the location parametrized by Rot) instead of the full sky
;             ** orthview only **
;        
;       HBOUND: scalar or vector of up to 3 elements.
;          For Hbound[i]>0, overplot the boundaries of Healpix pixels
;           for the resolution parameter Nside=hbound[i].
;           The first Nside will be plotted with solid lines, 
;           the second one (if any) with dashes, 
;           the third one (if any) with dots. Obviously, better results are
;           obtained for Hbounds elements in growing order.
;           Since 0-valued boundaries are not plotted, but used for linestyle
;           assignment, providing Hbound=[0,4] (or [0,0,4]) will
;           plot Nside=4 boundaries with dashes (resp. dots), while Hbound=4 would plot the same
;           boundaries with solid lines.
;
;       HELP : if set, the routine header is printed (by doc_library)
;             and nothing else is done
;
;	HIST_EQUAL : if not set, uses linear color mapping and 
;                     		puts the level 0 in the middle
;                     		of the color scale (ie, green for Blue-Red)
;				unless MIN and MAX are not symmetric
;	       	      if set,     uses a histogram equalized color mapping
;			(useful for non gaussian data field)
;                     (see also : LOG)
;
; 	HXSIZE: horizontal dimension (in cm) of the Hardcopy plot : Only for postscript printout
;    		** mollview : default = 26 cm ~ 10 in **
;               ** mollview : default = 15 cm         **
;    		(useful for large color printer)
;               (see also : PXSIZE)
;
;       IGLSIZE : character size of the input coordinates graticule labels in units of CHARSIZE
;             either scalar or 2-element vector (see GLSIZE)
;             default: 0 (ie, no labeling of graticules)
;             (see also: CHARSIZE, IGRATICULE)
;
;       IGRATICULE: if set, puts a graticule in the input coordinates
;          if both graticule and igraticule are set, these ones will
;          be represented by dashes
;
;	JPEG : string containing the name of a (lossless) .JPEG output
;	      if set to 0 or not set : no .JPEG done
;	      if set to 1            : output the plot in plot_XXX.jpeg
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;	      if set to a file name  : output the plot in that file 
;             (see also : CROP, GIF, PDF, PNG, PS and PREVIEW)
;
;       LATEX: if set to 1 or 2, enables LaTeX handling of character strings 
;         such as TITLEPLOT, SUBTITLE and UNITS
;        - if set to 2, together with PS or PDF outputs, these strings (and the graticule labels) 
;          will be processed by genuine LaTeX and inserted in the final PS or PDF file 
;         using psfrag package
;          (requires LaTeX and its color, geometry, graphicx and psfrag packages)
;           PFONTS settings will be ignored
;        - if set to 1, whatever the output is, LaTeX is partially emulated with TeXtoIDL routines
;          shipped with HEALPix (no extra requirements)
;           Can be used with PFONTS
;         default: 0 (not set)
;
; 	LOG: display the log of map (see also : HIST)
;         only applies to pixel with signal > 0.
;         (see OFFSET to offset signal)
;
;       MAP_OUT : name of the IDL variable that will contain
;         an un-altered projected map.
;         Unobserved pixels, and pixels outside the sphere take 
;       value !healpix.bad_value=-1.6375e30
;
; 	MAX : max value plotted, 
;		every data point larger than MAX takes the same color as MAX
;
; 	MIN : min value plotted, 
;		every data point smaller than MIN takes the same color as MIN
;
;	NESTED: specify that the online file is ordered in the nested scheme
;
; 	NOBAR : if set, the color bar (or the color wheel used for POLARIZATION=2)
;                     is hidden
;
;	NOLABELS : if set, no color bar label (min and max) is present
;
;	NOPOSITION : if set, the astronomical location of the map
;	        central point is not indicated
;               ** gnomview only **
;
;       NO_DIPOLE: if set to 1 (and GAL_CUT is not set) 
;                the best fit monopole *and* dipole over all valid pixels are removed
;                * if GAL_CUT is set to b>0, the best monopole and dipole fit is done on all valid
;                pixels with |galactic latitude|>b (in deg) and is removed from all
;                pixels
;             can not be used together with NO_MONOPOLE 
;             (see: GAL_CUT, NO_MONOPOLE)
;               ** mollview and orthview only **
;
;       NO_MONOPOLE: if set to 1 (and GAL_CUT is not set) 
;                the best fit monopole over all valid pixels is removed
;                * if GAL_CUT is set to b>0, the best monopole fit is done on all valid
;                pixels with |galactic latitude|>b (in deg) and is removed from all
;                pixels
;             can not be used together with NO_DIPOLE 
;             (see: GAL_CUT, NO_DIPOLE)
;               ** mollview and orthview only **
;
;       OFFSET: additive offset to apply to data (default = 0)
;               the data plotted is of the form FACTOR*(data + OFFSET)
;               can be used together with LOG
;               see also : FACTOR, LOG
;               Note : does NOT apply to polarization direction or amplitude
;               when POLARIZATION=3. Will apply to polarization amplitude when POLARIZATION=1.
;
; 	ONLINE: if set, you can put directly A HEALPIX VECTOR on File (and
;    		without header): useful when the vector is already
;    		available on line, and avoid to have to write it on disk
;    		just to be read by mollview
;		N.B. : the content of file_in is NOT altered in the
;		process
;               **  can not be used with /SAVE  **    *** OBSOLETE ***
;
;       OUTLINE : single structure, or set of structures, 
;                 each containing the coordinates of one outline to be overplotted.
;           Each structure should contain the following fields : 
;           - 'COORD' coordinate system (either, 'C', 'G', or 'E') of the contour
;           - 'RA'  or longitude coordinates (array)
;           - 'DEC' or lattitude coordinates (array of same size)
;           and can optionally contain the fields
;           - 'LINE[STYLE]' : +2 : black dashes
;                           +1 : black dots
;                            0 : black solid [default]
;                           -1 : black dots on white background
;                           -2 : black dashes on white background
;           - 'PSY[M]' symbol used to represent vertices of outline
;                    (same meaning as standard PSYM in IDL,
;                     if 9<=abs(psym)<=46, D. Fanning's CGSYMCAT symbols 
;                     definition will be used, for example psym=9 is an open circle)
;                    if <=0, the vertices are represented with the chosen symbols, and
;                        connected, by arcs of geodesics.
;                    if >0, only the vertices are shown
;                    (default = 0)
;           - 'SYM[SIZE]' symbol size (same meaning as SYMSIZE in IDL), (default = 1)
;          The line and symbol thickness are controlled (indirectly) via !P.THICK.
;          Outline can be either a single structure, or an array of structures,
;          or a structure of structures.
;           See for instance the function Outline_Earth to create a structure outlining
;          the Earth continents, rivers and/or countries.
;
;	PDF : string containing the name of a PDF output
;             if set to 0 or not set : output to screen
;	      if set to 1            : output the plot in plot_XXX.pdf
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;	      if set to a file name  : output the plot in that file 
;               (see: CROP, GIF, JPEG, PNG, PREVIEW and PS)
;
;      PFONTS:  2-element vector of integers [p0, p1] selecting the default IDL font 
;              of character strings such as theSubtitle, Title and Units.
;          p_0 must be in {-1,0,1} and selects the origin of the fonts among
;              -1:  Hershey Vector, 
;               0:  Device Specific, 
;               1:  True Type Fonts.
;          p_1 must be in {2,...,20} and selects the starting font of the character strings.
;              The font can be changed within each string with embedded formatting commands, 
;              as discussed on http://www.exelisvis.com/docs/Fonts_and_Colors.html.
;          default=[-1,6], corresponding to the Hershey vector font of type 'Complex Roman', 
;            and is equivalent to typing
;             !p.font=-1
;            and prepending the Subtitle, Title and Units strings with '!6'.
;          Note that PFONTS will be ignored if Latex=2 *and* PDF or PS are set.
;
;
;	PNG : string containing the name of a .PNG output
;	      if set to 0 or not set : no .PNG done
;	      if set to 1            : output the plot in plot_XXX.png
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;	      if set to a file name  : output the plot in that file 
;             (see also : CROP, GIF, JPEG, PDF, PNG, PS and PREVIEW)
;
;       POLARIZATION: 
;         if set to 0, no polarization information is plotted.
;
;         otherwise, and if the input data contains polarisation information
;             (ie, Stokes parameter Q and U for each pixel)
;
;         if set to 1 
;             the AMPLITUDE P = sqrt( U^2 + Q^2 ) of the polarisation is plotted
;
;         if set to 2 
;             the ANGLE phi = 0.5*ATAN(U/Q) of the polarisation is plotted
;             Note: the angles are color coded with a fixed color table (independant of Colt)
;
;         if set to 3 or [3, scale_factor, step_factor, thickness]
;             -the temperature is color coded (with a color table defined by Colt)
;             -and the polarisation is overplot as small rods.
;             Polarization can be a 4-element vector (the first element being 3).
;             The second element controls the average length of the rods (default=1), 
;             the third one controls their spacing (default=1),
;             while the fourth one controls their thickness 
;             (which also depends in a device dependent manner on !P.THICK, default=1).
;             Non positive values are replaced by 1.
;
;	PREVIEW : if set, there is a preview of the GIF, JPG, PDF, PostScript,
;	         PNG file being produced  
;                    (see: CROP, GIF, JPEG, PDF, PNG and PS)
;
;	PS :  string containing the name of a PS output
;             if set to 0 or not set : output to screen
;	      if set to 1            : output the plot in plot_XXX.ps
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;	      if set to a file name  : output the plot in that file 
;               (see: CROP, GIF, JPEG, PNG and PREVIEW)
;
; 	PXSIZE: number of horizontal screen_pixels / postscript_color_dots of the plot
;    		** mollview : default = 800, gnomview and cartview : default = 500 **
;    		(useful for high definition color printer)
;
; 	PYSIZE: number of vertical screen_pixels or postscript_color_dots of the plot
;    		default = PXSIZE
;    		(useful for high definition color printer)
;                ** gnomview only **
;
;       RETAIN: backing store for graphics windows in {0,1,2}. Default=2
;
;       RESO_ARCMIN: resolution of the gnomic map in arcmin
;       (default=1.5)
;                ** gnomview only **
;
; 	ROT :   vector with 1, 2 or 3 elements specifing the rotation angles in DEGREE
;               to apply to the map in the 'output' coordinate system (see coord)
;             = ( lon0, [lat0, rat0]) 
;               lon0 : longitude of the point to be put at the center of the plot
;		       the longitude increases Eastward, ie to the left of the plot 
;                      (unless flip is set)
; 		      =0 by default
;               lat0 : latitude of the point to be put at the center of the plot
; 		      =0 by default
;               rot0 : anti clockwise rotation to apply to the sky around the
;                     center (lon0, lat0) before projecting
;                     =0 by default
;
; 	SAVE: if set, tells that the file is in IDL saveset format, 
;    		the variable saved should be DATA 
;                 ** can not be used with /ONLINE **
;
;       SHADED: if set, the orthographic sphere is shaded, using a Phong model, to emulate 3D viewing.
;              The sphere is illuminated by isotropic ambiant light plus a single light source.
;                 ** Can NOT be used with GIF. **
;                   ** orthview only **
;
;       SILENT: if set, the code runs silently
;
;       SILHOUETTE: if set to a scalar or 2-element vector with silhouette[0] !=  0,
;          a silhouette is drawn around the map.
;          Its thickness is proportional to abs(silhouette[0])
;           (and also depends in a device-dependent manner on !P.THICK)
;          Its color is determined by abs(silhouette[1]) in [0,255] (defaulting to 0 = FG_COLOR)
;             ** mollview and orthview only **
;
;       STAGGER: scalar or 2 element vector.
;            - if stagger[0] is in ]0,2], 
;             3 copies of the same sphere centered at [-stagger[0], 0, stagger[0]]
;             (expressed in radius units) along the plot horizontal axis are
;             shown in ORTHOGRAPHIC projection
;             - stagger[1] (if defined), defines the angle of rotation (in degrees) applied
;               to the left and right partial spheres:
;             the lhs sphere is rotated downward by the angle provided, while the rhs one
;             is rotated upward. Rotations are swapped if FLIP is set.
;               ** orthview only **
;
; 	SUBTITLE : String containing the subtitle to the plot (see TITLEPLOT)
;
; 	TITLEPLOT : String containing the title of the plot, 
;     		if not set the title will be File (see SUBTITLE)
;
;       TRANSPARENT: some pixels are transparent in the produced PNG file
;            if set to 1: bad pixels (usually grey or BAD_COLOR) are transparent
;            if set to 2: background pixels (usually white or BG_COLOR) are transparent
;            if set to 3: all of the above
;            only valid with PNG
;
;       TRUECOLORS: if the input data is of the form [Npix,3] then the 3 fields
;            are respectively understood as {Red, Green, Blue} True Colors
;
;
;	UNITS : character string containing the units, to be put on the right
;		side of the color bar (see : NOBAR)
;
;       WINDOW: IDL window index (integer)
;                 if WINDOW < 0: virtual window: no visible window opened. Can be
;               used with PNG or GIF. The Z buffer will be used instead of the 
;               X server, allowing much faster production of the image over a slow network
;                 if WINDOW in [0,31]: the specified IDL window with index WINDOW is used
;               (or reused)
;                 if WINDOW > 31: a free (=unused) window with a random index > 31 will be
;               created and used : default

;	XPOS : The X position on the screen of the lower left corner
;	        of the window, in device coordinate
;
;	YPOS : The Y position on the screen of the lower left corner 
;               of the window, in device coordinate
;
; NOTES
; 	this routine doesn't use the IDL map_set because it is precisely bugged for 
; 	the mollweide projection (at least in IDL 4.0)
;
; SIDE EFFECTS
; 	this routine uses ghostview when PS and PREVIEW are selected 
;	or xv when GIF or PNG and PREVIEW are selected
;
; EXAMPLES
;       ;to plot the signal of the COBE-DMR 4 year map at 53 GHz
;       read_fits_sb, 'dmr_skymap_53a_4yr.fits', dmr53a, /merge  ; read it only one time
;       mollview, dmr53a, /online, 'Sig', /quad
;
;       ;to plot it in Galactic coordinate instead of Ecliptic
;       mollview, drm53a, /online, 'Sig', /quad, coord='g'
;
; COMMONS USED : view_data
;
; PROCEDURES USED: 
;       in the Healpix package :
;	  index_word, read_fits_sb, vec2pix_ring, vec2pix_nest, euler_matrix
;         see  http://www.tac.dk/~healpix
;       it also requires the IDL astro library
;         http://idlastro.gsfc.nasa.gov/homepage.html
;       and the COBE analysis software
;         http://www.gsfc.nasa.gov/astro/cobe/cgis.html
;
; MODIFICATION HISTORY:
; 	October 1997, Eric Hivon, TAC
; 	Nov, 5, 1997,  correction of some bugs for postscript output
; 	13-Nov-97, H. Dole, IAS: save and log keywords
;  	4-Dec-97, H. Dole, IAS: online keyword
; 	16-Dec-97, E.H, TAC: added pxsize, hxsize, subtitle, nobar
;	17-Dec-97, split the loop for projection, added nolabels
;	March-98, E.H. added UNITS keyword
;	April-May-98 E.H. added NESTED_ONLINE, XPOS, YPOS, NOPREVIEW
;       March-99     E.H. Caltech, improved the GIF output
;              modified to deal with structures
;              added Select, COORD, ROT, QUADCUBE  suppressed LON0
;       April-99     E.H. Caltech, improved graticule
;       Nov-99         added flip
;       Feb-00   added rmmonopole and dipole, changed common
;       March-00   changed to no_monopole and no_dipole, changed common
;       Sept-00    added polarisation plotting (Polarization)
;       June-02  : EH, Caltech. Hacked G. Giardino's polview into cartview
;       June-02    partial consolidation of gnomview/mollview/cartview
;       Jan-07    added WINDOW keyword
;       Jun-07:  edited doc header about default data to plot from cut sky file
;       Sep-07:  added /SILENT
;       Mar-08:  added GLSIZE and IGLSIZE
;       Apr-08:  can deal with cut sky data set without creating full sky map
;       Nov-08:  restore original color table and plot settings when exiting
;       May-09:  added /SHADED to orthview, implemented EXECUTE in orthview, fix
;              Min-Max for LOG, use Z buffer when window<0, added RETAIN keyword
;       Oct-09:  added /TRUECOLORS to all routines and MAP_OUT= to Gnomview
;       Apr-10:  accept array of structures in Outline; added MAP_OUT= to
;          Cartview and Mollview
;       Jan-12: added STAGGER to orthview; created azeqview; added JPEG to all
;       Sep-13: added BAD_COLOR, BG_COLOR, FG_COLOR
;       2015:   added PDF, LATEX and PFONTS
;       2017-05: addition of CUSTOMIZE, DEFAULT_SETTINGS, 
;                improvement of GLSIZE and IGLSIZE
;-

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

@viewcom ; define common
data_plot = 0 ; empty common array
; record original color table and PLOTS settings
record_original_settings, original_settings

loadsky                         ; cgis package routine, define rotation matrices
projection = 'MOLLWEIDE'
routine = 'mollview'

uroutine = strupcase(routine)
if keyword_set(help) then begin
    doc_library,'mollview'
    return
endif

if keyword_set(gif) then begin
    message_gif, code=routine, error=error_gif
    if (error_gif) then return
endif

if (n_params() lt 1 or n_params() gt 2) then begin
    PRINT, 'Wrong number of arguments in '+uroutine
    print,'Syntax : '
    print, uroutine+', File, [Select, ]'
    print,'              [ASINH=, '
    print,'              BAD_COLOR=, BG_COLOR=, '
    print,'              CHARSIZE=, CHARTHICK=, COLT=, COORD=, CROP=, '
    print,'              EXECUTE=, '
    print,'              FACTOR=, FG_COLOR=, FITS=, FLIP=, '
    print,'              GAL_CUT=, GIF=, GLSIZE=, GRATICULE=, '
    print,'              HBOUND=, HELP=, HIST_EQUAL=, HXSIZE=,'
    print,'              IGLSIZE=, IGRATICULE=,'
    print,'              JPEG=,'
    print,'              LOG=, '
    print,'              MAP_OUT=, MAX=, MIN=, '
    print,'              NESTED=, NOBAR=, NOLABELS=, '
    print,'              NO_DIPOLE=, NO_MONOPLE=, '
    print,'              OFFSET=, ONLINE=, OUTLINE=,'
    print,'              PDF=, PFONTS=, PNG=, POLARIZATION=, PREVIEW=, '
    print,'              PS=, PXSIZE=, PYSIZE=, QUADCUBE= ,'
    print,'              RETAIN=, ROT=,  '
    print,'              SAVE=, SILENT=, SUBTITLE=, '
    print,'              TITLEPLOT=, TRANSPARENT=, TRUECOLORS=, '
    print,'              UNITS=, WINDOW=, XPOS=, YPOS=]'
    print
    print,' Type '+uroutine+', /help '
    print,'   for an extended help'
    return
endif

IF (undefined(file_in)) then begin
    print,routine+': Undefined variable as 1st argument'
    return
endif
do_flip = keyword_set(flip)

if (!D.n_colors lt 4) then begin
    print,' : Sorry ... not enough colors ('+strtrim(string(!d.n_colors),2)+') available'
    return
endif

if (keyword_set(no_monopole) and keyword_set(no_dipole)) then begin
    print,routine+': choose either NO_MONOPOLE or NO_DIPOLE'
    print,'    (removal of best fit monopole only or best fit monopole+dipole)'
    return
endif

polar_type = 0
if keyword_set(polarization) then polar_type = polarization

loaddata_healpix, $
  file_in, select_in,$
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
  SAVE=save, ONLINE=online, NESTED=nested_online, UNITS=units, COORD=coord, FLIP=flip, $
  ROT=rot, QUADCUBE=quadcube, LOG=log, ERROR=error, $
  POLARIZATION=polarization, FACTOR=factor, OFFSET=offset, SILENT=silent, COMPRESS=1, PIXEL_LIST=pixel_list, $
  TRUECOLORS=truecolors, DATA_TC=data_tc
if error NE 0 then return

data2moll, $
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
  planmap, Tmax, Tmin, color_bar, planvec, vector_scale, $
  PXSIZE=pxsize, LOG=log, HIST_EQUAL=hist_equal, MAX=max_set, MIN=min_set, FLIP=flip,  $
  NO_DIPOLE=no_dipole, NO_MONOPOLE=no_monopole, UNITS=sunits, DATA_plot = data_plot, GAL_CUT=gal_cut, $
  POLARIZATION=polarization, SILENT=silent, PIXEL_LIST=pixel_list, ASINH=asinh, $
  TRUECOLORS=truecolors, DATA_TC=data_tc, MAP_OUT = arg_present(map_out) ? map_out : -1, ROT=rot, FITS=fits

proj2out, $
  planmap, Tmax, Tmin, color_bar, 0., title_display, $
  sunits, coord_out, do_rot, eul_mat, planvec, vector_scale, $
  CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, $
  HXSIZE=hxsize, NOBAR = nobar, NOLABELS = nolabels, PNG = png, PREVIEW = preview, PS=ps, PXSIZE=pxsize, $
  SUBTITLE = subtitle, TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
  POLARIZATION=polarization, OUTLINE=outline, /MOLL, FLIP=flip, COORD_IN=coord_in, IGRATICULE=igraticule, $
  HBOUND = hbound, WINDOW = window, EXECUTE=execute, SILENT=silent, GLSIZE=glsize, $
  IGLSIZE=iglsize, RETAIN=retain, TRUECOLORS=truecolors, TRANSPARENT=transparent, CHARTHICK=charthick, $
  JPEG=jpeg, BAD_COLOR=bad_color, BG_COLOR=bg_color, FG_COLOR=fg_color, PDF=pdf, $
  LATEX=latex, PFONTS=pfonts, CUSTOMIZE=CUSTOMIZE, DEFAULT_SETTINGS=default_settings, SILHOUETTE = silhouette


w_num = !d.window
; restore original color table and PLOTS settings
record_original_settings, original_settings, /restore


return
end

