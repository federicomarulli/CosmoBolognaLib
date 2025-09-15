
RULES :
------
 The permission to use and copy this software and its documentation, 
 without fee or royalty is limited to non-commercial purposes related to 
 PLANCK Surveyor project and provided that you agree to comply with
 the following copyright notice and statements,
 and that the same appear on ALL copies of the software and documentation.

 An appropriate acknowledgement has to be included in any
 publications based on work where the package has been used.

 Copyright 1997 by Eric Hivon and Krzysztof M. Gorski.
  All rights reserved.

CONTENTS :
--------
 this directory contains IDL routines to project data in HEALPIX tesselation,
 create maps, plot them, and I/O them in FITS format.
 (symbols :   (*) new,   (+) recommended,   (-) obsolete)

 README.txt   : this file

 examples/
*   coadd_raw_obs.pro	: procedure to transform a set of RA,DEC,MEASURE into a HEALPIX
			  signal  and  N_obs   maps (in galactic coordinate) 
			  for illustration purpose only
*   raw_obs.txt 	: ASCII file used for input of coadd_row_obs.pro

 fits/	      : fits I/O routines
*   extract_ttype.pro 	: function to extract the value of the keywords TTYPEi 
			  in an extension header
    fxbread.pro		: routine from the IDL ASTRONOMY USER'S LIBRARY 
                          corrected from a bug (in versions older than June 98)
			  if needed make sure that this one is used instead of the original one
			  (by setting correctly your !path variable in IDL)
			  (used by mollview)
-   read_fits_map.pro  	: procedure to read FITS file images and bintable
*+  read_fits_s.pro  	: procedure to read FITS file images, ascii and bintable
			  using structures as interface, 
			  much more flexible than write_fits_map
			  (used by mollview)
*   sxaddpar.pro	: routine from the IDL ASTRONOMY USER'S LIBRARY 
			  to write FITS keywords modified to deal with several COMMENT's
-   write_fits_map.pro 	: procedure to write bintable in FITS file
*+  write_fits_sb.pro   : procedure to write images and bintable in FITS file 
			  using structures as interface, 
			  much more flexible than write_fits_map

 misc/		: misc. IDL routines

 toolkit/      : set of routines for conversion to and from HEALPIX
*+  ang2pix_nest.pro	: procedure to convert (theta,phi) position on the sphere to 
                       	  pixel number for HEALPIX pixelisation and NESTED scheme
			  (previously ang_pix_nest.pro, used by mollview)
-   ang_pix_nest.pro    : obsolete, replaced by ang2pix_nest.pro

*+  ang2pix_ring.pro 	: procedure to convert (theta,phi) position on the sphere to 
                      	  pixel number for HEALPIX pixelisation and RING scheme
			  (previously ang_pix.pro, used by mollview)
-   ang_pix.pro    	: obsolete, replaced by ang2pix_ring.pro

*   init_pix2xy      	: initialize some arrays necessary to conversion routines

*   init_xy2pix	       	: initialize some arrays necessary to conversion routines
			  (used by mollview)
*   nest2ring.pro	: procedure to convert NESTED pixel number to RING pixel number

*   pix2ang_nest.pro 	: procedure to convert pixel number to 
			  (theta,phi) position on the sphere 
                          for HEALPIX pixelisation and NESTED scheme
*+  pix2ang_ring.pro  	: procedure to convert pixel number to 
			  (theta,phi) position on the sphere 
                          for HEALPIX pixelisation and RING scheme
			  (previously pix_ang.pro)
-   pix_ang.pro		: obsolete, replaced by pix2ang_ring.pro

*   pix2vec_nest.pro  	: procedure to convert pixel number to vector position on the sphere 
                          for HEALPIX pixelisation and NESTED scheme
*   pix2vec_ring.pro  	: procedure to convert pixel number to vector position on the sphere 
                          for HEALPIX pixelisation and RING scheme
*   ring2nest.pro	: procedure to convert RING pixel number to NESTED pixel number

*   vec2pix_nest.pro  	: procedure to convert vector position on the sphere to 
                          pixel number for HEALPIX pixelisation and NESTED scheme
*   vec2pix_ring.pro   	: procedure to convert vector position on the sphere to 
                          pixel number for HEALPIX pixelisation and RING scheme

 visu/        : visualisation  routines
    mollview.pro     	: Mollweide projection of the full sky
    mollcursor.pro      : show the latitude, longitude and pixel value on a mollweide map
    gnomview.pro        : Gnomic projection of a part of sky


TO GET MORE INFORMATION :
------------------------
 See ../html/idl_index.html with your Web brower
 or type doc_library,'routine name' under IDL


FEED BACK :
---------
 If you have any questions, comments or suggestions regarding 
 this package, please email them to  hivon@tac.dk

NOTE :
-----
 It requires the THE IDL ASTRONOMY USER'S LIBRARY 
 that can be found at
    http://idlastro.gsfc.nasa.gov/homepage.html
  latest version (jan 99) written for IDL V5.0 or newer
  a version for IDL 4 can also be found there

 It requires the COBE Analysis Software
   http://www.gsfc.nasa.gov/astro/cobe/cgis.html
  latest version = 4.0,      July 1998

Technical note (about colors in IDL)
---------------
It has been noted that, on various terminals
older versions of mollview could feature an unexpected color table
when executed immediately after opening the IDL session
(but following runs within the same session are just fine).

I don't know of patches working on all terminals for this
IDL strange behaviour with respect to colors.
Ask you local expert ..

update (march 1999)
The problem seems to be resolved on current version of mollview,
if not, let me know (hivon@tac.dk or efh@astro.caltech.edu)

Eric Hivon


