/**
Package implementing essential HEALPix functionality.<br>
healpix.essentials is a Java package aiming to combine the advantages
of other Java HEALPix libraries like
the gov.fnal.eag.healpix package by N. Kuropatkin
(<a href="http://home.fnal.gov/~kuropat/HEALPIX/PixTools.html">http://home.fnal.gov/~kuropat/HEALPIX/PixTools.html</a>),
Jan Kotek's org.asterope.healpix package
(<a href="https://github.com/jankotek/HEALPix">https://github.com/jankotek/HEALPix</a>),
and the healpix.core package developed at ESO.
<br>
Current features include:
<ul>
<li>close similarities with Healpix_Base_T class from Healpix C++, which allows
  simultaneous development and bug fixes for both.
<li>support for arbitrary positive Nside values in RING scheme; no longer
  limited to powers of 2
<li>maximum supported Nside value: 2^29
<li>significant performance improvements: most methods have been accelerated
  by integral factors, some by more than an order of magnitude.
<li>re-implementation of queryDisc and queryPolygon:
<ul>
  <li> query methods return RangeSet objects which allow much more compact
    storage of the result
  <li> new native query methods for NESTED ordering; these are slower than those
    for RING ordering, but much quicker than converting all pixels from a RING
    result to NESTED.
  <li> inclusive queries have been improved: several bugs were fixed, and the
    number of false positives in the result has been reduced. Users can now
    choose between quick inclusive queries returning more false positives,
    and slower ones returning fewer false positives.
</ul>
<li> the HealpixProc class offers a procedural (instead of object-oriented)
  interface to the HealpixBase functionality, which simplifies transition
  for users of the "Healpix" and "PixTools" classes.
  NOTE: this only works for Nside parameters which are powers of 2
<li> many bug fixes
<li> no external library dependencies, except for "nom.tam.fits" if FITS I/O is
  required
<li> the code base is thread-safe in the following sense:
<ul>
  <li> HealpixProc methods can be called concurrently
  <li> HealpixBase methods can be called concurrently on different objects
</ul>
</ul>
@copyright 2011, 2012 Max-Planck-Society
@author Martin Reinecke
*/
package healpix.essentials;
