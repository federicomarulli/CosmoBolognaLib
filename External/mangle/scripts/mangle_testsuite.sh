#!/bin/sh
# © M E C Swanson 2008
#
#This script tests the mangle installation and creates a tarball of
#output files for further examination.
#
#USAGE: source run_mangle_testsuite > test.log
mkdir test

echo "mangle environment variables are:"
echo "MANGLEBINDIR=$MANGLEBINDIR"
echo "MANGLESCRIPTSDIR=$MANGLESCRIPTSDIR"
echo "MANGLEDATADIR=$MANGLEDATADIR"
echo "PATH=$PATH"

echo "Checking if unformatted fortran files are readable."
cd ../masks/2df100k
../../bin/weight -z2dF100k ngp_fields.dat jnf || exit
../../bin/weight -z2dF100k sgp_fields.dat jnf || exit
cd ../2df230k
../../bin/weight -z2dF230k ngp_fields.dat jnf || exit
../../bin/weight -z2dF230k sgp_fields.dat jnf || exit
rm jnf

echo "Running mangle on 2qz10k mask ..."
cd ../2qz10k/
./2qz.sh
if [ -e 2qz_res4s.eps ]; then
    mv 2qz_res4s.eps ../../scripts/test
    if [ -e 2qz_north_res4s.eps ]; then
	mv 2qz_north_res4s.eps ../../scripts/test
    fi
    if [ -e 2qz_south_res4s.eps ]; then
	mv 2qz_south_res4s.eps ../../scripts/test
    fi
fi
mv 2qz_north_res4s.pol 2qz_south_res4s.pol  ../../scripts/test

rm 2qz_*

samp=dr6
echo "Running mangle on slice of SDSS $samp mask ..."
cd ../sdss
./sdss_quickstart.sh
if [ -e sdss_${samp}safe0_slice.eps ]; then 
    mv sdss_${samp}safe0_slice.eps ../../scripts/test
elif [ -e sdss_${samp}safe0_slice1.eps ] && [ -e sdss_${samp}safe0_slice2.eps ]; then
    mv sdss_${samp}safe0_slice1.eps sdss_${samp}safe0_slice2.eps ../../scripts/test
else
    mv sdss_${samp}safe0_slice.pol ../../scripts/test
fi
rm sdss_${samp}safe0_slice.*

echo "Trimming 2qz north mask with an icosahedron polygon ..."
cd ../../scripts/test
cp ../../masks/icosahedron/icosahedron.pol .
../../bin/poly2poly -J7,7 icosahedron.pol ico7.pol
../../scripts/trim_mask.sh 2qz_north_res4s.pol ico7.pol trimmed_mask.pol
rm icosahedron.pol ico7.pol

echo "Rasterizing 2qz north mask ... "
rm ../../masks/healpix/healpix_polys/nside16_p5s.pol
../../scripts/healpixrast.sh 2qz_north_res4s.pol 16 rasterized_mask.pol
../../scripts/healpixrast2fits.sh 2qz_north_res4s.pol 16 rasterized_mask.fits 16 rasterized_mask.gif
../../scripts/call ../../bin/fits2dat_binary.x 1 16 rasterized_mask.fits j2
echo healpix_weight 3072 > j1
cat j1 j2 > rasterized_mask1.dat
rm j1 j2 args.dat
../../bin/poly2poly rasterized_mask1.dat rasterized_mask1.pol

echo "Making pixelmaps of 2qz north mask ..."
cp ../../masks/2qz10k/azel.dat jazel
tail +2 jazel > azel.dat
../../scripts/make_pixelmaps.sh 2qz_north_res4s.pol azel.dat 0
rm azel.dat jazel

if which matlab >/dev/null 2>&1 ; then
    ../../bin/poly2poly -ol30 trimmed_mask.pol trimmed_mask.list
    ../../scripts/graphmask.sh trimmed_mask.list trimmed_mask.eps
    rm trimmed_mask.list*
    ../../bin/poly2poly -ol30 rasterized_mask.pol rasterized_mask.list
    ../../scripts/graphmask.sh rasterized_mask.list rasterized_mask.eps
    rm rasterized_mask.list*
    ../../bin/poly2poly -ol30 rasterized_mask1.pol rasterized_mask1.list
    ../../scripts/graphmask.sh rasterized_mask1.list rasterized_mask1.eps
    rm rasterized_mask1.list*
    for(( i=1; i<=4; i++ ))
      do
      ../../bin/poly2poly -ol30 2qz_north_res4s.pol.${i}s jpix.list
      graphmask.sh jpix.list pixmap${i}.eps
      rm jpix.list*
    done
elif which sm >/dev/null 2>&1 ; then
    ../../bin/poly2poly -og30 trimmed_mask.pol trimmed_mask.grph
    ../../scripts/graphmasksm.sh trimmed_mask.grph trimmed_mask.eps
    rm trimmed_mask.grph
    ../../bin/poly2poly -og30 rasterized_mask.pol rasterized_mask.grph
    ../../scripts/graphmasksm.sh rasterized_mask.grph rasterized_mask.eps
    rm rasterized_mask.grph
    ../../bin/poly2poly -og30 rasterized_mask1.pol rasterized_mask1.grph
    ../../scripts/graphmasksm.sh rasterized_mask1.grph rasterized_mask1.eps
    rm rasterized_mask1.grph
    for(( i=1; i<=4; i++ ))
      do
      ../../bin/poly2poly -og30 2qz_north_res4s.pol.${i}s jpix.grph
      ../../scripts/graphmasksm.sh jpix.grph pixmap${i}.eps
      rm jpix.grph*
    done
fi
cd ..
echo >&2 "mangle test suite complete!  Output files are in test.tar.gz"
mv test.log test
tar cfz test.tar.gz test

