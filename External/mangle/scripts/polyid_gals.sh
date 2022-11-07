#!/bin/sh
# © M E C Swanson 2008
#
#This script uses polyid to find in which polygons in the input mask
#each galaxy in the list lies.  The input galaxy file can have any number of columns,
#but the first 2 columns should be the RA and dec.
#output galaxy files for each resolution will attach the polygon number to the end of 
#each line in the galaxy file.
#
#By default, galaxies that are not in any polygon in the mask will be cut out of
#the output file.  To keep all galaxies in the output file, enter a 0 for the final
#command line argument.
#
#USAGE: polyid_gals.sh <input mask file> <input galaxy file> <output galaxy file> [<0 to turn cuts off>]
#EXAMPLES:
#polyid_gals.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat sdss_dr7safe0_polyidgals_cut.dat
#polyid_gals.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat sdss_dr7safe0_polyidgals_uncut.dat 0


if [ "$MANGLEBINDIR" = "" ] ; then
    MANGLEBINDIR="../../bin"
fi
if [ "$MANGLESCRIPTSDIR" = "" ] ; then
    MANGLESCRIPTSDIR="../../scripts"
fi
if [ "$MANGLEDATADIR" = "" ] ; then
    MANGLEDATADIR="../../masks"
fi

mask=$1
gals=$2
outgals=$3
cut=$4

#check command line arguments
if [ "$mask" = "" ] || [ "$gals" = "" ] || [ "$outgals" = "" ]; then
    echo >&2 "ERROR: enter the input polygon file and input and output galaxy files as command line arguments."
    echo >&2 "" 
    echo >&2 "USAGE: polyid_gals.sh <input mask file> <input galaxy file> <output galaxy file> [<0 to turn cuts off>]"
    echo >&2 "EXAMPLES:"
    echo >&2 "polyid_gals.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat sdss_dr7safe0_polyidgals_cut.dat"
    echo >&2 "polyid_gals.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat sdss_dr7safe0_polyidgals_uncut.dat 0"
    exit 1
fi

if [ "$cut" = "" ] ; then
    cut=1
fi
	  
#count number of columns in input galaxy file
head -n 1 $gals > jgals
numfields=`awk '{print NF}' jgals`
rm jgals
	  
#run polyid on galaxy file and assemble output galaxy files with polygon numbers
echo "Running polyid to find galaxies from $gals in $mask  ... "
$MANGLEBINDIR/polyid $mask $gals j1 || exit
tail +2 j1 > j2
awk '{print $3}' j2 > j3
paste $gals j3 > j4
if [ $cut -ne 0 ] ; then
    awk "NF == ($numfields+1) {print \$0}" j4 > $outgals
else
    cp j4 $outgals
fi
rm j1 j2 j3 j4

count=`wc -l < $outgals`
echo ""
echo "wrote $count galaxies to $outgals."

