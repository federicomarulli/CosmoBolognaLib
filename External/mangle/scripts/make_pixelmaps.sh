#!/bin/sh
# © M E C Swanson 2008
#Script to make pixelmaps of a mask at several resolutions 
#Pixelmaps give the average value of the mask's weight in each pixel,
#using one of mangle's pixelization schemes.
#pixelmaps are saved as $inputmask.$res$scheme, e.g. for 
#the input mask sdss_dr7safe0_res6d.pol, pixelmap files will be
#sdss_dr7safe0_res6d.pol.1d, sdss_dr7safe0_res6d.pol.2d, etc.
#
#if a galaxy file is given, this script also uses polyid to find in which pixels
#each galaxy in the list lies.  The input galaxy file can have any number of columns,
#but the first 2 columns should be the RA and dec.
#output galaxy files for each resolution will attach the pixel number to the end of 
#each line in the galaxy file that lies within a pixel in the pixelmap.  The naming
#scheme is the same as for the pixelmap files, e.g., for the input galaxy file
#sdss_dr7safe0_galdetails1.dat, the output galaxy files will be
#sdss_dr7safe0_galdetails1.dat.1d, sdss_dr7safe0_galdetails1.dat.2d, etc.
#
#By default, this script cuts out any pixels with average weights less than 0.5 for
#resolution 1 and average weights less than 0.8 for higher resolutions.  To run
#the script without these cuts, enter a 0 for the final command line argument.
#
#USAGE: make_pixelmaps.sh <input mask file> [<input galaxy file>] [<0 to turn cuts off>]
#EXAMPLES:
#make_pixelmaps.sh sdss_dr7safe0_res6d.pol 
#make_pixelmaps.sh sdss_dr7safe0_res6d.pol 0
#make_pixelmaps.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat
#make_pixelmaps.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat 0


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
cutweights=$3

#check command line arguments
if [ "$mask" = "" ] ; then
    echo >&2 "ERROR: enter the input polygon and galaxy files as command line arguments."
    echo >&2 "" 
    echo >&2 "USAGE: make_pixelmaps.sh <input mask file> [<input galaxy file>] [<0 to turn cuts off>]"
    echo >&2 "EXAMPLES:"
    echo >&2 "make_pixelmaps.sh sdss_dr7safe0_res6d.pol" 
    echo >&2 "make_pixelmaps.sh sdss_dr7safe0_res6d.pol 0"
    echo >&2 "make_pixelmaps.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat"
    echo >&2 "make_pixelmaps.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat 0"
    exit 1
fi

if [ "$gals" == 0 ] ; then
    cutweights=0
    gals=
fi
if [ "$cutweights" = "" ] ; then
    cutweights=1
fi

dresmax=6
dscheme="d"

head -n 100 $mask > jmaskhead

#grab pixelization info from input file
awk '/pixelization/{print $0}' < jmaskhead > jpix
resmax=`awk '{print substr($2, 1, length($2)-1)}' < jpix`
scheme=`awk '{print substr($2, length($2))}' < jpix`
rm jpix

#check if input file is snapped and balkanized
snapped=`awk '/snapped/{print $1}' < jmaskhead`
balkanized=`awk '/balkanized/{print $1}' < jmaskhead`
rm jmaskhead

#if input file is not pixelized, snapped, and balkanized, exit.
#if input file is pixelized to a fixed resolution, use it as is.
if [ "$resmax" = "" ]; then
    echo >&2 "ERROR: cannot make pixelmaps of an unpixelized mask."
    echo >&2 "Pixelize your mask using the pixelize function first."
    exit 1
elif [ "$resmax" = -1 ] ; then
    resmax=$dresmax
    echo "WARNING: making pixelmaps of a mask pixelized adaptively."
    echo "Attempting to make pixelmaps up to resolution $resmax."
    echo "To guarantee that pixelmaps can be made up to a given resolution,"
    echo "Pixelize your mask using a fixed resolution <r> by using 'pixelize -P<scheme>0,<r>'."
fi

#if input file isn't snapped, exit
if [ ! "$snapped" = "snapped" ]; then
    echo >&2 "ERROR: cannot make pixelmaps of an unsnapped mask."
    echo >&2 "Run snap on your mask before running this script,"
    echo >&2 "or add the 'snapped' keyword to your file if it is already snapped."
    exit 1
fi

#if input file isn't balkanized, exit
if [ ! "$balkanized" = "balkanized" ]; then
    echo >&2 "ERROR: cannot make pixelmaps of an unbalkanized mask."
    echo >&2 "Run balkanize on your mask before running this script,"
    echo >&2 "or add the 'balkanized' keyword to your file if it already consists"
    echo >&2 "of only non-overlapping polygons."
    exit 1
fi

for (( res=1; res<=$resmax; res++ ))
  do
  
  outmask=$mask.$res$scheme
  echo "Generating pixelmap for resolution $res ..."
  $MANGLEBINDIR/pixelmap -P${scheme}0,$res $mask $outmask || exit
  
#cut all pixels with weights less than 0.5 for resolution 1 or 0.8 for higher resolutions
  if [ $cutweights -ne 0 ] ; then
      if [ $res -le 1 ] ; then 
	  weightcut=0.5
      else
	  weightcut=0.8
      fi
      mv $outmask jp
      $MANGLEBINDIR/poly2poly -j$weightcut jp $outmask || exit
      echo "Pixels in $outmask have weights in the range defined in make_pixelmaps.sh."
      echo "To run make_pixelmaps.sh without these cuts, enter a 0 for the last command line argument."
      echo "EXAMPLE: make_pixelmaps.sh sdss_dr7safe0_res6d.pol sdss_dr7safe0_galdetails1.dat 0"
  fi
  

  if [ ! "$gals" == "" ] ; then
      if [ -e "$gals" ] ; then 
	  outgals=$gals.$res$scheme
	  $MANGLESCRIPTSDIR/polyid_gals.sh $outmask $gals $outgals
      else 
	  echo >&2 "ERROR: input galaxy file $gals not found!"
	  exit 1
      fi
  fi

done

