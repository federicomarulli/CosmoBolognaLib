#!/bin/sh
# © M E C Swanson 2008
#plot the angular mask described by a polygon file in graphics format
#using sm
#
#optional 3rd-6th arguments give the RA and Dec range to plot the mask
#If no range is given, range will be calculated automatically from the data.
#If only an RA range is given, Dec range will be calculated automatically
#to give equal scaling on x and y axes, centered on the range of all the data.
#
#Entering in "0 0 0 0" for the range will also give an automatic range,
#which can be used if you also wish to, e.g., enter a title. 
#Likewise entering "0 0" for the Dec range is equivalent to only specifying an RA range. 
#
#optional 7th argument gives a title for the graph
#optional 8th argument turns on drawing black outlines around each polygon
#
#if specifying the full path to a file (rather than just running in the directory 
#containing the file you want to plot), put the path in double-quotes as shown below.
#
#USAGE: graphmasksm.sh <infile> <outfile> [<ramin> <ramax>] [<decmin> <decmax>] [<title>] [<outlines>]
#EXAMPLES: 
#default range, no outlines: graphmasksm.sh "dr4/safe0/sdss_dr4safe0_mask.grph" "dr4/safe0/sdss_dr4safe0_mask.eps"
#defined range, title, no outlines: graphmasksm.sh sdss_slice.grph sdss_slice.eps -45 35 8 21 "SDSS slice"
#defined range, no title, outlines: graphmasksm.sh sdss_slice.grph sdss_slice.eps -45 35 8 21 " " on
#default range, no title, outlines: graphmasksm.sh sdss_slice.grph sdss_slice.eps 0 0 0 0 " " on
		
#WARNINGS: (a.k.a., reasons to use the matlab plotting script instead!)		
# -This script uses the opposite color scheme than the matlab script, i.e., 		
#  weight 1 = black, weight 0 = white	       
# -RA and Dec are treated as linear variables (no spherical projection).  This means that 		
#  ranges that span RA=0 won't plot properly.		
# -SM can only plot square eps files, so if you define both RA and Dec range, they should		
#  have the same length, otherwise your plot will appear stretched.		
# -SM has a line length limit of 1500 characters, which SDSS tends to have no trouble overloading.  		
#  If you get an error reading the file, try making your graphics file with lower precision and 
#  fewer points per 2pi, i.e., use poly2poly -og12 -p3 instead of poly2poly -og30.
#


if [ "$MANGLESCRIPTSDIR" = "" ] ; then
    MANGLESCRIPTSDIR="../../scripts"
fi

if [ ! "$7" = "" ] ; then
    title="\"$7\""
else
    title=""
fi

if [ "$8" = "on" ] ; then
    outlines=1
else
    outlines=""
fi

if [ $# -ge 2 ]; then
    if [ -e $1 ]; then 
	echo "Plotting $1 using sm ..."
	if [ -e $2 ]; then
	    rm $2
	fi
	sm -m $MANGLESCRIPTSDIR/graphmask.sm $1 $2 $3 $4 $5 $6 $title $outlines > sm.temp
    else
	echo >&2 "ERROR: file $1 not found."
	exit 1
    fi
else
    echo >&2 "USAGE: graphmasksm.sh <infile> <outfile> [<ramin> <ramax>] [<decmin> <decmax>] [<title>] [<outlines>]"
    echo >&2 "EXAMPLES:" 
    echo >&2 "default range, no outlines: graphmasksm.sh \"dr4/safe0/sdss_dr4safe0_mask.grph\" \"dr4/safe0/sdss_dr4safe0_mask.eps\""
    echo >&2 "defined range, title, no outlines: graphmasksm.sh sdss_slice.grph sdss_slice.eps -45 35 8 21 \"SDSS slice\""
    echo >&2 "defined range, no title, outlines: graphmasksm.sh sdss_slice.grph sdss_slice.eps -45 35 8 21 \" \" on"
    echo >&2 "default range, no title, outlines: graphmasksm.sh sdss_slice.grph sdss_slice.eps 0 0 0 0 \" \" on"
    exit 1
fi
if [ -e smerr.temp ] || [ ! -e $2 ]; then
    echo >&2 "ERROR: error in sm plotting script."
    echo >&2 "See log in sm.temp for more details."
    if [ -e smerr.temp ]; then
	rm smerr.temp
    fi
    exit 1
fi
rm sm.temp
echo all done!
