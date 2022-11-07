#!/bin/sh
# © M E C Swanson 2008
#plot the angular mask described by a polygon file in .list format
#using the matlab mapping toolbox
#
#optional 3rd-6th arguments give the RA and Dec range to plot the mask
#if no range is given, the full sky will be plotted
#To plot range that includes RA=0, use a negative number for the "low"
#end, e.g. for a 10 degree strip covering RA=355 degrees to RA=5 degrees,
#use a range of -5 5.
#Entering in "0 0 0 0" for the range will also give the default full sky behavior,
#which can be used if you also wish to, e.g., enter a title. 
#
#optional 7th argument gives a title for the graph
#optional 8th argument turns on drawing black outlines around each polygon
#
#if specifying the full path to a file (rather than just running in the directory 
#containing the file you want to plot), put the path in double-quotes as shown below.
#
#USAGE: graphmask.sh <infile> <outfile> [<ramin> <ramax> <decmin> <decmax>] [<title>] [<outlines>]
#EXAMPLES: 
#fullsky, no outlines: graphmask.sh "dr4/safe0/sdss_dr4safe0_mask.list" "dr4/safe0/sdss_dr4safe0_mask.eps"
#defined range, title, no outlines: graphmask.sh sdss_slice.list sdss_slice.eps -45 35 8 21 "SDSS slice"
#defined range, no title, outlines: graphmask.sh sdss_slice.list sdss_slice.eps -45 35 8 21 "" on
#fullsky, no title, outlines: graphmask.sh sdss_slice.list sdss_slice.eps 0 0 0 0 "" on

if [ "$MANGLESCRIPTSDIR" = "" ] ; then
    MANGLESCRIPTSDIR="../../scripts"
fi

if [ $# -eq 2 ]; then
    if [ -e $1 ]; then 
	if [ -e $2 ]; then
	    rm $2
	fi
	matlab -nodisplay -r addpath\(\'$MANGLESCRIPTSDIR\'\)\;graphmask\(\'$1\',\'$2\'\)
    else
	echo >&2 "ERROR: file $1 not found."
	exit 1
    fi
elif [ $# -ge 6 ]; then
    if [ -e $1 ]; then 
	if [ -e $2 ]; then
	    rm $2
	fi
	echo "addpath('$MANGLESCRIPTSDIR'); graphmask('$1','$2',[$3,$4,$5,$6],'$7','$8')" > jgraphtemp.m
	matlab -nodisplay < jgraphtemp.m
	rm jgraphtemp.m
    else
	echo >&2 "ERROR: file $1 not found."
	exit 1
    fi
else
    echo >&2 "USAGE: graphmask.sh <infile> <outfile> [<ramin> <ramax> <decmin> <decmax>] [<title>] [<outlines>]"
    echo >&2 "EXAMPLES:" 
    echo >&2 "fullsky, no outlines: graphmask.sh \"dr4/safe0/sdss_dr4safe0_mask.list\" \"dr4/safe0/sdss_dr4safe0_mask.eps\""
    echo >&2 "defined range, title, no outlines: graphmask.sh sdss_slice.list sdss_slice.eps -45 35 8 21 \"SDSS slice\""
    echo >&2 "defined range, no title, outlines: graphmask.sh sdss_slice.list sdss_slice.eps -45 35 8 21 \"\" on"
    echo >&2 "fullsky, no title, outlines: graphmask.sh sdss_slice.list sdss_slice.eps 0 0 0 0 \"\" on"
    exit 1
fi
if [ -e matlabexit.temp ] || [ ! -e $2 ]; then
    echo >&2 "ERROR: error in matlab plotting script."
    if [ -e matlabexit.temp ]; then
	rm matlabexit.temp
    fi
    exit 1
fi
echo all done!
