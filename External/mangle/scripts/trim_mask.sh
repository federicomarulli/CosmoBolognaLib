#! /bin/sh
# © M E C Swanson 2008
#script to trim a mask so that it only includes polygons in a specified region
#USAGE: trim_mask.sh <mask> <trimmer> <outfile> [<pixelization arguments>]
#EXAMPLE:trim_mask.sh mask.pol trimmer.pol trimmed_mask.pol
# You can also optionally provide arguments to the pixelize function
# in the 4th argument:
#EXAMPLE: trim_mask.sh mask.pol trimmer.pol trimmed_mask.pol -Ps0,6


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
trimmer=$2
outfile=$3
pixargs=$4

#check command line arguments:
if [ "$mask" = "" ] || [ "$trimmer" = "" ] || [ "$outfile" = "" ] ; then
    echo >&2 "ERROR: enter the input polygon file, a polygon defining the region you" 
    echo >&2 "want to trim to, and the output polygon file as command line arguments."
    echo >&2 "" 
    echo >&2 "USAGE: trim_mask.sh <mask> <trimmer> <outfile> [<pixelization arguments>]"
    echo >&2 "EXAMPLE: trim_mask.sh mask.pol trimmer.pol trimmed_mask.pol"
    echo >&2 "EXAMPLE w/ optional argument to pixelize: trim_mask.sh mask.pol trimmer.pol trimmed_mask.pol -Ps0,6"   
    exit 1
fi

#if no argument for pixelize is given, pixelize to resolution 6 with the simple scheme:
if [ "$pixargs" = "" ] ; then
    pixargs="-Ps0,6"
fi

head -n 100 $mask > jmaskhead
head -n 100 $trimmer > jtrimmerhead

#grab pixelization info from input files
awk '/pixelization/{print $0}' < jmaskhead > jpix
res1=`awk '{print substr($2, 1, length($2)-1)}' < jpix`
scheme1=`awk '{print substr($2, length($2))}' < jpix`
awk '/pixelization/{print $0}' < jtrimmerhead > jpix
res2=`awk '{print substr($2, 1, length($2)-1)}' < jpix`
scheme2=`awk '{print substr($2, length($2))}' < jpix`
rm jpix jmaskhead jtrimmerhead

#if input files are pixelized, make sure they are consistent:
if [ ! "$res1" = "" ] && [ ! "$res2" = "" ] ; then
    if [ $res1 -eq -1 ] || [ $res2 -eq -1 ] ; then
	echo >&2 "ERROR: cannot trim a mask pixelized adaptively."
	echo >&2 "Pixelize both your mask and your trimmer polygon(s) using a fixed resolution," 
	echo >&2 "e.g. -Ps0,8, and try again."
	exit 1
    fi
    if [ $res1 -ne $res2 ] ; then
	echo >&2 "ERROR: mask polygons and trimmer polygons must be pixelized to the same resolution."
	echo >&2 "Pixelize both your mask and your trimmer polygon(s) using a fixed resolution," 
	echo >&2 "e.g. -Ps0,8, and try again."
	exit 1
    fi    
#if input files are unpixelized, pixelize them:
elif [ ! "$res1" = "" ] && [ "$res2" = "" ] ; then
    if [ $res1 -eq -1 ] ; then
	echo >&2 "ERROR: cannot trim a mask pixelized adaptively."
	echo >&2 "Pixelize both your mask and your trimmer polygon(s) using a fixed resolution," 
	echo >&2 "e.g. -Ps0,8, and try again."
	exit 1
    fi
    echo "$MANGLEBINDIR/pixelize -P${scheme1}0,$res1 $trimmer trimmer_pix"
    $MANGLEBINDIR/pixelize -P${scheme1}0,$res1 $trimmer trimmer_pix || exit
    trimmer="trimmer_pix"
elif [ "$res1" = "" ] && [ ! "$res2" = "" ] ; then
    if [ $res2 -eq -1 ] ; then
	echo >&2 "ERROR: cannot trim a mask pixelized adaptively."
	echo >&2 "Pixelize both your mask and your trimmer polygon(s) using a fixed resolution," 
	echo >&2 "e.g. -Ps0,8, and try again."
	exit 1
    fi
    echo "$MANGLEBINDIR/pixelize  -P${scheme2}0,$res2 $mask mask_pix"
    $MANGLEBINDIR/pixelize  -P${scheme2}0,$res2 $mask mask_pix || exit
    mask="mask_pix"
else
    echo "$MANGLEBINDIR/pixelize $pixargs $trimmer trimmer_pix"
    $MANGLEBINDIR/pixelize $pixargs $trimmer trimmer_pix || exit
    echo "$MANGLEBINDIR/pixelize $pixargs $mask mask_pix"
    $MANGLEBINDIR/pixelize $pixargs $mask mask_pix || exit
    trimmer="trimmer_pix"
    mask="mask_pix"
fi

#check if input file is snapped
snapped=`awk '/snapped/{print $1}' < $mask`

#if mask file isn't snapped, snap it
if [ ! "$snapped" = "snapped" ]; then
    echo "Snapping $mask ..."
    mv $mask jp
    $MANGLEBINDIR/snap jp $mask  || exit
    rm jp
fi

#find the complement of the trimmer polygon
echo "$MANGLESCRIPTSDIR/find_complement.sh $trimmer trimmer_comp"
$MANGLESCRIPTSDIR/find_complement.sh $trimmer trimmer_comp || exit

#set the complement to have zero weight
echo 0 > jw0
echo "$MANGLEBINDIR/weight -zjw0 trimmer_comp jcomp0"
$MANGLEBINDIR/weight -zjw0 trimmer_comp jcomp0 || exit
rm jw0

#balkanize the mask with the zero-weight complement, thereby trimming off 
#everything that does not lie within the trimmer polygon 
echo "$MANGLEBINDIR/balkanize $mask jcomp0 jb"
$MANGLEBINDIR/balkanize $mask jcomp0 jb || exit
rm jcomp0

#unify to get rid of zero weight polygons
echo "$MANGLEBINDIR/unify jb $outfile"
$MANGLEBINDIR/unify jb $outfile || exit
rm jb

echo "Polygons of $1 trimmed by $2 written to ${outfile}."

#clean up

rm trimmer_comp
if [ -e trimmer_pix ] ; then 
    rm trimmer_pix
fi
if [ -e mask_pix ] ; then 
    rm mask_pix
fi
