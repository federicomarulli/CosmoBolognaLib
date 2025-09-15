#! /bin/csh -f
#
# =================================
# For HEALPix internal usage only
# =================================
#
# updated 2007-May 15 to deal with 2008Mar07 version
# 2008-03-18: compare files 
# 2009-10-08: deal with 2009Aug20 version
# 2010-05-26: deal with 2010May24 version
# 2012-02-20: deal with 2012Feb01 version
# 2013-02-08: deal with 2013Jan28 version
# 2014-01-06: deal with latest version
# 2016-08-22: detect duplicates
# Note: do not use github because it does not include the Coyote directory
#
# directory with new routines
set dir = /Applications/itt/external_contributions/astron/pro/

# list of routines in current directory
set list = `ls -x *.pro`

set CP = '/bin/cp -p'
#  set CP = 'echo'

set DIFF = '/usr/bin/diff -q'

set nf = 0
set gf = 0
set df = 0
set nf_list = ()
set df_list = ()

# look for duplicates
set duplicates = `find $dir -name \*.pro | xargs basename | sort -n | uniq -d`
set ndup = $#duplicates
if ($ndup > 0) then
    echo
    echo 'Error: found '$ndup' duplicates in '$dir
    foreach file ($duplicates)
	echo ' *****      '$file
	find $dir -name $file -print
    end
    echo 'Aborting!'
    exit
endif

# update routines in current directory
foreach file ($list)
	echo $file
	set fullname = `find $dir -name $file -print`
	if ($#fullname > 1) then
	   echo '-----------------------------'
	   echo 'Duplicate files: '
	   echo $fullname
	   sdiff -s $fullname[1] $fullname[2]
	   echo '-----------------------------'
	   set fullname = $fullname[1]
        endif
	if ($fullname == '') then
	   echo $file' NOT FOUND !'
	   @ nf++
	   set nf_list = ($nf_list $file)
	else
	  #echo $fullname
	  #$DIFF $file $fullname
	  set notsame = `$DIFF $file $fullname | wc -l`
	  if ($notsame == 1) then
	     echo $file':  Files are different!'
	     @ df++
	     set df_list = ($df_list $file)
	  endif
	  $CP $fullname .
	   @ gf++
	endif	  
end
echo
echo
echo 'updated current directory with files from '${dir}
echo
date
echo 'Updated successfully   ' $gf ' files,'
echo 'out of which ' $df ' are different:'
echo $df_list
echo '(beware that some new dependencies may not be satisfied).'
echo 'Not found: ' $nf ' files :'
echo $nf_list

exit

