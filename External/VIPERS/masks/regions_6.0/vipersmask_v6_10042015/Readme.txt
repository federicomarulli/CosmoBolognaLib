VIPERS pointing masks v6

10/04/2015 added two quadrants W1P125Q3,W1P149Q3
09/04/2015 v6.0 update

contact: Ben Granett ben.granett@brera.inaf.it

~~~ summary ~~~

These are region files describing the VIPERS survey geometry matching
the v6.0 database release.

The pre-imaging with the region polygons overdrawn can be viewed on
these pages:
 http://www.brera.inaf.it/~ben/vipersmask6/index_W1.html
 http://www.brera.inaf.it/~ben/vipersmask6/index_W4.html

What's new for v6: 
* All polygons have been regenerated.  
* In most cases nothing should change from previous releases.  
* Pointings that have been reobserved have been updated.  
* Polygons have been drawn for the lower two rows of W1.

Notes:
There are VIPERS targets outside of these masks but they should
all have redshift flag <= 1.

The 9xx pointings are re-observations.  I have merged the masks so that 
the pointing names are unchanged from previous releases.  The re-observations
filled two previously missing quadrants W1P125Q3 and W1P149Q3.

W1P928 -> W1P128
W1P925 -> W1P125 
W1P949 -> W1P149

W4P916 -> W4P016
W4P925 -> W4P025
W4P929 -> W4P029


The following quadrants have regions with large fraction
of failed targets
  W4P040Q1
  W4P075Q2


~~~ list of files ~~~

pointings_6.0.txt  -- list of pointings and quadrants in v6
vipers_W1_v6.0.reg -- concatenation of polygons for W1
vipers_W4_v6.0.reg -- concatenation of polygons for W4
regions/           -- subdirectory containing region files for every quadrant.
