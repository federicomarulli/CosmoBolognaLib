VIPERS photometric masks

contact: ben.granett@brera.inaf.it

History
13/11/2012 -- This version is tagged "Samhain"
12/11/2012 -- added one rectangle to photo_W1_ccd.poly.reg
04/11/2012 -- version 1


What's included: 

1. There are 5 rectangles in W1 of poor photometric quality identified by
high photo-z chi2 and discrepancies between T07 and VIPERS photometry.

2. USNO stars 

The masks around Stars brighter than USNO B=11 were adjusted by hand.
At fainter magnitudes a scaled template was used consisting of a
circle and cross pattern.

For 11 <= B < 17 the following magnitude-radius relations were applied.

  B < 15.19: 
       log10(Radius) = -2.60*log10(B) + 2.33
  B >= 15.19:
       log10(Radius) = -6.55*log10(B) + 6.99

  If there is no B magnitude, R is used with these relations:
  R < 14.28:
       log10(Radius) = -2.52*log10(R) + 2.13
  R >= 14.28:
       log10(Radius) = -6.44*log10(R) + 6.66

The crosses are added using a lookup table, below. The columns are as
follows.  m1 and m2 define the edges of a magnitude bin (B mag).  The
numbers y and x define the lengths of the vertical and horizontal
rectangles in the cross in arcminutes.  The vertical length is y*2 and
the horizontal length is x*2.  The width of the rectangles is 2
arcsec.  The value -1 indicates no cross pattern.

#m1 m2 y x
5 6 -1 -1
6 7 3 3
7 8 4 2
8 9 4 2
9 10 2 2
10 11 2 2
11 12 1 1
12 13 1 1
13 14 .5 .5
14 15 .5 .5
15 17 .3 .3


