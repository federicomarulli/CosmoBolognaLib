# ===================================================================
# Example code: how to construct a catalogue of extragalactic objects
# ===================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib 
import cblCatalogue as cbl
from CosmoBolognaLib import StringVector as sv

# define a std::vector of input files (in this case with dim=1) 
file_cat = "cat.dat"
file_cat_vec = sv(1, file_cat)

# create a catalogue of galaxies, retrieving their comoving coordinates from an input file 
catalogue = cbl.Catalogue(cbl.ObjectType__Galaxy_, cbl.CoordinateType__comoving_, file_cat_vec)
catalogue2 = cbl.Catalogue(catalogue)

# print the coordinates of the first object 
print("The coordinates of the first galaxy in the catalogue are:", catalogue[0].xx(), catalogue[0].yy(), catalogue[0].zz())
