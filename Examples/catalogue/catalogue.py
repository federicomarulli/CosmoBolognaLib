# ===================================================================
# Example code: how to construct a catalogue of extragalactic objects
# ===================================================================

### import the CosmoBolognaLib ###
import CosmoBolognaLib as cbl
from CosmoBolognaLib import EnumTypes as et
from CosmoBolognaLib import StringVector as sv

### define a vector of input files (in this case with dim=1) ###
file_cat = "cat.dat"
file_cat_vec = sv(1, file_cat)

### create a catalogue of galaxies, retrieving their comoving coordinates from an input file ###
catalogue = cbl.Catalogue(et._Galaxy_, et._comovingCoordinates_, file_cat_vec)
catalogue2 = cbl.Catalogue(catalogue)

### print the coordinates of the first object ###
print "The coordinates of the first galaxy in the catalogue are:", catalogue[0].xx(), catalogue[0].yy(), catalogue[0].zz()
