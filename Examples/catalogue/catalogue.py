import CosmoBolognaLib as cbl
from CosmoBolognaLib import EnumTypes as et
from CosmoBolognaLib import StringVector as sv

file_cat = "cat.dat"
file_cat_vec = sv(1, file_cat)

catalogue = cbl.Catalogue(et._Galaxy_, et._comovingCoordinates_, file_cat_vec)
catalogue2 = cbl.Catalogue (catalogue)

