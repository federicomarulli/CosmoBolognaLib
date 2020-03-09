# ==========================================
# Example code: how to subsample a catalogue
# ==========================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib 
import CosmoBolognaLib as cbl
from CosmoBolognaLib import StringVector as sv

import time

# define a std::vector of input files (in this case with dim=1) 
file_cat = "cat.dat"
file_cat_vec = sv(1, file_cat)

# create a catalogue of galaxies, retrieving their comoving coordinates from an input file 
catalogue = cbl.Catalogue(cbl.ObjectType__Galaxy_, cbl.CoordinateType__comoving_, file_cat_vec)

# Define a mask by creating a new class that inherits from cbl object MaskObject
# This must have a constructor (__init__) and a callable function (__call__)
# The callable method must take in input an object and return a boolean
class PYMaskObject(cbl.MaskObject):

    def __init__(self):
        cbl.MaskObject.__init__(self)

    def __call__(self, obj):
        return (obj.xx()>-70) & (obj.xx()<10)

mask = PYMaskObject()

# Apply the mask
print("Start")
start = time.time()
cat2 = catalogue.sub_catalogue(mask)
stop = time.time()
print(stop-start)
print("Stop")

# Verify the mask works by comparing with a simpler version of the same operation
cat3 = catalogue.sub_catalogue(cbl.Var__X_, -70, 10)
print("Full sample: ", catalogue.nObjects())
print("User-defined mask: ", cat2.nObjects())
print("Min-Max cut:  ", cat3.nObjects())


