# =================================================================
# Example code: how to create a mesh grid of a catalogue of objects
# =================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib 
import CosmoBolognaLib as cbl

# load the input object catalogue
file_cat = "../input/halo_catalogue.txt"
var_names_tracers = [cbl.Var__X_, cbl.Var__Y_, cbl.Var__Z_]
columns_tracers = [1, 2, 3]
tracers_catalogue = cbl.Catalogue(cbl.ObjectType__Halo_, cbl.CoordinateType__comoving_, cbl.VarCast(var_names_tracers), columns_tracers, [file_cat], 0)

# store the mean particle separation of the catalogue
mps = tracers_catalogue.mps()
        
# generate the mesh grid
MG = cbl.MeshGrid(tracers_catalogue.var(cbl.Var__X_), tracers_catalogue.var(cbl.Var__Y_), tracers_catalogue.var(cbl.Var__Z_), 2*mps, [])

# test the closeObjects function, which returns the indices of the objects within a specified distance range from a point
print("Clos object test:")
Xcoord = 450.
Ycoord = 500.
Zcoord = 500.
closeObj = MG.closeObjects(Xcoord, Ycoord, Zcoord, 15.0*mps, 12.5*mps)
print("Number of objects between 12.5 and 15 MPS: " + str(len(list(closeObj))))
