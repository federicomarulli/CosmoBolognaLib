# ======================================================
# Example code: how to divide a catalogue in sub-regions
# ======================================================

import CosmoBolognaLib as cbl
from CosmoBolognaLib import StringVector as sv
import numpy as np

# read the catalogue from file
file_cat = sv(["../clustering/input/cat.dat"])

catalogue = cbl.Catalogue(cbl.ObjectType__Galaxy_, cbl.CoordinateType__observed_, file_cat)

# set the regions as angular cells of same area
nCells_Ra = 5
nCells_Dec = 2
cbl.set_ObjectRegion_RaDec(catalogue, nCells_Ra, nCells_Dec);

# plot the angular coordinates of the sample with different colors for
# each region

import matplotlib.pyplot as plt

region_list = np.array(catalogue.region_list())

ra, dec, reg = np.array(catalogue.var(cbl.Var__RA_)),\
               np.array(catalogue.var(cbl.Var__Dec_)),\
               np.array(catalogue.var(cbl.Var__Region_))

for rr in region_list:
    ww = np.where(reg==rr)
    plt.plot(ra[ww], dec[ww], '.')
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$\delta$")
plt.show(block=False)
