## ====================================================
## Example code: how to estimate the halo mass function
## ====================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the halo mass function
mf = cbl.MassFunction(cosmology)

print("The halo mass function at M=1e13 Msun/h and z=1, predicted by the Press & Schechter model is:", mf.mass_function(1.e13, 1., "PS", "EisensteinHu"))
