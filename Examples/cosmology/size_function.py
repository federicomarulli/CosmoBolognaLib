## ====================================================
## Example code: how to estimate the void size function
## ====================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the void size function
sf = cbl.SizeFunction(cosmology)

print("The void size function at R=1 Mpc/h and z=1, predicted by the Sheth and van de Weygaert (2004) model is:", sf.size_function(1., 1., "SvdW", 1.))

