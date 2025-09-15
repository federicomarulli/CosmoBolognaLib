## ======================================================
## Example code: how to estimate the function Gamma_{123}
## ======================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the Cosmology module of the CosmoBolognaLib
import cblCosmology as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the function Gamma_{123}
threept = cbl.ThreePointCorrelation(cosmology)

print("The non-local contribution to the three-point correlation function at r1=r2=1 Mpc/h is:", threept.Gamma_3PCF(1., 1., 1., [1., 1.], [1., 1.]))
