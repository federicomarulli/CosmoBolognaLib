## ===============================================
## Example code: how to estimate the mass variance
## ===============================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the mass variance
SG = cbl.Sigma(cosmology)

print("The mass variance at M=1e15 Msun and z=1 is:", '%.2f' % SG.sigma2M([1.e15], "EisensteinHu", 1.)[0])
