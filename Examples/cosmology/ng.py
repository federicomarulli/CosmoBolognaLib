## ==========================================
## Example code: how to estimate the skewness
## ==========================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the skewness
ng = cbl.NG(cosmology)

print("The skewness at M=1e13 Msun/h is:", ng.skewness(1.e13, "EisensteinHu"))
