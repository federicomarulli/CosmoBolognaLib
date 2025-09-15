## =======================================================================
## Example code: how to estimate the matter two-point correlation function
## =======================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")
cosmology.set_parameter("sigma8",0.8111)

# create the object used for the computation of the matter two-point correlation function
PX = cbl.PkXi(cosmology)

print("The matter two-point correlation function at r=50 Mpc/h and z=1 is:", '%.2f' % PX.xi_matter(50, "EisensteinHu", False, 1.))
