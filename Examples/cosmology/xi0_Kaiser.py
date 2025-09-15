## ==============================================================================
## Example code: how to compute the z-space matter monopole from the Kaiser model
## ==============================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib and numpy modules
import CosmoBolognaLib as cbl
import numpy as np

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the matter power spectrum
PXS = cbl.PkXizSpace(cosmology)
print("xi0(r=1 Mpc/h, fsigma8=1, bsigma8=1, z=1) =", PXS.xi0_Kaiser(1., 1., 1., "EisensteinHu", 1.))
