## ===========================================================================
## Example code: how to compute the real-space matter nonlinear power spectrum
## ===========================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib and numpy modules
import CosmoBolognaLib as cbl
import numpy as np

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the matter power spectrum
PXN = cbl.PkXiNonLinear(cosmology)
kk = np.logspace(-4., np.log10(500.), 100)
print("P_delta-delta(k=1e-4, z=1) = ", PXN.Pk_DeltaDelta(kk, 1., "EisensteinHu")[0])
