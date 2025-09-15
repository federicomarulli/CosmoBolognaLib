## ============================================
## Example code: how to estimate f(z)*sigma8(z) 
## ============================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib modules 
import cblCosmology as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# set the object to manage redshift-space distortion computations
rsd = cbl.RSD(cosmology)

# estimate f*sigma8 at z=1 with CAMB
print("f*sigma8(z=1) =", rsd.fsigma8(1., "EisensteinHu"))
