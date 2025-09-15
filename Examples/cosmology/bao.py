## =============================================================
## Example code: how to estimate the sound horizon at drag epoch
## =============================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the sound horizon at drag epoch
bao = cbl.BAO(cosmology)

print("The sound horizon at drag epoch is:", bao.rs("EisensteinHu"))

