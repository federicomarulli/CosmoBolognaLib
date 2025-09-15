## ===========================================
## Example code: how to estimate the halo bias
## ===========================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the bias 
bias = cbl.Bias(cosmology)

print("The halo bias at M=1e13 Msun/h and z=1, predicted by the Sheth & Tormen (1999) model, is:", bias.bias_halo(1.e13, 1., "ST99", "EisensteinHu"))


