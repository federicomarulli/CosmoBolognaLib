## ====================================================
## Example code: how to estimate the Cosmic Mach Number
## ====================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the Cosmic Mach Number
vel = cbl.Velocities(cosmology)

print("The Cosmic Mach Number at r=2 Mpc/h and z=1 is:", '%.2f' % vel.CMN(2., 2., "EisensteinHu", 1.))
