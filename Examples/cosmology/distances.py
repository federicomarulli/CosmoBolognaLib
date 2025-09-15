# ==============================================================
# Example code: how to convert redshifts into comoving distances
# ==============================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the Cosmology module of the CosmoBolognaLib
import cblCosmology as cbl

# define a cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# compute the comoving distance at z=1
dc = cosmology.D_C(1)
print('the comoving distance at z=1 is', '%.2f' % dc, 'Mpc/h')

# in Mpc units
cosmology.set_unit(False)
dc = cosmology.D_C(1)
print('the comoving distance at z=1 is', '%.2f' % dc, 'Mpc')
