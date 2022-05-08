# ==============================================================
# Example code: how to convert redshifts into comoving distances
# ==============================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib
import CosmoBolognaLib as cbl

# define a cosmological model, using default parameters
cosm = cbl.Cosmology(cbl.CosmologicalModel__Planck18_)

# compute the comoving distance at z=1
dc = cosm.D_C(1)
print('the comoving distance at z=1 is', '%.2f' % dc, 'Mpc/h')

# in Mpc units
cosm.set_unit(False)
dc = cosm.D_C(1)
print('the comoving distance at z=1 is', '%.2f' % dc, 'Mpc')
