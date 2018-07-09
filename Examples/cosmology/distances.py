# ==============================================================
# Example code: how to convert redshifts into comoving distances
# ==============================================================


### import cosmological functions ###

import CosmoBolognaLib as cbl


### define a cosmological model, using default parameters ###

cosm = cbl.Cosmology(cbl.CosmologicalModel__Planck15_)


### compute the comoving distance at z=1 ###

dc = cosm.D_C(1)

print 'the comoving distance at z=1 is', '%.2f' % dc, 'Mpc/h'

