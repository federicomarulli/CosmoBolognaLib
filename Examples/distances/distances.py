# ==============================================================
# Example code: how to convert redshifts into comoving distances
# ==============================================================


### import cosmological functions ###

from CosmoBolognaLib import Cosmology 


### define a cosmological model, using default parameters ###

cosm = Cosmology()


### compute the comoving distance at z=1 ###

dc = cosm.D_C(1)

print 'the comoving distance at z=1 is', '%.2f' % dc, 'Mpc/h'

