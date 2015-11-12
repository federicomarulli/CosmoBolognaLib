# ================================================
# How to convert redshifts into comoving distances
# ================================================

### import Python modules for scientific computing and plotting ###

import numpy as np
import matplotlib.pyplot as plt


### import the CosmoBolognaLib ###

from CosmoBolognaLib import Cosmology 


### create an object of class Cosmology with default cosmological parameters ###

cosm = Cosmology()


### set the redshifts ###

redshift = np.linspace(0., 2., 100)


### compute the comoving distances ###

dc = [ cosm.D_C(zz) for zz in redshift ]


### plot the results ###

plt.rc("text", usetex=True)
plt.xlabel(r"$z$")
plt.ylabel(r"$D_C(z)$")
plt.plot(redshift, dc, label="$\Omega_{M}=$"+str(cosm.Omega_matter())+" $\Omega_{DE}=$"+str(cosm.Omega_DE()))
plt.legend(loc="lower right")
plt.show()
