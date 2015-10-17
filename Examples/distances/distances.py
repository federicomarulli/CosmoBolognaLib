# ================================================
# How to convert redshifts into comoving distances
# ================================================

### author of the Python wrapper: Alfonso Veropalumbo ###


### import Python modules for scientific computing and plotting ###

import numpy as np
import matplotlib.pyplot as plt


### import the CosmoBolognaLib ###

from CosmoBolognaLib import pyCosmologyCBL


### set Omega_matter and Omega_DE for three test cosmologies ###

c1 = {"Omega_matter" : 0.2, "Omega_DE" : 0.8}
c2 = {"Omega_matter" : 0.3, "Omega_DE" : 0.7}
c3 = {"Omega_matter" : 0.4, "Omega_DE" : 0.6}


### create three objects of class Cosmology ###

cosmo1 = pyCosmologyCBL(c1)
cosmo2 = pyCosmologyCBL(c2)
cosmo3 = pyCosmologyCBL(c3)


### set the redshifts ###

z = np.linspace(0., 2., 100)


### compute the comoving distances in the three cosmologies ###

dc1 = [ cosmo1.D_C(zz) for zz in z ]
dc2 = [ cosmo2.D_C(zz) for zz in z ]
dc3 = [ cosmo3.D_C(zz) for zz in z ]


### plot the results ###

plt.rc("text", usetex=True)
plt.xlabel(r"$z$")
plt.ylabel(r"$D_C(z)$")
plt.plot(z, dc1, label=r"$\Omega_{M}=0.2$, $\Omega_{DE}=0.8$")
plt.plot(z, dc2, label=r"$\Omega_{M}=0.3$, $\Omega_{DE}=0.7$")
plt.plot(z, dc3, label=r"$\Omega_{M}=0.4$, $\Omega_{DE}=0.6$")
plt.legend(loc="lower right")
plt.show()
