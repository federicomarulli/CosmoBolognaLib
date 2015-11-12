# =======================================================================
# How to compute power spectrum and two point correlation function models
# =======================================================================


### import Python modules for scientific computing and plotting ###

import numpy as np
import matplotlib.pyplot as plt


### import the CosmoBolognaLib ###

from CosmoBolognaLib import Cosmology


### create an object of class Cosmology ###

cosmo = Cosmology()


### define k and r for the computation of dark matter power spectrum and  two point correlation correlation function ###

kk = np.logspace(-3, 0, 100)
rr = np.linspace(1., 100, 50)


### Compute the power spectrum using CAMB  ###

PkCAMB = np.asarray([cosmo.Pk(kk[i], "CAMB", 0, 0.2) for i in range(len(kk))])


### Compute the two point correlation function using CAMB ###

xiCAMB = [cosmo.xi_DM(rr[i], "CAMB", 0.2, "test", 0) for i in range(len(rr))]


### Plot the results ###

plt.figure(1)
plt.loglog(kk, PkCAMB, label="CAMB")
plt.xlabel(r"$k$")
plt.ylabel(r"$P(k)$")
plt.legend(loc="upper right")

plt.figure(2)
plt.xlabel(r"$r \, \, [\mathrm{Mpc} h^{-1}]$")
plt.ylabel(r"$\xi(r)$")
plt.loglog(rr, xiCAMB, label="CAMB Linear")
plt.legend(loc="upper right")
plt.show()
