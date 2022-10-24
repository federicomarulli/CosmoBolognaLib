# =====================================================================================
# Example code: how to compute power spectrum and two-point correlation function models
# =====================================================================================


# import Python modules for scientific computing and plotting 
import numpy as np
import os
import matplotlib.pyplot as plt

# import the CosmoBolognaLib
import CosmoBolognaLib as cbl

# create an object of class Cosmology
cosmo = cbl.Cosmology(cbl.CosmologicalModel__Planck15_)

# define k and r for the computation of dark matter power spectrum and  two point correlation correlation function
kk = np.logspace(-3, 0, 100)
rr = np.linspace(1., 100, 50)

# compute the power spectrum using CAMB
PkCAMB = cosmo.Pk_matter(kk, "CAMB", False, 0.2)

# compute the two point correlation function using CAMB
xiCAMB = np.asarray([cosmo.xi_matter(rr[i], "CAMB", False, 0.2) for i in range(len(rr))])


# plot the results

plt.figure(1)
plt.loglog(kk, PkCAMB, label="CAMB linear")
plt.xlabel(r"$k \, \, [h\, \mathrm{Mpc}^{-1}]$")
plt.ylabel(r"$P(k) \, \, [(h^{-1}\mathrm{Mpc})^3]$")
plt.legend(loc="upper right")

plt.figure(2)
plt.xlabel(r"$r \, \, [\mathrm{Mpc}\, h^{-1}]$")
plt.ylabel(r"$\xi(r)$")
plt.loglog(rr, xiCAMB, label="CAMB linear")
plt.legend(loc="upper right")

plt.show(block=False)
