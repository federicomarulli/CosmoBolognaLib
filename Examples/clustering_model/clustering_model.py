# =======================================================================
# How to compute power spectrum and two point correlation function models
# =======================================================================

### author of the Python wrapper: Alfonso Veropalumbo ###

### import Python modules for scientific computing and plotting ###

import numpy as np
import matplotlib.pyplot as plt

### import the CosmoBolognaLib ###

from CosmoBolognaLib import pyCosmologyCBL

### set Omega_matter and Omega_DE for three test cosmologies ###

cosmopar = { "Omega_matter" : 0.3, "Omega_DE" : 1-0.3, "hh" : 0.68}  


### create an object of class Cosmology ###

cosmo = pyCosmologyCBL(cosmopar)


### define the sigma8 parameter ###

cosmo.set_sigma8(0.8)

### define k and r for the computation of dark matter power spectrum and  two point correlation correlation function ###

kk = np.logspace(-3,0,100)
r = np.linspace(1.,200,50)

### Compute the power spectrum using CAMB and EisensteinHu method ###

PkCAMB = np.asarray([cosmo.Pk(kk[i],"CAMB",0,0.2) for i in range(len(kk))])
PkEH = np.asarray([cosmo.Pk(kk[i],"EisensteinHu",0,0.2) for i in range(len(kk))])

### Compute the two point correlation function using CAMB and EisensteinHu method ###

xiEH = [cosmo.xi_DM(r[i],"EisensteinHu",0.2,NL=0) for i in range(len(r))]
xiCAMB = [cosmo.xi_DM(r[i],"CAMB",0.2,NL=0) for i in range(len(r))]
xiCAMBNL = [cosmo.xi_DM(r[i],"CAMB",0.2) for i in range(len(r))]


### Plot the results ###
plt.loglog(kk,PkEH,label="EisensteinHu")
plt.loglog(kk,PkCAMB,label="CAMB")
plt.xlabel(r"$k$")
plt.ylabel(r"$P(k)$")
plt.legend(loc="upper right")
plt.show()


plt.xlabel(r"$r \, \, [\mathrm{Mpc} h^{-1}]$")
plt.ylabel(r"$\xi(r)$")
plt.loglog(r,xiEH,label="EisensteinHu")
plt.loglog(r,xiCAMB,label="CAMB Linear")
plt.loglog(r,xiCAMBNL,label="CAMB NonLinear")
plt.legend(loc="upper right")
plt.show()
