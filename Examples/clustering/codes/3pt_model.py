# =======================================================================
# Example code: how to compute the three-point correlation function model
# =======================================================================

import numpy as np
import matplotlib.pyplot as plt
import CosmoBolognaLib as cbl
from CosmoBolognaLib import DoubleVector as dv
import os

''' Define the cosmology '''
cosmo = cbl.Cosmology(cbl.CosmologicalModel__Planck15_)

''' Compute Pk '''
redshift = 1
kk = np.logspace(-4, 2, 200)
Pk_matter = np.array( cosmo.Pk_matter(kk, "CAMB", False, redshift) )

''' Parameters for 3pt signal '''
rr = dv(np.linspace(1., 300, 200))
theta = np.linspace(0, np.pi, 100)
r1, r2 = 20, 40

''' Slepian model '''
zeta_DM_s = np.array( cosmo.zeta_DM(r1, r2, theta, "Slepian", kk, Pk_matter) )
q_DM_s = np.array( cosmo.Q_DM(r1, r2, theta, "Slepian", kk, Pk_matter) )

''' Barriga-Gatzagnaga model '''
zeta_DM_bg = np.array( cosmo.zeta_DM(r1, r2, theta, "BarrigaGatzanaga", kk, Pk_matter) )
q_DM_bg = np.array( cosmo.Q_DM(r1, r2, theta, "BarrigaGatzanaga", kk, Pk_matter) )

'''Plot the results'''
plt.xlabel(r"$\theta/\pi$")
plt.ylabel(r"$\zeta(r1, r2, \theta)$")
plt.plot(theta/np.pi, zeta_DM_s, '-k', label=r"Slepian et al 2016")
plt.plot(theta/np.pi, zeta_DM_bg, '--r', label=r"Barriga-Gaztanaga 2002")
plt.legend(loc="best")
plt.show(block=False)
plt.xlabel(r"$\theta/\pi$")
plt.ylabel(r"$Q(r1, r2, \theta)$")
plt.plot(theta/np.pi, q_DM_s, '-k', label=r"Slepian et al 2016")
plt.plot(theta/np.pi, q_DM_bg, '--r', label=r"Barriga-Gaztanaga 2002")
plt.legend(loc="best")
plt.show(block=False)

''' NonLocal Term '''
q_DM_NL = np.array( cosmo.Q_nonLocal(r1, r2, theta, kk, Pk_matter) )

plt.xlabel(r"$\theta/\pi$")
plt.ylabel(r"$Q_{NonLocal}(r1, r2, \theta)$")
plt.plot(theta/np.pi, q_DM_NL, '-k')
plt.show(block=False)


