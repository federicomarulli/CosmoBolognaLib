# =======================================================================
# Example code: how to compute the three-point correlation function model
# =======================================================================

import numpy as np
import matplotlib.pyplot as plt
import CosmoBolognaLib as cbl
from CosmoBolognaLib import DoubleVector as dv
import os

# set the CosmoBolognaLib and the current directories
cbl.SetDirs(os.getcwd()+"/../../../", os.getcwd()+"/")

''' Define the cosmology '''
cosmo = cbl.Cosmology(cbl.CosmologicalModel__Planck15_)

''' Compute Pk '''
redshift = 1
kk = np.logspace(-4, 2, 200)
Pk_DM = np.array( [cosmo.Pk (_kk, "CAMB", False, redshift) for _kk in kk] )

''' Parameters for 3pt signal '''
rr = dv(np.linspace(1., 300, 200))
theta = np.linspace(0, np.pi, 100)
r1, r2 = 20, 40

''' Slepian model '''
zeta_DM_s = np.array( cosmo.zeta_DM(r1, r2, theta, "Slepian", kk, Pk_DM) )
q_DM_s = np.array( cosmo.Q_DM(r1, r2, theta, "Slepian", kk, Pk_DM) )

''' Barriga-Gatzagnaga model '''
zeta_DM_bg = np.array( cosmo.zeta_DM(r1, r2, theta, "BarrigaGatzanaga", kk, Pk_DM) )
q_DM_bg = np.array( cosmo.Q_DM(r1, r2, theta, "BarrigaGatzanaga", kk, Pk_DM) )

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
q_DM_NL = np.array( cosmo.Q_nonLocal(r1, r2, theta, kk, Pk_DM) )

plt.xlabel(r"$\theta/\pi$")
plt.ylabel(r"$Q_{NonLocal}(r1, r2, \theta)$")
plt.plot(theta/np.pi, q_DM_NL, '-k')
plt.show(block=False)


