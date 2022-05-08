#
# This file is part of the MGCAMB code
#
# MGCAMB v3.0
#
#   authors:
#       Alex Zucca: azucca@sfu.ca

# This file contains a script that generates a set of CAMB parameters file for testing the GRtransition and the different masses of neutrinos

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab as pil

from matplotlib import rc
rc('font', **{'size': 10, 'family': 'serif', 'serif': ['Computer Modern Roman']})
rc('text', usetex=True)


mnu     = 0.05
grtrans = 0.001
root_gr = '../results/results_consistency/LCDM_GR_mnu_0p'+str(mnu)[2:]
clsGR_file = root_gr+'_scalCls.dat'
mpkGR_file = root_gr+'_matterpower.dat'
clsGR = np.loadtxt(clsGR_file)
mpkGR = np.loadtxt(mpkGR_file)
root = '../results/results_consistency/LCDM_mg_mnu_0p'+str(mnu)[2:]+'_grt_0p'+str(grtrans)[2:]
cls_file = root + '_scalCls.dat'
mpk_file = root + '_matterpower.dat'
cls = np.loadtxt(cls_file)
mpk = np.loadtxt(mpk_file)

delta_mpk = (mpk[:,1] - mpkGR[:,1])/mpkGR[:,1]
plt.semilogx(mpk[:,0], delta_mpk)
plt.xlabel(r'$k$ [h / Mpc]')
plt.ylabel(r'$\Delta P(k) / P(k)^{\rm CAMB}$')
plt.ylim(-2e-3,2e-3)
plt.xlim(1e-4,1)
pil.savefig('mpk_offset.pdf', bbox_inches = 'tight')
#plt.show()







