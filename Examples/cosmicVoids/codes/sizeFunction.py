# ==========================================================================
# Example code: how to compute the theoretical size function of cosmic voids
# ==========================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib #
import CosmoBolognaLib as cbl
import numpy as np

# set the CosmoBolognaLib and the current directories
cbl.SetDirs("../../../", "./")

# define a cosmological model, using default parameters #
cosm = cbl.Cosmology()

# Minimum and maximum of effective void radii
R_min = 1.
R_max = 30.

# number of radii at which the size function is computed
n_val = 10

# list of effective void radii with logarithmic binning
RR = np.logspace(np.log10(R_min), np.log10(R_max), n_val , endpoint=True)

# redshift of the sample
zz = 0.

# effective bias of the mass tracers
b_eff = 1.

print('The size function at z = %.2g'%zz ,' is :')

# compute the  Sheth & van de Weygaert size function 
for i in range(n_val):
    print('%.6e' %cosm.size_function(RR[i], zz, "SvdW", b_eff), '(h/Mpc)^3 at R = %.3g' %RR[i], ' Mpc/h ')
