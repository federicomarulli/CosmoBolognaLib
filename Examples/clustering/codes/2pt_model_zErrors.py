# =======================================================================================================================
# Example code: how to compute the two-point correlation function model taking into account the effect of redshift errors
# =======================================================================================================================

# import the CosmoBolognaLib and the other required Python modules
import CosmoBolognaLib as cbl
from CosmoBolognaLib import DoubleVector as dv
import numpy as np
import os
import matplotlib.pyplot as plt

# set the cosmological model 
cosmology = cbl.Cosmology(cbl.CosmologicalModel__Planck18_)

# compute the linear growth rate at z=1
redshift = 1.
linear_growth_rate = cosmology.linear_growth_rate(redshift, 1.)

# set the galaxy bias
bias = 2

# set the wave vector and comoving scale vectors
kk = np.logspace(-4, 2, 500)
rr = np.linspace(1, 150, 100)

# compute the matter power spectrum
Pk = cosmology.Pk_matter(kk, "CAMB", False, redshift)
interpPk = cbl.FuncGrid(dv(kk), dv(Pk), "Spline")

# set the redshift errors
SigmaZ = [1.e-4, 1.e-3, 1.e-2, 1.e-1]


# plot the results

figure = plt.figure()
ax = figure.add_subplot(111)
ax.set_xlabel(r"$s [\mathrm{Mpc}/h]$")
ax.set_ylabel(r"$\xi(s)$")

ax.set_xscale("log")
ax.set_yscale("log")

for ss in SigmaZ:
   SigmaS = cbl.cc*ss*(1.+redshift)/cosmology.HH(redshift)
   xi = cbl.damped_Xi(dv(rr), bias, linear_growth_rate, SigmaS, dv(kk), interpPk)

   ax.plot(rr, xi, label=r"$\sigma_z = %g$"%ss)

ax.legend(loc="best")

plt.show(block=False) # set block=True to see the plot
