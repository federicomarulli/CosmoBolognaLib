# import python modules
import numpy as np
import os
import matplotlib.pyplot as plt

# import the CosmoBolognaLib 
import CosmoBolognaLib as cbl 
from CosmoBolognaLib import DoubleVector as dv

# set the cosmological model, with default parameters
cosmology = cbl.Cosmology()

# compute the dark matter power spectrum
kk = np.logspace(-4, 2, 200)
Pk = cosmology.Pk_matter(kk, "CAMB", False, 0.2)

# get correlation function from fftlog: dir is the transformation
# direction, mu is the order of the Bessel function (see the
# documentation for other options)
dir = 1
mu = 0
rr = np.linspace(1., 200, 100)
xi = np.array(cbl.transform_FFTlog(dv(rr), dir, dv(kk), dv(Pk), mu))

# plot results
plt.plot(rr, xi*rr*rr)
plt.xlabel(r"$s$ $[$Mpc$h^{-1}]$")
plt.ylabel(r"$\xi(s)\cdot s^2$ $[$Mpc$^2h^{-2}]$")
plt.plot(rr, xi*rr*rr, '-')

plt.show(block=False)
