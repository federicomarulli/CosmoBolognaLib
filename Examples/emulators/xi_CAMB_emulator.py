## =======================================================================================
## Example code: how to compute the matter two-point correlation function with an emulator
## =======================================================================================

import numpy as np
import matplotlib.pyplot as plt
import CosmoBolognaLib as cbl

# create a Cosmology object
cosmology = cbl.LCDM("Planck18")
PX = cbl.PkXi(cosmology)

# define the spatial scales and the redshift
rr = np.linspace(30., 120., 100)
zz = 0.

# compute the 2pt correlation function with CAMB
xi_r2_matter = [PX.xi_matter(RR, "CAMB", False, zz)*RR*RR for RR in rr]

# Create an Emulator object
emu = cbl.Emulator("xi_r2_CAMB")

# compute the 2pt correlation function with CAMB
xi_r2_emu = emu.model([cosmology.Omega_matter(), 1-cosmology.Omega_matter(), cosmology.scalar_amp(), zz], rr)

# Plot results
plt.plot(rr, xi_r2_matter, label="CAMB")
plt.plot(rr, xi_r2_emu, label="Emulated CAMB", ls="--")

plt.xlabel(r"$r \ [h^{-1} \ \mathrm{Mpc}]$")
plt.ylabel(r"$\xi(r) \ r^2 \ $")

plt.legend(loc="best")
plt.show(block=False)
