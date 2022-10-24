# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

import CosmoBolognaLib as cbl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

cosm = cbl.Cosmology(cbl.CosmologicalModel__Planck18_)
MM = np.logspace(12, 15, 20)
redshift = 0.5

model_MF = "Tinker"
f_R0 = [0., 1.e-4, 1.e-5]
label = ["LCDM", "fR4", "fR5"]
color = ["k", "r", "g"]
ls = ["-", "--", ":"]

fig = plt.figure(figsize=(10,9))
outer_grid = gridspec.GridSpec(1,1, wspace=0.13, hspace=0.1)
inner_grid = gridspec.GridSpecFromSubplotSpec(2,1, subplot_spec=outer_grid[0], height_ratios=[2.5,1.], wspace=0., hspace=0.)

ax = plt.Subplot(fig, inner_grid[0])

MF = []

for n in range(len(f_R0)):
    print ("\n\x1b[0;30;43m"+" I'm computing the mass function the model "+label[n]+"... \x1b[0m")
    MF.append([cosm.mass_function_fR(MM[i], redshift, model_MF, f_R0[n], True) for i in range(len(MM))])
    ax.plot(MM, MF[n], label=label[n], color=color[n], ls=ls[n])

print("\n")

ax.legend(loc='best', fontsize=20)
ax.set_xticklabels([])
ax.get_yaxis().set_tick_params(which='both', direction='in', labelsize=15)
ax.set_ylabel("$\ \mathrm{d}\Phi/\mathrm{d}M\\ [(h^{-1} \ \mathrm{M_{\odot}})^3]$", fontsize=22)
ax.set_xscale('log')
ax.set_yscale('log')
fig.add_subplot(ax)

ax_r = plt.Subplot(fig, inner_grid[1])

ax_r.axhline(y=0, color=color[0])

for n in range(len(f_R0)-1):
    ax_r.plot(MM, (np.array(MF[n+1])-np.array(MF[0]))/MF[0]*100, color=color[n+1], ls=ls[n+1])

ax_r.set_ylim([-80., 80.])
ax_r.set_xscale('log')
ax_r.get_xaxis().set_tick_params(which='both', direction='in', labelsize=15)
ax_r.get_yaxis().set_tick_params(which='both', direction='in', labelsize=15)
ax_r.set_xlabel("$M \ [h^{-1} \ M_\odot]$", fontsize=22)
ax_r.set_ylabel("% diff", fontsize=22)
fig.add_subplot(ax_r)

plt.show(block=False)
