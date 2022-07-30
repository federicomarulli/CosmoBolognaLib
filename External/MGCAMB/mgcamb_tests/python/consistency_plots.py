#
# This file is part of the MGCAMB code
#
# MGCAMB v3.0
#
#   authors:
#       Alex Zucca: azucca@sfu.ca

# This file contains a script that generates a two panel plot for the consistency checks.
# The first panel gives the CMB maximum relative error w.r.t default CAMB, the second
# panel gives the maximum relative error of the matter power spectrum

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab as pil
from itertools import cycle

# font as in latex
from matplotlib import rc
rc('font', **{'size': 16, 'family': 'serif', 'serif': ['Computer Modern Roman']})
rc('text', usetex=True)

# this is the list of grtrans and mnu we wish to consider
grtrans_list    = [0.0001, 0.001, 0.005, 0.01, 0.05, 0.1]
mnu_list        = [0.05,   0.1,   0.15,  0.2,  0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

# lists of errors and quantities for the plots
cls_error_list  = []
cls_stdev_list  = []
cls_max_list    = []
mpk_error_list  = []
mpk_stdev_list  = []
mpk_max_list    = []

#----------------------------------------------
# first generate the arrays of the errors
for mnu in mnu_list:
    #print('doing mnu:', mnu)
    # first, load the GR files
    root_gr     = '../results/results_consistency/LCDM_GR_mnu_0p'+str(mnu)[2:]
    clsGR_file  = root_gr+'_scalCls.dat'
    mpkGR_file  = root_gr+'_matterpower.dat'
    clsGR       = np.loadtxt(clsGR_file)
    mpkGR       = np.loadtxt(mpkGR_file)
    
    cls_error_this_mnu      = []
    cls_stdev_this_mnu      = []
    mpk_error_this_mnu      = []
    mpk_stdev_this_mnu      = []
    cls_max_err_this_mnu    = []
    mpk_max_err_this_mnu    = []
    
    for grtrans in grtrans_list:
        
        root = '../results/results_consistency/LCDM_mg_mnu_0p'+str(mnu)[2:]+'_grt_0p'+str(grtrans)[2:]
        
        cls_file = root + '_scalCls.dat'
        mpk_file = root + '_matterpower.dat'
        
        cls = np.loadtxt(cls_file)
        mpk = np.loadtxt(mpk_file)

        # temperature avg and stdev
        #this_cls_avg   = np.mean(np.abs((cls[:,1] - clsGR[:,1])**2/clsGR[:,1]**2))
        this_cls_avg    = np.mean(np.abs((cls[:,1] - clsGR[:,1])/clsGR[:,1]) )
        this_cls_var    = np.var( np.abs((cls[:,1] - clsGR[:,1])**2/clsGR[:,1]**2))
        this_cls_max    = np.max( np.abs((cls[:,1] - clsGR[:,1])/clsGR[:,1]))
        this_cls_stdev  = np.sqrt(this_cls_var)
        cls_error_this_mnu.append(this_cls_avg)
        cls_stdev_this_mnu.append(this_cls_stdev)
        cls_max_err_this_mnu.append(this_cls_max)

        # matter power spectrum avg and stdev
        #this_mpk_avg   = np.mean(np.abs((mpk[:,1] - mpkGR[:,1])**2/mpkGR[:,1]**2))
        this_mpk_avg    = np.mean(np.abs((mpk[:,1] - mpkGR[:,1])/mpkGR[:,1]))
        this_mpk_var    = np.var(np.abs((mpk[:,1] - mpkGR[:,1]**2)/mpkGR[:,1]**2))
        this_mpk_stdev  = np.sqrt(this_mpk_var)
        this_mpk_max    = np.max(np.abs((mpk[:,1] - mpkGR[:,1])/mpkGR[:,1]))
        mpk_error_this_mnu.append(this_mpk_avg)
        mpk_stdev_this_mnu.append(this_mpk_stdev)
        mpk_max_err_this_mnu.append(this_mpk_max)


    # append the errors and their stdevs
    cls_error_list.append(cls_error_this_mnu)
    cls_stdev_list.append(cls_stdev_this_mnu)
    cls_max_list.append(cls_max_err_this_mnu)
    mpk_error_list.append(mpk_error_this_mnu)
    mpk_stdev_list.append(mpk_stdev_this_mnu)
    mpk_max_list.append(mpk_max_err_this_mnu)



# prepare the figures
fig         = plt.figure(figsize=(16,6))
colors      = ['C'+str(i) for i in range(len(grtrans_list))]
marker_list = ['v', '+', 's', 'o', '^', '*']
markers     = cycle(marker_list)
alphas_list = [0.75]
alphas      = cycle(alphas_list)

# first panel (CMB temperatures fluctuations)
ax1         = plt.subplot(1,2,1)
print('Cls err:')
for (i, m, a) in zip(range(len(grtrans_list)),markers, alphas):
    
    this_error = []
    this_std   = []
    
    for j in range(len(mnu_list)):
        #this_error.append(cls_error_list[j][i])
        if (cls_max_list[j][i] < 1e-7 ):
            this_error.append(1e-7)
        else:
            this_error.append(cls_max_list[j][i])
        this_std.append(cls_stdev_list[j][i])
    
    ax1.semilogy(mnu_list, this_error, marker = m, linestyle = '', color = colors[i], label = r'$a_{\rm trans}$ = '+str(grtrans_list[i]))
    #ax1.semilogy(mnu_list, this_error, marker = m, linestyle = '', color = 'black', label = r'$a_{\rm trans}$ = '+str(grtrans_list[i]), alpha = a)
    print('grtrans',grtrans_list[i])
    for (mnu,err) in zip(mnu_list, this_error):
        print('mnu: {:.5e},  err: {:.5e}'.format(mnu, err))

ax1.legend(bbox_to_anchor=(0.5, -0.25), ncol=3, loc='center left')
ax1.set_xlim(0,0.55)
#.set_ylim(1e-7,1e-4)
ax1.set_xlabel(r'$m_{\nu}$ [eV]')
ax1.set_ylabel(r'$\max | \Delta C_{\ell} / C_{\ell}^{\rm CAMB} | $')

# second panel (matter power spectrum fluctuations)
ax2     = plt.subplot(1,2,2)
markers = cycle(marker_list)
alphas  = cycle(alphas_list)

print('P(k) err:')
for (i, m, a) in zip(range(len(grtrans_list)),markers, alphas):
    
    this_error = []
    this_std   = []
    
    for j in range(len(mnu_list)):
        #this_error.append(mpk_error_list[j][i])
        this_error.append(mpk_max_list[j][i])
        this_std.append(mpk_stdev_list[j][i])
    
    ax2.semilogy(mnu_list, this_error,  marker = m, linestyle = '', color = colors[i], label = r'$a_{\rm trans}$ = '+str(grtrans_list[i]))
    #ax2.semilogy(mnu_list, this_error,  marker = m, linestyle = '', color = 'black', label = r'$a_{\rm trans}$ = '+str(grtrans_list[i]), alpha = a)

    print('grtrans',grtrans_list[i])
    for (mnu,err) in zip(mnu_list, this_error):
        print('mnu: {:.5e},  err: {:.5e}'.format(mnu, err))


ax2.set_xlim(0,0.55)
ax2.set_xlabel(r'$m_{\nu}$ [eV]')
ax2.set_ylabel(r'$ \max | \Delta P(k) / P(k)^{\rm CAMB} |  $')
#ax2.set_ylim(1e-6,1e-2)


# save figure
pil.savefig('cls_mpk_consistency_max.pdf', bbox_inches='tight')
#pil.savefig('cls_mpk_consistency_max.png', bbox_inches='tight')







