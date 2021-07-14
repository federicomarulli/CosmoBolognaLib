import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pil


gr = np.loadtxt('LCDM_std_mnu0p060_scalCls.dat')
mg1 = np.loadtxt('LCDM_mg_mnu0p060_gt_0p1_scalCls.dat')
mg2 = np.loadtxt('LCDM_mg_mnu0p060_gt_0p01_scalCls.dat')
mg3 = np.loadtxt('LCDM_mg_mnu0p060_gt_0p001_scalCls.dat')
mg4 = np.loadtxt('LCDM_mg_mnu0p060_gt_0p0001_scalCls.dat')

plt.title(r'$\mu,\gamma$, parametrization')
plt.semilogx(mg1[:,0], (mg1[:,1] - gr[:,1])/gr[:,1],linestyle = '-', label = 'GRtrans = 0.1')
plt.semilogx(mg2[:,0], (mg2[:,1] - gr[:,1])/gr[:,1],linestyle = '--', label = 'GRtrans = 0.01')
plt.semilogx(mg3[:,0], (mg3[:,1] - gr[:,1])/gr[:,1],linestyle = '-.', label = 'GRtrans = 0.001')
plt.semilogx(mg4[:,0], (mg4[:,1] - gr[:,1])/gr[:,1],linestyle = ':', label = 'GRtrans = 0.0001')

plt.legend(loc='upper right')

plt.show()
