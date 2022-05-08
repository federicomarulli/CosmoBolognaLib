import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pil

mg = np.loadtxt('test_LCDM_mg_scalCls.dat')
gr = np.loadtxt('test_LCDM_std_scalCls.dat')

plt.semilogx(gr[:,0], gr[:,1], color = 'C0', label = r'GR')
plt.semilogx(mg[:,0], mg[:,1], color = 'C1', linestyle = '--', label = r'$\mu-\gamma$')

plt.legend(loc='upper right')

plt.show()
