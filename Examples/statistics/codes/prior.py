import numpy as np
import matplotlib.pyplot as plt
from CosmoBolognaLib import Prior 
from CosmoBolognaLib import EnumTypes

prior_mean, prior_sigma = 0., 1.
prior_limits = [-1, 1.]
prior = Prior(EnumTypes._GaussianPrior_, [prior_mean, prior_sigma], prior_limits)

nExtr = 2000
sample = np.array([prior.sample(i) for i in range(nExtr)])
psample = np.array([prior(ss) for ss in  sample])

plt.hist(sample, 20, normed=True)
plt.plot(sample, psample, '.')
plt.show()
