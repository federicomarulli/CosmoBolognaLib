import numpy as np
import matplotlib.pyplot as plt
import CosmoBolognaLib as cbl

xmin = -10.
xmax = 10.
mean = -1.
sigma = 0.1

prior = cbl.PriorDistribution(cbl.DistributionType__Gaussian_, [mean, sigma], xmin, xmax)

nExtr = 2000
sample = np.array([prior.sample(i) for i in range(nExtr)])
psample = np.array([prior(ss) for ss in sample])

plt.hist(sample, 20, normed=True)
plt.plot(sample, psample, '.')
plt.show(block=False)
