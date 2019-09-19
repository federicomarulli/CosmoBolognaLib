# ===================================================================
# Example code: how to perform a Bayesian fit to a set of data points
# with a generic model; to run this python script, first compile the
# wrapper with: make modelpy
# ===================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import Python modules for scientific computing 
import os
import numpy as np

# import Python modules for plotting
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# import the CosmoBolognaLib modules 
import CosmoBolognaLib as cbl 

# import the wrapper that returns an object of type Model1D
from modelpy import getModel1D

# set the input/output file/directories 
dir_input = "../input/"
dir_output = "../output/"
file_data = "data.dat"
file_output_start = "model_starting_values.dat"
file_output_bestfit = "model_bestfit.dat"

# construct the dataset by reading an input file 
data = cbl.Data1D(dir_input+file_data)

# set the stuff used to construct the model: here an object of class cosmology, just as an example 
cosmology = cbl.Cosmology()

# set the model to construct the likelihood
model = getModel1D(cosmology)
    
# define a Gaussian likelihood
likelihood = cbl.Likelihood(data, model, cbl.LikelihoodType__Gaussian_Error_)

# set the model parameter
valA, valB, valC = 1., 1., 1. 
start = [valA, valB]

# fix the parameter C, i.e. parameter[2], to valC=1
likelihood.parameters().fix(2, valC) 

# write the model
likelihood.write_model(dir_output, file_output_start, start)

# limits for A and B
minA, maxA = -10., 10.
minB, maxB = -10., 10.
limits = [ [minA, maxA] , [minB, maxB] ]

# maximize the likelihood and write the output
likelihood.maximize(start, limits)
likelihood.write_model_at_bestfit(dir_output, file_output_bestfit)

# construct the priors 
prior_A = cbl.PriorDistribution(cbl.DistributionType__Uniform_, minA, maxA)
prior_B = cbl.PriorDistribution(cbl.DistributionType__Uniform_, minB, maxB)
prior_C = cbl.PriorDistribution(cbl.DistributionType__Constant_, valC)
prior_distributions = cbl.PriorDistributionPtrVector([prior_A, prior_B, prior_C]) # 

# construct the posterior
posterior = cbl.Posterior(prior_distributions, likelihood, 696)

# sample the posterior (starting the MCMC chain from the maximum of the posterior to speed up the chain convergence)
nwalkers = 10;
chain_size = 5000;
posterior.initialize_chains(chain_size, nwalkers, 1.e-5, [valA, valB])
posterior.sample_stretch_move(2)

# show the median MCMC values of the four parameters on screen
print ("\n")
for i in range(posterior.parameters().nparameters()):
    print("Posterior median of %s = %g\n"%(posterior.parameters().name(i), posterior.parameters().bestfit_value(i)))

# show all the MCMC statistics on screen
burn_in = 0
thin = 1
posterior.show_results(burn_in, thin)

# store the chain ouputs
posterior.write_results("../output/", "chains_linear_relation", burn_in, thin)

# store the best-fit model
posterior.write_model_from_chain("../output/", "model_from_chain.dat", [], [], burn_in, thin)

# plot the results
nPar = posterior.parameters().nparameters()

try:
    import corner
except ImportError:
    print("This script requires the corner.py package.")
    print("Please install it if you are interested in producing corner plots")
    print("https://corner.readthedocs.io/en/latest/pages/custom.html")
    pass

# create the corner plot (for simplicity, the fixed parameters are thrown away)

# read the chains
burn_in = 10
thin = 2
table = np.array( [posterior.parameters().parameter_chain_values(i, burn_in, thin) for i in range(nPar)] )
names = [ posterior.parameters().name(i) for i in range(nPar)]
status = [posterior.parameters().status(i) for i in range(nPar)]

# select only non-fixed parameters
from itertools import compress
notFixed = list(map(lambda s : s!="FIXED", status))
notFixedNames = list(compress(names, notFixed))
notFixedChains = np.array(list(compress(table, notFixed))).T

figure = corner.corner(notFixedChains, labels=notFixedNames)
plt.show(block=False)
