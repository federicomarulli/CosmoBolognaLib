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

try:
    import corner
except ImportError:
    print("This script requires the corner.py package.")
    print("Please install it if you are interested in producing corner plots")
    print("https://corner.readthedocs.io/en/latest/pages/custom.html")
    pass


# set the input/output file/directories 
dir_input = "../input/"
dir_output = "../output/"
file_data1 = "data1.dat"
file_data2 = "data2.dat"

def go (data, model, nwalkers, chain_size, results_name, model_name):

    # define a Gaussian likelihood
    likelihood = cbl.Likelihood(data, model, cbl.LikelihoodType__Gaussian_Error_)

    # set the model parameter
    valA, valB, valC = 1., 1., 1. 

    # limits for A and B
    minA, maxA = -10., 10.
    minB, maxB = -10., 10.

    # construct the priors 
    prior_A = cbl.PriorDistribution(cbl.DistributionType__Uniform_, minA, maxA)
    prior_B = cbl.PriorDistribution(cbl.DistributionType__Uniform_, minB, maxB)
    prior_C = cbl.PriorDistribution(cbl.DistributionType__Constant_, valC)
    prior_distributions = cbl.PriorDistributionPtrVector([prior_A, prior_B, prior_C])  

    # construct the posterior
    posterior = cbl.Posterior(prior_distributions, likelihood, 696)

    # sample the posterior (starting the MCMC chain from the maximum of the posterior to speed up the chain convergence)
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
    posterior.write_results("../output/", results_name, burn_in, thin)

    # store the best-fit model
    posterior.write_model_from_chain("../output/", model_name, [], [], burn_in, thin)

    return posterior

def plot_contours(posterior, burn_in, thin, figure, color):

    # read the chains
    nPar = posterior.parameters().nparameters()
    table = np.array( [posterior.parameters().parameter_chain_values(i, burn_in, thin) for i in range(nPar)] )
    names = [ posterior.parameters().name(i) for i in range(nPar)]
    status = [posterior.parameters().status(i) for i in range(nPar)]
    weights = posterior.weight(burn_in, thin)

    # select only non-fixed parameters
    from itertools import compress
    notFixed = list(map(lambda s : s!="FIXED", status))
    notFixedNames = list(compress(names, notFixed))
    notFixedChains = np.array(list(compress(table, notFixed))).T

    corner.corner(notFixedChains, \
                  weights=1./np.array(weights),\
                  labels=notFixedNames,\
                  truths=(1, 2, 1), truth_color='black', \
                  plot_datapoints=False, plot_density=False, levels=(1-np.exp(-0.5), 1-np.exp(-2), 1-np.exp(-4.5)),\
                  color=color, fill_contours=True, hist_kwargs={"density" : True, "color":color}, fig=figure)


# Set parameters for posterior sampling
nwalkers = 100
chain_size = 5000
burn_in = 100
thin = 10

# set the stuff used to construct the model: here an object of class cosmology, just as an example 
cosmology = cbl.Cosmology()
   
# construct the dataset by reading an input file 
data1 = cbl.Data1D(dir_input+file_data1)

# set the model to construct the likelihood
model1 = getModel1D(cosmology, 0)
posterior1 = go(data1, model1, nwalkers, chain_size, "model1", "model1.dat")

# construct the dataset by reading an input file 
data2 = cbl.Data1D(dir_input+file_data2)

# set the model to construct the likelihood
model2 = getModel1D(cosmology, 1)
posterior2 = go(data2, model2, nwalkers, chain_size, "model2", "model2.dat")

#Plot the contours

fig, axes = plt.subplots(3, 3, figsize=(15, 15))

plot_contours(posterior2, burn_in, thin, fig, "r")
plot_contours(posterior1, burn_in, thin, fig, "b")

#Do the importance sampling
posterior1.importance_sampling("../output/", "model2_chain.dat", nwalkers)

# store the chain ouputs
posterior1.write_results("../output/", "model_1+2_importance_sampling")

plot_contours(posterior1, burn_in, thin, fig, "g")

plt.show()
