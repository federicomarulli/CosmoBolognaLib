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
import pandas as pd
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
    valA, valB, valC = 1., 2., 1. 

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
    
    # store the chain ouputs
    posterior.write_results("../output/", results_name)
    
    return posterior



def plot_contours (posterior, burn_in, thin, figure, axes, color, label):

    # read the chains
    nPar = posterior.parameters().nparameters()
    table = np.array( [posterior.parameters().parameter_chain_values(i, burn_in, thin) for i in range(nPar)] )
    names = [posterior.parameters().name(i) for i in range(nPar)]
    status = [posterior.parameters().status(i) for i in range(nPar)]
    weights = posterior.weight(burn_in, thin)

    # select only non-fixed parameters
    from itertools import compress
    notFixed = list(map(lambda s : s!="FIXED", status))
    notFixedNames = list(compress(names, notFixed))
    notFixedChains = np.array(list(compress(table, notFixed))).T

    corner.corner(notFixedChains, \
                  weights=np.array(weights),\
                  labels=notFixedNames,\
                  truths=(1, 2., 4.), truth_color='black', \
                  plot_datapoints=False, show_titles=True, smooth=True, plot_density=False, levels=(1-np.exp(-0.5), 1-np.exp(-2), 1-np.exp(-4.5)),\
                  color=color, fill_contours=True, hist_kwargs={"density" : True, "color":color}, fig=figure)

    axes[0][1].plot(1,1,color, label=label)
    axes[0][1].legend(loc="best", fontsize=11)



# function that joines two chains: one after the other
def join_probes(posterior,dir_input,file_data_1, file_data_2, nwalkers):

    file_name_1 = dir_input+file_data_1
    file_name_2 = dir_input+file_data_2
    imp_chain_1 = np.loadtxt(file_name_1)
    imp_chain_2 = np.loadtxt(file_name_2)
    samples_1 = pd.DataFrame(imp_chain_1[:,:])
    samples_2 = pd.DataFrame(imp_chain_2[:,:])
    joined_samples = pd.concat([samples_1,samples_2])
    joined_file_name = "joined_impsamp.dat"
    open(dir_input+joined_file_name,"w").close()
    joined_samples.to_csv(dir_input+joined_file_name,header=None,index=None,sep='\t',mode="w")

    joined_posterior = posterior
    joined_posterior.read_chain(dir_input, joined_file_name,nwalkers)

    return joined_posterior



# set parameters for posterior sampling
nwalkers = 100
chain_size = 500
burn_in = 0
thin = 1

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

# do the importance sampling and plot the contours

fig, axes = plt.subplots(3, 3, figsize=(15, 15))

plot_contours(posterior1, burn_in, thin, fig, axes, "b", "Posterior 1")
plot_contours(posterior2, burn_in, thin, fig, axes, "y", "Posterior 2")

# do the importance sampling
posterior1.importance_sampling("../output/", "model2_chain.dat")

# store the chain ouputs
posterior1.write_results("../output/", "model_1+2_importance_sampling")

# plot the contours
plot_contours(posterior1, burn_in, thin, fig, axes, "g", "imp. sampling 1+2")

# do the importance sampling
posterior2.importance_sampling("../output/", "model1_chain.dat")

# store the chain outputs
posterior2.write_results("../output/", "model_2+1_importance_sampling")

# plot the contours
plot_contours(posterior2, burn_in, thin, fig, axes, "r", "imp. sampling 2+1")

file_impsamp1 = "model_1+2_importance_sampling_chain.dat"
file_impsamp2 = "model_2+1_importance_sampling_chain.dat"

# join chains from importance sampling
joined_posterior = join_probes(posterior1, dir_output, file_impsamp1, file_impsamp2, 1)

# plot the contours
plot_contours(joined_posterior,burn_in, thin, fig, axes, "k", "imp. sampling concat.")

plt.show()
