# =================================================================================================
# Example code: how to measure the angle-averaged two-point correlation function, i.e. the monopole 
# =================================================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import Python modules for scientific computing 
import os
import numpy as np

# import the CosmoBolognaLib modules 
import CosmoBolognaLib as cbl 

# define the cosmological model, with default parameters 
cosmology = cbl.Cosmology(cbl.CosmologicalModel__Planck15_)

# Input/Output files and directories 
file_catalogue = ("../input/cat.dat",)
dir_output = "../output/"
dir_pairs = dir_output+"pairs/"
dir_random_cat = dir_output
os.system("mkdir -p "+dir_output+" "+dir_pairs)

# read the input galaxy catalogue (with polar coordinates: RA, Dec, redshift)  
print("I'm reading the input catalogue...")
catalogue = cbl.Catalogue(cbl.ObjectType__Galaxy_, cbl.CoordinateType__observed_, file_catalogue, cosmology)

# construct the random catalogue (with cubic geometry) 
print("I'm creating the catalogue...")
N_R = 1 # random/object ratio
random_catalogue = cbl.Catalogue(cbl.RandomType__createRandom_box_, catalogue, N_R)


# measure the monopole of the two-point correlation function 
  
# binnig parameters #
rMin = 5.   # minimum separation 
rMax = 20.  # maximum separation 
nbins = 5  # number of bins
shift = 0.5 # spatial shift used to set the bin centre 

# create the object used to measure the two-point correlation function
TwoP = cbl.TwoPointCorrelation1D_monopole(catalogue, random_catalogue, cbl.BinType__logarithmic_, rMin, rMax, nbins, shift)

# measure the two-point correlation function
TwoP.measure(cbl.ErrorType__Poisson_, dir_pairs)

# store the output data
file_xi = "xi.dat"
TwoP.write(dir_output, file_xi)

