# =================================================================
# Example code: how to measure the three-point correlation function
# =================================================================

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
HOME = os.getenv("HOME")
file_catalogue = (HOME+"/CosmoBolognaLib/Examples/clustering/input/cat.dat",)
dir_output = HOME+"/CosmoBolognaLib/Examples/clustering/output/"
dir_triplets = dir_output+"triplets/"
dir_random_cat = dir_output
os.system("mkdir -p "+dir_output+" "+dir_triplets)

# read the input galaxy catalogue (with polar coordinates: RA, Dec, redshift)  
print("I'm reading the input catalogue...")
catalogue = cbl.Catalogue(cbl.ObjectType__Galaxy_, cbl.CoordinateType__observed_, file_catalogue, cosmology)

# construct the random catalogue (with cubic geometry) 
print("I'm creating the catalogue...")
N_R = 2 # random/object ratio
random_catalogue = cbl.Catalogue(cbl.RandomType__createRandom_box_, catalogue, N_R)

# construct the sub-regions used for jackknife and bootstrap
print("I'm constructing the sub-regions used for jackknife and bootstrap...")
nx = 3
ny = 3
nz = 3
cbl.set_ObjectRegion_SubBoxes(catalogue, random_catalogue, nx, ny, nz)


# measure the three-point correlation function 
  
# binnig parameters 
side_s = 20.    # 1st side of the triangle
side_u = 2.     # ratio between the 1st and 2nd sides of the triangle (u*s)
perc = 0.0225   # tolerance
nbins = 15      # number of bins

# create the object used to measure the three-point correlation function #
ThreeP = cbl.ThreePointCorrelation_comoving_reduced(catalogue, random_catalogue, cbl.TripletType__comoving_theta_, side_s, side_u, perc, nbins)

# measure the two-point correlation function 
ThreeP.measure(cbl.ErrorType__Jackknife_,dir_triplets, dir_output)

# store the output data 
file_zeta = "zeta_JK.dat"
ThreeP.write(dir_output, file_zeta, True)
file_Q_cov = "Q_cov.dat"
ThreeP.write_covariance(dir_output, file_Q_cov)
