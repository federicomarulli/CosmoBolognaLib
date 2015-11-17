# =================================================
# How to measure the two-point correlation function
# =================================================

### import Python modules for scientific computing ###

import os
import numpy as np


### import the CosmoBolognaLib ###

from CosmoBolognaLib import Cosmology
from CosmoBolognaLib import Galaxy
from CosmoBolognaLib import Catalogue
from CosmoBolognaLib import random_catalogue_box
from CosmoBolognaLib import TwoPointCorrelation


### set the cosmological parameters and create the object 'cosmology' ###

cosmology = Cosmology()


### Input/Output files and directories ###
  
file_catalogue = "../input/cat.dat"
dir_output = "../output/"
dir_pairs = dir_output+"pairs/"
dir_random_cat = dir_output

os.system("mkdir -p "+dir_output+" "+dir_pairs)


### read the input catalogue and create the object 'catalogue' ### 

print ("I'm reading the input catalogue...")

ra, dec, redshift = np.genfromtxt(file_catalogue, usecols=(0,1,2), unpack=True)

galaxy = [ Galaxy(ra[i], dec[i], redshift[i], cosmology) for i in range(len(ra))]

catalogue = Catalogue(galaxy)


### construct the random catalogue and create the object 'random_catalogue' ###

N_R = 1.

nRandom = catalogue.nObjects()*N_R

random_catalogue = random_catalogue_box(catalogue, nRandom, dir_random_cat)


### measure the two-point correlation function ###
  
## parameters of the method ##
rMIN = 1.
rMAX = 50.
logbinSize = 0.05
binSize = 0.5
cosSize = 0.02

## create the object used to measure the two-point correlation function ##
TwoP = TwoPointCorrelation(catalogue, random_catalogue)

## set the parameters ##
TwoP.setParameters(rMIN, rMAX, logbinSize, binSize, cosSize)

## measure the two-point correlation function ##
TwoP.measure_xi(dir_pairs, 1)
  
## store the output data ##
TwoP.write_xi(dir_output)

