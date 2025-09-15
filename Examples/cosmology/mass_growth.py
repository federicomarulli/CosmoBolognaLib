## ===========================================================================
## Example code: how to compute the cumulative formation redshift distribution
## ===========================================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib module
import CosmoBolognaLib as cbl

# set the cosmological model, using default parameters
cosmology = cbl.LCDM("Planck18")

# create the object used for the computation of the formation redshift distribution
mg = cbl.MassGrowth(cosmology)

ww = 1.	       # rescaled variable w as in Lacey and Coles (1993)
ff = 0.5       # assembled fraction
author = "GTS" # Giocoli et al. (2012)

print("The cumulative formation redshift distribution at w=1 and f=0.5, predicted by Giocoli et al. (2012) model is:", mg.cumPw(ww, ff, author))
