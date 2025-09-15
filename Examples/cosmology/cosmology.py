## =============================================================
## Example code: how to set and use different cosmologica models
## =============================================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib modules 
import cblCosmology as cbl


### set various cosmological models, with Planck18 parameters ###

# flat Lambda-cold dark matter model
LCDM = cbl.LCDM("Planck18")

# flat Lambda-cold dark matter model with massive neutrinos
LnuCDM = cbl.LnuCDM("Planck18")

# flat w-cold dark matter model
w0CDM = cbl.w0CDM("Planck18")
    
# non-flat Lambda-cold dark matter model
kLCDM = cbl.kLCDM("Planck18")

# flat w-cold dark matter model with massive neutrinos
w0nuCDM = cbl.w0nuCDM("Planck18")
    
# flat w0wa-cold dark matter model
w0waCDM = cbl.w0waCDM("Planck18")
    
# non-flat Lambda-cold dark matter model with massive neutrinos
kLnuCDM = cbl.kLnuCDM("Planck18")

# non-flat w-cold dark matter model
kw0CDM = cbl.kw0CDM("Planck18")

# flat w0wa-cold dark matter model with massive neutrinos
w0wanuCDM = cbl.w0wanuCDM("Planck18")

# non-flat w-cold dark matter model with massive neutrinos
kw0nuCDM = cbl.kw0nuCDM("Planck18")
    
# non-flat w0wa-cold dark matter 
kw0waCDM = cbl.kw0waCDM("Planck18")
    
# non-flat w0wa-cold dark matter model with massive neutrinos
kw0wanuCDM = cbl.kw0wanuCDM("Planck18")

# vector containing all the different cosmological models
cosmology = [LCDM, LnuCDM, w0CDM, kLCDM, w0nuCDM, w0waCDM, kLnuCDM, kw0CDM, w0wanuCDM, kw0nuCDM, kw0waCDM, kw0wanuCDM]


### print the model names and the values of the cosmological model parameters ###

for i in range(len(cosmology)):
    print("\n", cosmology[i].model())
    cosmology[i].print_parameters()
               
