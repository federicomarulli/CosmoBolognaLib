# ============================================
# Example code: how to manage a parameter file
# ============================================

# to ensure compatibility in Python versions 2.x and 3.x
from __future__ import print_function

# import the CosmoBolognaLib modules 
import CosmoBolognaLib as cbl
from CosmoBolognaLib import StringVector as sv


### How to modify an existing parameter file ###

# read an existing parameter file
parfile = cbl.ParameterFile("../../External/CAMB/inifiles/params.ini")

# get a key value
output_root = parfile.get_key("output_root", "default_root")
print("output_root =",output_root)

# set the value for the existing output_root key
parfile.set_key("output_root", "tessssst")

# add a new key
parfile["my_parameter"].push_back("666")

# write the new parameter file
parfile.write("new_params.ini")


### How to create a new parameter file ###

# construct an empty ParameterFile object
parfile = cbl.ParameterFile()

# add some keys
parfile["redshift"] = ["0.0", "0.1"]
parfile["name"] = ["tessssst"]

# write the parameter file
parfile.write("params_from_scratch.ini")
