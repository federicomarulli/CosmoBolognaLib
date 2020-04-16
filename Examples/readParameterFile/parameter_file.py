import CosmoBolognaLib as cbl
from CosmoBolognaLib import StringVector as sv

# Basic - Construct ParameterFile object by reading an existing parameter file
parfile = cbl.ParameterFile("../../External/CAMB/params.ini")

# get key value
output_root = parfile.get_key("output_root", "default_root")
print(output_root)

# set value for existing key
parfile.set_key("output_root", "ciao")

# add a key
parfile["ciao"].push_back("asdf")

# write the parameter file
parfile.write("new_params.ini")

# Advanced - Construct an empty ParameterFile object
parfile = cbl.ParameterFile()

#Add some keys
parfile["redshift"] = ["0.0", "0.1"]
parfile["name"] = ["tessssst"]

# write the parameter file
parfile.write("from_scratch_params.ini")
