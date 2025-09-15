#!/usr/bin/python3

import os
import sys


####################################################

### function to check the output of os.system ###

def check_sys (out):
    if (out>0):
        sys.exit(0)

        
### function to check if valgrind got any errors ###

def check_valgrind ():
    if "ERROR SUMMARY: 0 errors from 0" in open("valgrind_output.dat").read():
        print("\n valgrind: OK! \n")
        os.remove("valgrind_output.dat")
    else:
        print("\n Error from vagrind! \n")
        print("".join(open("valgrind_output.dat")))
        os.remove("valgrind_output.dat")
        sys.exit(0)

        
### function to run the example code (with or withoud valgrind) ###

flags_Valgrind = "--track-origins=yes" #"--leak-check=full"
error_msg = "\n Error: please select a valid programming language among C++, python and python2! Read the file Instructions.txt for more details \n"

def com (ex, language, input):
    if not ("valgrind" in sys.argv):
        if (language=="C++"):
            out = os.system("./"+ex+"  "+input)
            check_sys(out)
        elif (language=="python2"):
            out = os.system("python2 "+ex+"  "+input)
            check_sys(out)
        elif (language=="python"):
            out = os.system("python "+ex+"  "+input)
            check_sys(out)
        else:
            print(error_msg)
            sys.exit(0)
    else:
        if (language=="C++"):
            out = os.system("valgrind "+flags_Valgrind+" --log-file=\"valgrind_output.dat\" ./"+ex+"  "+input)
            check_sys(out)
            check_valgrind()
        elif (language=="python2"):
            out = os.system("python2 "+ex+"  "+input)
            check_sys(out)
        elif (language=="python"):
            out = os.system("python "+ex+"  "+input)
            check_sys(out)             
        else:
            print(error_msg)
            sys.exit(0)

        
### function to check if the example is ok ###

def check (dir, file, language, input=""):
    if (("ALL" in sys.argv and language!="python2") or (file in sys.argv and "C++" in sys.argv and language=="C++") or (file in sys.argv and "python2" in sys.argv and language=="python2") or (file in sys.argv and "python" in sys.argv and language=="python") or ("all" in sys.argv and "python2" in sys.argv and language=="python2") or ("all" in sys.argv and "python" in sys.argv and language=="python") or ("all" in sys.argv and "C++" in sys.argv and language=="C++")):
        print("\n\n-----> Testing", file, "<-----\n")
        os.chdir(cwd+"/Examples/"+dir)
        com(file, language, input)

        
####################################################


os.system("clear")
cwd = os.getcwd()


### compile the CBL ###

if not ("nocompile" in sys.argv):
    comm = "make purgeALL && make FLAGS=\"-O0 -g\" && "
    if ("python" in sys.argv):
        comm = comm+"make allExamples PY=python FLAGS=\"-O0 -g\""
    else:
        comm = comm+"make allExamples PY=python2 FLAGS=\"-O0 -g\""
    out = os.system(comm)
    check_sys(out)


####################################################


### check the C++ examples -> check(directory, example_executable, language, input) ###

check("vectors", "vectors", "C++") 

check("eigen", "eigen", "C++") 

check("randomNumbers", "randomNumbers", "C++") 
check("randomNumbers", "randomNumbers_custom", "C++") 
check("randomNumbers", "correlated_samples", "C++") 

check("histogram", "histogram", "C++") 

check("wrappers", "integration_gsl", "C++")
check("wrappers", "minimisation_gsl", "C++") 
#check("wrappers", "integration_cuba", "C++") # seg fault with valgrind!!!
check("wrappers", "fits", "C++") 

check("covsample", "covsample", "C++") 

check("cosmology", "LCDM", "C++") 
check("cosmology", "cosmology", "C++") 
check("cosmology", "distances", "C++") 
#check("cosmology", "model_cosmology", "C++") 
check("cosmology", "growth_rate", "C++")
check("cosmology", "Pk_BoltzmannSolver", "C++")
#check("cosmology", "velocities", "C++")
check("cosmology", "three_point_correlation", "C++")
check("cosmology", "size_function", "C++")
#check("cosmology", "ng", "C++")
check("cosmology", "bias", "C++")
check("cosmology", "mass_function", "C++")
#check("cosmology", "fsigma8", "C++")
check("cosmology", "bao", "C++")
check("cosmology", "mass_growth", "C++")
check("cosmology", "Pk_DeltaDelta", "C++")
check("cosmology", "xi0_Kaiser", "C++")
check("cosmology", "sigma", "C++")
check("cosmology", "xi_matter", "C++")

check("data", "data1D", "C++")

check("statistics/codes", "prior", "C++")
#check("statistics/codes", "fit", "C++") 
#check("statistics/codes", "sampler", "C++")

check("meshGrid/codes", "meshGridTest", "C++") 

check("catalogue", "catalogue", "C++") 
#check("HOD/codes", "catalogueHOD", "C++")

check("lognormal/codes", "lognormal", "C++") 

check("numberCounts/codes", "numberCounts", "C++")
check("numberCounts/codes", "numberCounts_errors", "C++")

check("clustering/codes", "2pt_monopole", "C++") 
check("clustering/codes", "2pt_monopole_errors", "C++")
check("clustering/codes", "2pt_multipoles", "C++")
check("clustering/codes", "2pt_2D", "C++")
check("clustering/codes", "2pt_projected", "C++")
check("clustering/codes", "2pt_angular", "C++")
#check("clustering/codes", "power_spectrum_angular", "C++") 
check("clustering/codes", "3pt", "C++")
check("clustering/codes", "3pt_multipoles", "C++")

check("clustering/codes", "model_2pt_monopole_BAO", "C++")
check("clustering/codes", "model_2pt_monopole_RSD", "C++")
check("clustering/codes", "model_2pt_projected", "C++") 
check("clustering/codes", "model_2pt_2D", "C++") 
check("clustering/codes", "model_2pt_multipoles", "C++")
check("clustering/codes", "model_power_spectrum_multipoles", "C++")
#check("clustering/codes", "model_power_spectrum_angular", "C++") 
check("clustering/codes", "model_3pt", "C++") 

check("cosmicVoids/codes", "cleanVoidCatalogue", "C++")
check("cosmicVoids/codes", "modelling_VoidAbundances", "C++") 

check("parameterFile", "readParameterFile", "C++")

check("emulators", "xi_CAMB_emulator", "C++") 

check("maps", "maps", "C++") 

####################################################


### compile the CBL in python ###

if not ("nopy" in sys.argv):
    os.chdir(cwd)
    out = os.system("cd Examples/statistics/codes && make modelpy && cd ../../.. && make cleanpy")
    check_sys(out)
    if ("ALL" in sys.argv or "python" in sys.argv):
        os.system("make python PY=python3")
    elif ("python2" in sys.argv):
        os.system("make python PY=python2")   

        
####################################################


### check the python examples -> check(directory, example_executable, language, input) ###

check("data", "table.py", "python2")
check("data", "table.py", "python")

check("funcGrid", "funcgrid_bspline.py", "python2")
check("funcGrid", "funcgrid_bspline.py", "python")

check("wrappers", "fft_fftlog.py", "python2")
check("wrappers", "fft_fftlog.py", "python")

check("cosmology", "cosmology.py", "python2") 
check("cosmology", "cosmology.py", "python")

check("cosmology", "distances.py", "python2") 
check("cosmology", "distances.py", "python")

#check("cosmology", "velocities.py", "python2") 
#check("cosmology", "velocities.py", "python")

check("cosmology", "three_point_correlation.py", "python2") 
check("cosmology", "three_point_correlation.py", "python")

check("cosmology", "size_function.py", "python2") 
check("cosmology", "size_function.py", "python")

check("cosmology", "ng.py", "python2") 
check("cosmology", "ng.py", "python")

check("cosmology", "bias.py", "python2") 
check("cosmology", "bias.py", "python")

check("cosmology", "mass_function.py", "python2") 
check("cosmology", "mass_function.py", "python")

#check("cosmology", "mass_function_fR.py", "python2") 
#check("cosmology", "mass_function_fR.py", "python")

check("cosmology", "bao.py", "python2") 
check("cosmology", "bao.py", "python")

#check("cosmology", "fsigma8.py", "python2") 
#check("cosmology", "fsigma8.py", "python")

check("cosmology", "mass_growth.py", "python2") 
check("cosmology", "mass_growth.py", "python")

check("cosmology", "Pk_DeltaDelta.py", "python2") 
check("cosmology", "Pk_DeltaDelta.py", "python")

check("cosmology", "xi0_Kaiser.py", "python2") 
check("cosmology", "xi0_Kaiser.py", "python")

check("cosmology", "sigma.py", "python2") 
check("cosmology", "sigma.py", "python")

check("cosmology", "xi_matter.py", "python2") 
check("cosmology", "xi_matter.py", "python")

check("statistics/codes", "prior.py", "python2")
check("statistics/codes", "prior.py", "python")

#check("statistics/codes", "fit.py", "python2")
#check("statistics/codes", "fit.py", "python")

check("meshGrid/codes", "meshGridTest.py", "python2")
check("meshGrid/codes", "meshGridTest.py", "python")

#check("catalogue", "divide_catalogue.py", "python2")
#check("catalogue", "divide_catalogue.py", "python")

check("catalogue", "catalogue.py", "python2")
check("catalogue", "catalogue.py", "python")

check("catalogue", "mask_catalogue.py", "python2")
check("catalogue", "mask_catalogue.py", "python")

check("clustering/codes", "2pt_monopole.py", "python2")
check("clustering/codes", "2pt_monopole.py", "python")

check("clustering/codes", "2pt_model.py", "python2")
check("clustering/codes", "2pt_model.py", "python")

check("clustering/codes", "2pt_model_zErrors.py", "python2")
check("clustering/codes", "2pt_model_zErrors.py", "python")

check("clustering/codes", "3pt_model.py", "python2")
check("clustering/codes", "3pt_model.py", "python")

check("clustering/codes", "3pt.py", "python2")
check("clustering/codes", "3pt.py", "python")

check("cosmicVoids/codes", "cleanVoidCatalogue.py", "python2", "../input/parameter_file.ini") 
check("cosmicVoids/codes", "cleanVoidCatalogue.py", "python", "../input/parameter_file.ini")

check("parameterFile", "parameter_file.py", "python2") 
check("parameterFile", "parameter_file.py", "python") 

check("emulators", "xi_CAMB_emulator.py", "python2") 
check("emulators", "xi_CAMB_emulator.py", "python") 


####################################################


### create the documentation ###

if not ("nodoc" in sys.argv):
    os.chdir(cwd)
    out = os.system("make documentation")
    check_sys(out)
    print("".join(open("Doc/WARNING_LOGFILE")))
