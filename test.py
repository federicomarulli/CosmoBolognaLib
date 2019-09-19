#!/usr/bin/python

import os
import sys


####################################################

### check the output of os.system ###

def check_sys (out):
    if (out>0):
        sys.exit(0)

### check if valgrind got any errors ###

def check_valgrind ():
    if "ERROR SUMMARY: 0 errors from 0" in open("valgrind_output.dat").read():
        print "\n valgrind: OK! \n"
        os.remove("valgrind_output.dat")
    else:
        print "\n Error from vagrind! \n"
        print "".join(file("valgrind_output.dat"))
        os.remove("valgrind_output.dat")
        sys.exit(0)

        
### run the example code (with or withoud valgrind) ###

flags_Valgrind = "--track-origins=yes" #"--leak-check=full"

def com (ex, language, input):
    if not ("valgrind" in sys.argv):
        if (language=="C++"):
            out = os.system("./"+ex+"  "+input)
            check_sys(out)
        elif (language=="python"):
            out = os.system("python "+ex+"  "+input)
            check_sys(out)
        elif (language=="python3"):
            out = os.system("python3 "+ex+"  "+input)
            check_sys(out)
        else:
            print "\n Error: language not allowed! \n"
            sys.exit(0)
    else:
        if (language=="C++"):
            out = os.system("valgrind "+flags_Valgrind+" --log-file=\"valgrind_output.dat\" ./"+ex+"  "+input)
            check_sys(out)
            check_valgrind()
        elif (language=="python"):
            out = os.system("python "+ex+"  "+input)
            check_sys(out)
        elif (language=="python3"):
            out = os.system("python3 "+ex+"  "+input)
            check_sys(out)             
        else:
            print "\n Error: language not allowed! \n"
            sys.exit(0)

        
### check if the example is ok ###

def check (dir, file, language, input=""):
    if (("all" in sys.argv and language!="python3") or file in sys.argv or ("python" in sys.argv and language=="python") or ("python3" in sys.argv and language=="python3") or ("C++" in sys.argv and language=="C++")):
        print "\n\n-----> Testing ", file, "<-----\n"
        os.chdir(cwd+"/Examples/"+dir)
        com(file, language, input)

        
####################################################


os.system("clear")
cwd = os.getcwd()


### compile the CBL ###

if not ("nocompile" in sys.argv):
    out = os.system("make clean && make FLAGS=\"-O0 -g\" && make allExamples FLAGS=\"-O0 -g\"")
    check_sys(out)
    

####################################################


### check the C++ examples -> check(directory, example_executable, language, input) ###

check("vectors", "vectors", "C++") 

check("eigen", "eigen", "C++") 

check("randomNumbers", "randomNumbers", "C++") 
check("randomNumbers", "randomNumbers_custom", "C++") 
check("randomNumbers", "correlated_samples", "C++") 

check("wrappers", "integration_gsl", "C++")
check("wrappers", "minimisation_gsl", "C++") 
#check("wrappers", "integration_cuba", "C++") # seg fault with valgrind!!!
check("wrappers", "fits", "C++") 

check("covsample", "covsample", "C++") 

check("cosmology", "cosmology", "C++") 
check("cosmology", "fsigma8", "C++")
check("cosmology", "distances", "C++") 
check("cosmology", "model_cosmology", "C++") 

check("data", "data1D", "C++")

check("statistics/codes", "prior", "C++")
check("statistics/codes", "fit", "C++") 
check("statistics/codes", "sampler", "C++")

check("catalogue", "catalogue", "C++") 

check("numberCounts/codes", "numberCounts", "C++")
check("numberCounts/codes", "numberCounts_errors", "C++")

check("clustering/codes", "2pt_monopole", "C++") 
check("clustering/codes", "2pt_monopole_errors", "C++")
check("clustering/codes", "2pt_multipoles", "C++")
check("clustering/codes", "2pt_2D", "C++")
check("clustering/codes", "2pt_projected", "C++")
check("clustering/codes", "2pt_angular", "C++")
check("clustering/codes", "3pt", "C++")
check("clustering/codes", "3pt_multipoles", "C++")

check("clustering/codes", "model_2pt_monopole_BAO", "C++")
check("clustering/codes", "model_2pt_monopole_RSD", "C++")
check("clustering/codes", "model_2pt_projected", "C++") 
check("clustering/codes", "model_2pt_2D", "C++") 
check("clustering/codes", "model_2pt_multipoles", "C++") 
check("clustering/codes", "model_3pt", "C++") 

check("cosmicVoids/codes", "sizeFunction", "C++") 
check("cosmicVoids/codes", "cleanVoidCatalogue", "C++")
check("cosmicVoids/codes", "modelling_VoidAbundances", "C++") 

check("readParameterFile", "readParameterFile", "C++") 


####################################################


### compile the CBL in python ###

if not ("nopy" in sys.argv):
    os.chdir(cwd)
    out = os.system("make clean && make cleanpy && make FLAGS=\"-O0 -g\" && make CAMB && make python")
    check_sys(out)
    
    
####################################################


### check the python examples -> check(directory, example_executable, language, input) ###

check("wrappers", "fft_fftlog.py", "python")
check("wrappers", "fft_fftlog.py", "python3")

check("cosmology", "distances.py", "python") 
check("cosmology", "distances.py", "python3") 

check("statistics/codes", "prior.py", "python")
check("statistics/codes", "prior.py", "python3")

check("statistics/codes", "fit.py", "python")
check("statistics/codes", "fit.py", "python3")

check("catalogue", "catalogue.py", "python")
check("catalogue", "catalogue.py", "python3")

check("clustering/codes", "2pt_monopole.py", "python")
check("clustering/codes", "2pt_monopole.py", "python3")

check("clustering/codes", "2pt_model.py", "python")
check("clustering/codes", "2pt_model.py", "python3")

check("clustering/codes", "2pt_model_zErrors.py", "python")
check("clustering/codes", "2pt_model_zErrors.py", "python3")

check("clustering/codes", "3pt_model.py", "python")
check("clustering/codes", "3pt_model.py", "python3")

check("clustering/codes", "3pt.py", "python")
check("clustering/codes", "3pt.py", "python3")

check("cosmicVoids/codes", "sizeFunction.py", "python")
check("cosmicVoids/codes", "sizeFunction.py", "python3")

check("cosmicVoids/codes", "cleanVoidCatalogue.py", "python", "../input/parameter_file.ini") 
check("cosmicVoids/codes", "cleanVoidCatalogue.py", "python3", "../input/parameter_file.ini")


####################################################


### create the documentation ###

if not ("nodoc" in sys.argv):
    os.chdir(cwd)
    out = os.system("make documentation")
    check_sys(out)
    print "".join(file("Doc/WARNING_LOGFILE"))
