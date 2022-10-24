CXX = g++

FLAGS0 = -std=c++11 -fopenmp
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.4.0/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2.1/

FLAGS_LIB = -Wl,-rpath,$(dirLib) -L$(dirLib) -lWRAP_LIB
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ1 = integration_gsl.o
OBJ2 = minimisation_gsl.o
OBJ3 = integration_cuba.o
OBJ4 = fits.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

integration_gsl: $(OBJ1) 
	$(CXX) $(OBJ1) -o integration_gsl $(FLAGS_LIB)

minimisation_gsl: $(OBJ2) 
	$(CXX) $(OBJ2) -o minimisation_gsl $(FLAGS_LIB) 

integration_cuba: $(OBJ3) 
	$(CXX) $(OBJ3) -o integration_cuba $(FLAGS_LIB) 

fits: $(OBJ4) 
	$(CXX) $(OBJ4) -o fits $(FLAGS_LIB) 

clean:
	rm -f *.o integration_gsl minimisation_gsl integration_cuba fits catalogue_with_weights.fits *~ \#* temp* core*

integration_gsl.o: integration_gsl.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c integration_gsl.cpp

minimisation_gsl.o: minimisation_gsl.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c minimisation_gsl.cpp 

integration_cuba.o: integration_cuba.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c integration_cuba.cpp 

fits.o: fits.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c fits.cpp 

