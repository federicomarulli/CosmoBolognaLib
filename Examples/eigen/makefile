CXX = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -march=native -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/
dirCUBA = $(dirLib)External/Cuba-4.2/
dir_CCfits = $(dirLib)External/CCfits/include
dir_Eigen = $(dirLib)External/eigen-3.3.4/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNCGRID -lRAN -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ1 = eigen.o
OBJ2 = spherical_harmonics.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

eigen: $(OBJ1) 
	$(CXX) $(OBJ1) -o eigen $(FLAGS_LIB) 

spherical_harmonics: $(OBJ2) 
	$(CXX) $(OBJ2) -o spherical_harmonics $(FLAGS_LIB) 

clean:
	rm -f *.o eigen *~ \#* temp* core*

eigen.o: eigen.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c eigen.cpp 

spherical_harmonics.o: spherical_harmonics.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c spherical_harmonics.cpp 
