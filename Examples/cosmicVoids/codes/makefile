CXX = g++

FLAGS0 = -std=c++17 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.4.0/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2.1/

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL
FLAGS_INC = -I$(HOME)/include/  -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits) 

OBJ1 = cleanVoidCatalogue.o
OBJ2 = modelling_VoidAbundances.o
OBJ3 = BitVF_LaZeVo.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif


cleanVoidCatalogue: $(OBJ1) 
	$(CXX) $(OBJ1) -o cleanVoidCatalogue $(FLAGS_LIB)

modelling_VoidAbundances: $(OBJ2) 
	$(CXX) $(OBJ2) -o modelling_VoidAbundances $(FLAGS_LIB)

BitVF_LaZeVo: $(OBJ3) 
	$(CXX) $(OBJ3) -o BitVF_LaZeVo $(FLAGS_LIB)


clean:
	rm -f *.o cleanVoidCatalogue modelling_VoidAbundances BitVF_LaZeVo *~ \#* temp* core*


cleanVoidCatalogue.o: cleanVoidCatalogue.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c cleanVoidCatalogue.cpp

modelling_VoidAbundances.o: modelling_VoidAbundances.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c modelling_VoidAbundances.cpp

BitVF_LaZeVo.o: BitVF_LaZeVo.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c BitVF_LaZeVo.cpp

