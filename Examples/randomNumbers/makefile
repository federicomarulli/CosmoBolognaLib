CXX = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.3.7/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lKERNEL -lRAN -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits) 

OBJ1 = randomNumbers.o
OBJ2 = randomNumbers_custom.o
OBJ3 = correlated_samples.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

randomNumbers: $(OBJ1) 
	$(CXX) $(OBJ1) -o randomNumbers $(FLAGS_LIB) 

randomNumbers_custom: $(OBJ2) 
	$(CXX) $(OBJ2) -o randomNumbers_custom $(FLAGS_LIB) 

correlated_samples: $(OBJ3) 
	$(CXX) $(OBJ3) -o correlated_samples $(FLAGS_LIB) 

clean:
	rm -f *.o randomNumbers randomNumbers_custom correlated_samples *~ \#* temp* core*

randomNumbers.o: randomNumbers.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c randomNumbers.cpp 

randomNumbers_custom.o: randomNumbers_custom.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c randomNumbers_custom.cpp 

correlated_samples.o: correlated_samples.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c correlated_samples.cpp 
