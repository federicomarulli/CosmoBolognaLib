CXX = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.3.7/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ1 = cosmology.o 
OBJ2 = fsigma8.o
OBJ3 = distances.o 
OBJ4 = model_cosmology.o
OBJ5 = Pk_dynamical_DE.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

cosmology: $(OBJ1) 
	$(CXX) $(OBJ1) -o cosmology $(FLAGS_LIB) 

fsigma8: $(OBJ2) 
	$(CXX) $(OBJ2) -o fsigma8 $(FLAGS_LIB) 

distances: $(OBJ3) 
	$(CXX) $(OBJ3) -o distances $(FLAGS_LIB) 

model_cosmology: $(OBJ4) 
	$(CXX) $(OBJ4) -o model_cosmology $(FLAGS_LIB)

Pk_dynamical_DE: $(OBJ5) 
	$(CXX) $(OBJ5) -o Pk_dynamical_DE $(FLAGS_LIB)

clean:
	rm -f *.o cosmology fsigma8 distances model_cosmology Pk_dynamical_DE *~ \#* temp* core*

cosmology.o: cosmology.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c cosmology.cpp

fsigma8.o: fsigma8.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fsigma8.cpp 

distances.o: distances.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c distances.cpp 

model_cosmology.o: model_cosmology.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_cosmology.cpp

Pk_dynamical_DE.o: Pk_dynamical_DE.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c Pk_dynamical_DE.cpp
