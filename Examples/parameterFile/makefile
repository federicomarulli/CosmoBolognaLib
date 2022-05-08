CXX = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror 

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.4.0/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2.1/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lPARF
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ = readParameterFile.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

readParameterFile: $(OBJ) 
	$(CXX) $(OBJ) -o readParameterFile $(FLAGS_LIB)

clean:
	rm -f *.o readParameterFile *~ \#* temp* core* params_from_scratch.ini new_params.ini


readParameterFile.o: readParameterFile.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c readParameterFile.cpp

