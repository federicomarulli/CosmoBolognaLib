CXX = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/
dirCUBA = $(dirLib)External/Cuba-4.2/
dir_CCfits = $(dirLib)External/CCfits/include
dir_Eigen = $(dirLib)External/eigen-3.3.4/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNCGRID -lRAN 
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ = eigen.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

eigen: $(OBJ) 
	$(CXX) $(OBJ) -o eigen $(FLAGS_LIB) 

clean:
	rm -f *.o eigen *~ \#* temp* core*

eigen.o: eigen.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c eigen.cpp 
