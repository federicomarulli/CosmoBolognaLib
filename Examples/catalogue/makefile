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

OBJ = catalogue.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)

        ES = dylib
endif

catalogue: $(OBJ) 
	$(CXX) $(OBJ) -o catalogue $(FLAGS_LIB)

clean:
	rm -f *.o catalogue *~ \#* temp* core*


catalogue.o: catalogue.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c catalogue.cpp

