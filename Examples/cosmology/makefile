C = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/eigen-3.3.4/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ1 = cosmology.o 
OBJ2 = fsigma8.o
OBJ3 = distances.o 
OBJ4 = model_cosmology.o 

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

cosmology: $(OBJ1) 
	$(C) $(OBJ1) -o cosmology $(FLAGS_LIB) 

fsigma8: $(OBJ2) 
	$(C) $(OBJ2) -o fsigma8 $(FLAGS_LIB) 

distances: $(OBJ3) 
	$(C) $(OBJ3) -o distances $(FLAGS_LIB) 

model_cosmology: $(OBJ4) 
	$(C) $(OBJ4) -o model_cosmology $(FLAGS_LIB) 

clean:
	rm -f *.o cosmology fsigma8 distances model_cosmology *~ \#* temp* core*

cosmology.o: cosmology.cpp makefile $(dirLib)*.$(ES)
	$(C) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c cosmology.cpp

fsigma8.o: fsigma8.cpp makefile $(dirLib)*.$(ES)
	$(C) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fsigma8.cpp 

distances.o: distances.cpp makefile $(dirLib)*.$(ES)
	$(C) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c distances.cpp 

model_cosmology.o: model_cosmology.cpp makefile $(dirLib)*.$(ES)
	$(C) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_cosmology.cpp
