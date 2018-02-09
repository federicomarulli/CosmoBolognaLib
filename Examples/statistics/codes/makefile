C = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirEH = $(dirLib)External/EH/
dir_Eigen = $(dirLib)External/eigen-3.3.4/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirEH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits) 

OBJ1 = prior.o
OBJ2 = fit.o
OBJ3 = sampler.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

prior: $(OBJ1) 
	$(C) $(OBJ1) -o prior $(FLAGS_LIB) $(FLAGS_GSL)

fit: $(OBJ2) 
	$(C) $(OBJ2) -o fit $(FLAGS_LIB) $(FLAGS_GSL)

sampler: $(OBJ3) 
	$(C) $(OBJ3) -o sampler $(FLAGS_LIB) $(FLAGS_GSL)


clean:
	rm -f *.o *~ \#* temp* core* prior fit sampler

prior.o: prior.cpp makefile $(dirLib)*.$(ES) 
	$(C) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c prior.cpp

fit.o: fit.cpp makefile $(dirLib)*.$(ES) 
	$(C) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fit.cpp

sampler.o: sampler.cpp makefile $(dirLib)*.$(ES) 
	$(C) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c sampler.cpp
