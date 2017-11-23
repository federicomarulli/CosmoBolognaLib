CC = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirEH = $(dirLib)External/EH/
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirEH) -I$(dirCUBA) 

OBJ1 = prior.o
OBJ2 = fit.o
OBJ3 = sampler.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

prior: $(OBJ1) 
	$(CC) $(OBJ1) -o prior $(FLAGS_LIB) $(FLAGS_GSL)

fit: $(OBJ2) 
	$(CC) $(OBJ2) -o fit $(FLAGS_LIB) $(FLAGS_GSL)

sampler: $(OBJ3) 
	$(CC) $(OBJ3) -o sampler $(FLAGS_LIB) $(FLAGS_GSL)


clean:
	rm -f *.o *~ \#* temp* core* prior fit sampler

prior.o: prior.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c prior.cpp

fit.o: fit.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fit.cpp

sampler.o: sampler.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c sampler.cpp
