CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirNum = $(PWD)/../../../../Numerical/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirNum) -I$(dirEH) 

OBJ1 = prior.o
OBJ2 = fit.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

prior: $(OBJ1) 
	$(CC) $(OBJ1) -o prior $(FLAGS_LIB) $(FLAGS_GSL)

fit: $(OBJ2) 
	$(CC) $(OBJ2) -o fit $(FLAGS_LIB) $(FLAGS_GSL)

clean:
	rm *.o *~ \#* temp* core* prior fit -f

prior.o: prior.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c prior.cpp

fit.o: fit.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fit.cpp
