CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT 
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirEH) -isystem $(HOME)/Numerical/

OBJ1 = prior.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

prior: $(OBJ1) 
	$(CC) $(OBJ1) -o prior $(FLAGS_LIB) 

clean:
	rm -f *.o *~ \#* temp* core* prior

prior.o: prior.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c prior.cpp
