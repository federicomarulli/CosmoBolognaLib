CC = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirEH = $(dirLib)External/EH/
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM -lCM -lCAT
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirEH) -I$(dirCUBA)

OBJ = catalogue.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)

        ES = dylib
endif

catalogue: $(OBJ) 
	$(CC) $(OBJ) -o catalogue $(FLAGS_LIB)

clean:
	rm -f *.o catalogue *~ \#* temp* core*


catalogue.o: catalogue.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c catalogue.cpp

