CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM -lCM -lCAT
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirEH)

OBJ = catalogue.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

catalogue: $(OBJ) 
	$(CC) $(OBJ) -o catalogue $(FLAGS_LIB)

clean:
	rm -f *.o catalogue *~ \#* temp* core*


catalogue.o: catalogue.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c catalogue.cpp

