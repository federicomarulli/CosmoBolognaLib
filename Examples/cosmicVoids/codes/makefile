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

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM -lCM -lCAT
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirEH) -I$(dirCUBA) 

OBJ1 = sizeFunction.o
OBJ2 = cleanVoidCatalogue.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

sizeFunction: $(OBJ1) 
	$(CC) $(OBJ1) -o sizeFunction $(FLAGS_LIB)

cleanVoidCatalogue: $(OBJ2) 
	$(CC) $(OBJ2) -o cleanVoidCatalogue $(FLAGS_LIB)

clean:
	rm -f *.o sizeFunction cleanVoidCatalogue *~ \#* temp* core*


sizeFunction.o: sizeFunction.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c sizeFunction.cpp

cleanVoidCatalogue.o: cleanVoidCatalogue.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c cleanVoidCatalogue.cpp

