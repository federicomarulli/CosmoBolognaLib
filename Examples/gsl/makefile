CC = g++

FLAGS = -O3 -std=c++11

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirEH) 

OBJ1 = integration.o
OBJ2 = minimisation.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

integration: $(OBJ1) 
	$(CC) $(OBJ1) -o integration $(FLAGS_LIB) 

minimisation: $(OBJ2) 
	$(CC) $(OBJ2) -o minimisation $(FLAGS_LIB) 

clean:
	rm -f *.o integration minimisation *~ \#* temp* core*

integration.o: integration.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c integration.cpp 

minimisation.o: minimisation.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c minimisation.cpp 
