CC = g++

FLAGS = -O3 -std=c++11

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lCOSM 
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirEH) -isystem $(HOME)/Numerical/

OBJ = distances.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

distances: $(OBJ) 
	$(CC) $(OBJ) -o distances $(FLAGS_LIB) 

clean:
	rm -f *.o distances *~ \#* temp* core*

distances.o: distances.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c distances.cpp 
