CC = g++

FLAGS = -O3 -std=c++11

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirEH) 

OBJ = randomNumbers.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

randomNumbers: $(OBJ) 
	$(CC) $(OBJ) -o randomNumbers $(FLAGS_LIB) 

clean:
	rm -f *.o randomNumbers *~ \#* temp* core*

randomNumbers.o: randomNumbers.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c randomNumbers.cpp 
