CC = g++

FLAGS = -O3 -std=c++11

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO)

OBJ = vectors.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

vectors: $(OBJ) 
	$(CC) $(OBJ) -o vectors $(FLAGS_LIB) 

clean:
	rm -f *.o vectors *~ \#* temp* core*

vectors.o: vectors.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c vectors.cpp 
