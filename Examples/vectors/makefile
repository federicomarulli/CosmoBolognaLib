CC = g++

FLAGS = -O3 -std=c++11

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirNum = $(PWD)/../../../Numerical/
dirEH = $(dirLib)Cosmology/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirHO) -I$(dirNum) -I$(dirEH) 

OBJ = vectors.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

vectors: $(OBJ) 
	$(CC) $(OBJ) -o vectors $(FLAGS_LIB) 

clean:
	rm *.o vectors *~ \#* temp* core* -f

vectors.o: vectors.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c vectors.cpp 
