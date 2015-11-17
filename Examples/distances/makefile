CC = g++

FLAGS = -O3 -std=c++11

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirNum = $(PWD)/../../../Numerical/
dirEH = $(dirLib)Cosmology/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lCOSM 
FLAGS_INC = -I$(dirH) -I$(dirHO) -I$(dirNum) -I$(dirEH) 

OBJ = distances.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

distances: $(OBJ) 
	$(CC) $(OBJ) -o distances $(FLAGS_LIB) 

clean:
	rm *.o distances *~ \#* temp* core* -f

distances.o: distances.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c distances.cpp 
