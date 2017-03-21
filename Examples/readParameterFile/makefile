CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lREADP
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirEH)

OBJ = readParameterFile.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

readParameterFile: $(OBJ) 
	$(CC) $(OBJ) -o readParameterFile $(FLAGS_LIB)

clean:
	rm -f *.o readParameterFile *~ \#* temp* core*


readParameterFile.o: readParameterFile.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c readParameterFile.cpp

