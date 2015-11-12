CC = g++

FLAGS = -O3 -std=c++11

DIR = $(PWD)/../../../

dirLib = $(DIR)CosmoBolognaLib/
varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirNum = $(DIR)Numerical/
dirEH = $(dirLib)Cosmology/EH/

FLAGS_LIB = -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lCOSM 
FLAGS_INC = -I$(dirH) -I$(dirHO) -I$(dirNum) -I$(dirEH) 

OBJ = distances.o 

distances: $(OBJ) 
	$(CC) $(OBJ) -o distances $(FLAGS_LIB) 

clean:
	rm *.o distances *~ \#* temp* core* -f

distances.o: distances.cpp makefile $(dirLib)*.so
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c distances.cpp 
