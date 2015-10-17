CC = g++

FLAGS = -O3 -std=c++11

DIR = $(PWD)/../../../

dirLib = $(DIR)CosmoBolognaLib/
varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirNum = $(DIR)Numerical/
dirEH = $(dirLib)Cosmology/EH/

FLAGS_LIB = -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(dirH) -I$(dirHO) -I$(dirNum) -I$(dirEH) 

OBJ = vectors.o 

vectors: $(OBJ) 
	$(CC) $(OBJ) -o vectors $(FLAGS_LIB) 

clean:
	rm *.o vectors *~ \#* temp* core* -f

vectors.o: vectors.cpp makefile $(dirLib)*.so
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c vectors.cpp 
