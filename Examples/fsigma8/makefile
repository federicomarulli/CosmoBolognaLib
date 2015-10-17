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

OBJ = fsigma8.o 

fsigma8: $(OBJ) 
	$(CC) $(OBJ) -o fsigma8 $(FLAGS_LIB) 

clean:
	rm *.o fsigma8 *~ \#* temp* core* -f

fsigma8.o: fsigma8.cpp makefile $(dirLib)*.so
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fsigma8.cpp 
