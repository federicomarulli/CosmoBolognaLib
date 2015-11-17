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

OBJ = fsigma8.o 

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

fsigma8: $(OBJ) 
	$(CC) $(OBJ) -o fsigma8 $(FLAGS_LIB) 

clean:
	rm *.o fsigma8 *~ \#* temp* core* -f

fsigma8.o: fsigma8.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fsigma8.cpp 
