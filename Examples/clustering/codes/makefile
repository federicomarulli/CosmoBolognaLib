CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP -lMODEL -lGLOB
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirEH) -isystem $(HOME)/Numerical/ 

OBJ1 = 2pt_monopole.o
OBJ2 = 2pt_monopole_errors.o
OBJ3 = 2pt_2D.o
OBJ4 = 2pt_projected.o
OBJ5 = 2pt_angular.o
OBJ6 = 3pt.o
OBJ7 = modelBias_2pt_projected.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

2pt_monopole: $(OBJ1) 
	$(CC) $(OBJ1) -o 2pt_monopole $(FLAGS_LIB)

2pt_monopole_errors: $(OBJ2) 
	$(CC) $(OBJ2) -o 2pt_monopole_errors $(FLAGS_LIB) 

2pt_2D: $(OBJ3) 
	$(CC) $(OBJ3) -o 2pt_2D $(FLAGS_LIB)

2pt_projected: $(OBJ4) 
	$(CC) $(OBJ4) -o 2pt_projected $(FLAGS_LIB) 

2pt_angular: $(OBJ5) 
	$(CC) $(OBJ5) -o 2pt_angular $(FLAGS_LIB) 

3pt: $(OBJ6) 
	$(CC) $(OBJ6) -o 3pt $(FLAGS_LIB)

modelBias_2pt_projected: $(OBJ7) 
	$(CC) $(OBJ7) -o modelBias_2pt_projected $(FLAGS_LIB)

clean:
	rm -f *.o 2pt_monopole 2pt_monopole_errors 2pt_2D 2pt_projected 2pt_angular 3pt modelBias_2pt_projected *~ \#* temp* core*


2pt_monopole.o: 2pt_monopole.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_monopole.cpp

2pt_monopole_errors.o: 2pt_monopole_errors.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_monopole_errors.cpp

2pt_2D.o: 2pt_2D.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_2D.cpp

2pt_projected.o: 2pt_projected.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_projected.cpp

2pt_angular.o: 2pt_angular.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_angular.cpp

3pt.o: 3pt.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 3pt.cpp 

modelBias_2pt_projected.o: modelBias_2pt_projected.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c modelBias_2pt_projected.cpp
