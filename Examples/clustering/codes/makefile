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
OBJ7 = model_2pt_monopole_BAO.o
OBJ8 = model_2pt_monopole_RSD.o
OBJ9 = model_2pt_projected.o
OBJ10 = model_2pt_2D.o

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

model_2pt_monopole_BAO: $(OBJ7) 
	$(CC) $(OBJ7) -o model_2pt_monopole_BAO $(FLAGS_LIB)

model_2pt_monopole_RSD: $(OBJ8) 
	$(CC) $(OBJ8) -o model_2pt_monopole_RSD $(FLAGS_LIB)

model_2pt_projected: $(OBJ9) 
	$(CC) $(OBJ9) -o model_2pt_projected $(FLAGS_LIB)

model_2pt_2D: $(OBJ10) 
	$(CC) $(OBJ10) -o model_2pt_2D $(FLAGS_LIB)

clean:
	rm -f *.o 2pt_monopole 2pt_monopole_errors 2pt_2D 2pt_projected 2pt_angular 3pt model_2pt_monopole_BAO model_2pt_monopole_RSD model_2pt_projected model_2pt_2D *~ \#* temp* core*


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

model_2pt_monopole_BAO.o: model_2pt_monopole_BAO.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_monopole_BAO.cpp

model_2pt_monopole_RSD.o: model_2pt_monopole_RSD.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_monopole_RSD.cpp

model_2pt_projected.o: model_2pt_projected.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_projected.cpp

model_2pt_2D.o: model_2pt_2D.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_2D.cpp
