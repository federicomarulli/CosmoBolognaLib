CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirNum = $(PWD)/../../../../Numerical/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP -lMODEL -lGLOB
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirNum) -I$(dirEH) 

OBJ1 = 2pt_monopole.o
OBJ2 = 2pt_2D.o
OBJ3 = 2pt_jackknife.o
OBJ4 = 3pt.o 
OBJ5 = 2pt_projected_jackknife.o
OBJ6 = modelBias_2pt_projected.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

2pt_monopole: $(OBJ1) 
	$(CC) $(OBJ1) -o 2pt_monopole $(FLAGS_LIB)

2pt_2D: $(OBJ2) 
	$(CC) $(OBJ2) -o 2pt_2D $(FLAGS_LIB)

2pt_jackknife: $(OBJ3) 
	$(CC) $(OBJ3) -o 2pt_jackknife $(FLAGS_LIB) 

3pt: $(OBJ4) 
	$(CC) $(OBJ4) -o 3pt $(FLAGS_LIB) 

2pt_projected_jackknife: $(OBJ5) 
	$(CC) $(OBJ5) -o 2pt_projected_jackknife $(FLAGS_LIB) 

modelBias_2pt_projected: $(OBJ6) 
	$(CC) $(OBJ6) -o modelBias_2pt_projected $(FLAGS_LIB) 

clean:
	rm *.o 2pt_monopole 2pt_2D 2pt_jackknife 3pt 2pt_projected_jackknife modelBias_2pt_projected *~ \#* temp* core* -f


2pt_monopole.o: 2pt_monopole.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_monopole.cpp

2pt_2D.o: 2pt_2D.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_2D.cpp

2pt_jackknife.o: 2pt_jackknife.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_jackknife.cpp

3pt.o: 3pt.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 3pt.cpp 

2pt_projected_jackknife.o: 2pt_projected_jackknife.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_projected_jackknife.cpp

modelBias_2pt_projected.o: modelBias_2pt_projected.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c modelBias_2pt_projected.cpp

