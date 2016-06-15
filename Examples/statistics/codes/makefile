CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirHM = $(dirLib)Headers/Models/
dirNum = $(PWD)/../../../../Numerical/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT 
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirHO) -I$(dirHM) -I$(dirNum) -I$(dirEH) 

OBJ1 = prior.o
OBJ2 = fit_linear_relation.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

prior: $(OBJ1) 
	$(CC) $(OBJ1) -o prior $(FLAGS_LIB) $(FLAGS_GSL)

fit_linear_relation: $(OBJ2) 
	$(CC) $(OBJ2) -o fit_linear_relation $(FLAGS_LIB) $(FLAGS_GSL)

clean:
	rm *.o *~ \#* temp* core* prior fit_linear_relation -f

prior.o: prior.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c prior.cpp

fit_linear_relation.o: fit_linear_relation.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fit_linear_relation.cpp
