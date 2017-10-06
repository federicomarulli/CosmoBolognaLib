CC = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirEH = $(dirLib)External/EH/
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirEH) -I$(dirCUBA) 

OBJ1 = integration_gsl.o
OBJ2 = integration_cuba.o
OBJ3 = minimisation.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

integration_gsl: $(OBJ1) 
	$(CC) $(OBJ1) -o integration_gsl $(FLAGS_LIB) 

integration_cuba: $(OBJ2) 
	$(CC) $(OBJ2) -o integration_cuba $(FLAGS_LIB) 

minimisation: $(OBJ3) 
	$(CC) $(OBJ3) -o minimisation $(FLAGS_LIB) 

clean:
	rm -f *.o integration_gsl integration_cuba minimisation *~ \#* temp* core*

integration_gsl.o: integration_gsl.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c integration_gsl.cpp 

integration_cuba.o: integration_cuba.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c integration_cuba.cpp 

minimisation.o: minimisation.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c minimisation.cpp 
