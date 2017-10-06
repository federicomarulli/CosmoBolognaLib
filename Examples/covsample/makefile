CC = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirCUBA)
FLAGS_GSL = -lgsl -lgslcblas -lm 

OBJ = covsample.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

covsample: $(OBJ) 
	$(CC) $(OBJ) -o covsample $(FLAGS_LIB) $(FLAGS_GSL)

clean:
	rm -f *.o covsample *~ \#* temp* core*

covsample.o: covsample.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c covsample.cpp 
