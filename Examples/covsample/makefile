CC = g++

FLAGS = -O3 -std=c++11

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO)
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
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c covsample.cpp 
