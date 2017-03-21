CC = g++

FLAGS = -O3 -std=c++11

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lCOSM
FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dirH) -I$(dirHO) -I$(dirEH)

OBJ = fsigma8.o 

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

fsigma8: $(OBJ) 
	$(CC) $(OBJ) -o fsigma8 $(FLAGS_LIB) 

clean:
	rm -f *.o fsigma8 *~ \#* temp* core*

fsigma8.o: fsigma8.cpp makefile $(dirLib)*.$(ES)
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fsigma8.cpp 
