CXX = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/eigen-3.3.4/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ1 = 2pt_monopole.o
OBJ2 = 2pt_monopole_errors.o
OBJ3 = 2pt_2D.o
OBJ4 = 2pt_projected.o
OBJ5 = 2pt_angular.o
OBJ6 = 2pt_multipoles.o
OBJ7 = 3pt.o
OBJ8 = 3pt_multipoles.o
OBJ9 = model_2pt_monopole_BAO.o
OBJ10 = model_2pt_monopole_RSD.o
OBJ11 = model_2pt_projected.o
OBJ12 = model_2pt_2D.o
OBJ13 = model_2pt_multipoles.o
OBJ14 = model_3pt.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

all: 2pt_monopole 2pt_monopole_errors 2pt_2D 2pt_projected 2pt_angular 2pt_multipoles 3pt 3pt_multipoles model_2pt_monopole_BAO model_2pt_monopole_RSD model_2pt_projected model_2pt_2D model_2pt_multipoles model_3pt

2pt_monopole: $(OBJ1) 
	$(CXX) $(OBJ1) -o 2pt_monopole $(FLAGS_LIB)

2pt_monopole_errors: $(OBJ2) 
	$(CXX) $(OBJ2) -o 2pt_monopole_errors $(FLAGS_LIB) 

2pt_2D: $(OBJ3) 
	$(CXX) $(OBJ3) -o 2pt_2D $(FLAGS_LIB)

2pt_projected: $(OBJ4) 
	$(CXX) $(OBJ4) -o 2pt_projected $(FLAGS_LIB) 

2pt_angular: $(OBJ5) 
	$(CXX) $(OBJ5) -o 2pt_angular $(FLAGS_LIB) 

2pt_multipoles: $(OBJ6) 
	$(CXX) $(OBJ6) -o 2pt_multipoles $(FLAGS_LIB) 

3pt: $(OBJ7) 
	$(CXX) $(OBJ7) -o 3pt $(FLAGS_LIB)

3pt_multipoles: $(OBJ8) 
	$(CXX) $(OBJ8) -o 3pt_multipoles $(FLAGS_LIB)

model_2pt_monopole_BAO: $(OBJ9) 
	$(CXX) $(OBJ9) -o model_2pt_monopole_BAO $(FLAGS_LIB)

model_2pt_monopole_RSD: $(OBJ10) 
	$(CXX) $(OBJ10) -o model_2pt_monopole_RSD $(FLAGS_LIB)

model_2pt_projected: $(OBJ11) 
	$(CXX) $(OBJ11) -o model_2pt_projected $(FLAGS_LIB)

model_2pt_2D: $(OBJ12) 
	$(CXX) $(OBJ12) -o model_2pt_2D $(FLAGS_LIB)

model_2pt_multipoles: $(OBJ13) 
	$(CXX) $(OBJ13) -o model_2pt_multipoles $(FLAGS_LIB)

model_3pt: $(OBJ14) 
	$(CXX) $(OBJ14) -o model_3pt $(FLAGS_LIB)

clean:
	rm -f *.o 2pt_monopole 2pt_monopole_errors 2pt_2D 2pt_projected 2pt_angular 2pt_multipoles 3pt 3pt_multipoles model_2pt_monopole_BAO model_2pt_monopole_RSD model_2pt_projected model_2pt_2D model_2pt_multipoles model_3pt *~ \#* temp* core*


2pt_monopole.o: 2pt_monopole.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_monopole.cpp

2pt_monopole_errors.o: 2pt_monopole_errors.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_monopole_errors.cpp

2pt_2D.o: 2pt_2D.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_2D.cpp

2pt_projected.o: 2pt_projected.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_projected.cpp

2pt_angular.o: 2pt_angular.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_angular.cpp

2pt_multipoles.o: 2pt_multipoles.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_multipoles.cpp

3pt.o: 3pt.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 3pt.cpp 

3pt_multipoles.o: 3pt_multipoles.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 3pt_multipoles.cpp 

model_2pt_monopole_BAO.o: model_2pt_monopole_BAO.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_monopole_BAO.cpp

model_2pt_monopole_RSD.o: model_2pt_monopole_RSD.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_monopole_RSD.cpp

model_2pt_projected.o: model_2pt_projected.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_projected.cpp

model_2pt_2D.o: model_2pt_2D.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_2D.cpp

model_2pt_multipoles.o: model_2pt_multipoles.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_2pt_multipoles.cpp

model_3pt.o: model_3pt.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_3pt.cpp 
