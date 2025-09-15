CXX = g++

FLAGS0 = -std=c++17 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.4.0/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2.1/

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ1 = 2pt_monopole.o
OBJ2 = 2pt_monopole_errors.o
OBJ3 = 2pt_2D.o
OBJ4 = 2pt_projected.o
OBJ5 = 2pt_angular.o
OBJ6 = 2pt_multipoles.o
OBJ7 = power_spectrum_angular.o
OBJ8 = 3pt.o
OBJ9 = 3pt_multipoles.o
OBJ10 = model_2pt_monopole_BAO.o
OBJ11 = model_2pt_monopole_RSD.o
OBJ12 = model_2pt_projected.o
OBJ13 = model_2pt_2D.o
OBJ14 = model_2pt_multipoles.o
OBJ15 = model_power_spectrum_angular.o
OBJ16 = model_3pt.o
OBJ17 = model_power_spectrum_multipoles.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

all: 2pt_monopole 2pt_monopole_errors 2pt_2D 2pt_projected 2pt_angular 2pt_multipoles power_spectrum_angular 3pt 3pt_multipoles model_2pt_monopole_BAO model_2pt_monopole_RSD model_2pt_projected model_2pt_2D model_2pt_multipoles model_power_spectrum_angular model_3pt model_power_spectrum_multipoles

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

power_spectrum_angular: $(OBJ7) 
	$(CXX) $(OBJ7) -o power_spectrum_angular $(FLAGS_LIB)

3pt: $(OBJ8) 
	$(CXX) $(OBJ8) -o 3pt $(FLAGS_LIB)

3pt_multipoles: $(OBJ9) 
	$(CXX) $(OBJ9) -o 3pt_multipoles $(FLAGS_LIB)

model_2pt_monopole_BAO: $(OBJ10) 
	$(CXX) $(OBJ10) -o model_2pt_monopole_BAO $(FLAGS_LIB)

model_2pt_monopole_RSD: $(OBJ11) 
	$(CXX) $(OBJ11) -o model_2pt_monopole_RSD $(FLAGS_LIB)

model_2pt_projected: $(OBJ12) 
	$(CXX) $(OBJ12) -o model_2pt_projected $(FLAGS_LIB)

model_2pt_2D: $(OBJ13) 
	$(CXX) $(OBJ13) -o model_2pt_2D $(FLAGS_LIB)

model_2pt_multipoles: $(OBJ14) 
	$(CXX) $(OBJ14) -o model_2pt_multipoles $(FLAGS_LIB)

model_power_spectrum_angular: $(OBJ15) 
	$(CXX) $(OBJ15) -o model_power_spectrum_angular $(FLAGS_LIB)

model_3pt: $(OBJ16) 
	$(CXX) $(OBJ16) -o model_3pt $(FLAGS_LIB)

model_power_spectrum_multipoles: $(OBJ17) 
	$(CXX) $(OBJ17) -o model_power_spectrum_multipoles $(FLAGS_LIB)

clean:
	rm -f *.o 2pt_monopole 2pt_monopole_errors 2pt_2D 2pt_projected 2pt_angular 2pt_multipoles power_spectrum_angular 3pt 3pt_multipoles model_2pt_monopole_BAO model_2pt_monopole_RSD model_2pt_projected model_2pt_2D model_2pt_multipoles model_power_spectrum_angular model_3pt model_power_spectrum_multipoles *~ \#* temp* core* valgrind_output.dat


2pt_monopole.o: 2pt_monopole.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c 2pt_monopole.cpp

2pt_monopole_errors.o: 2pt_monopole_errors.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c 2pt_monopole_errors.cpp

2pt_2D.o: 2pt_2D.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c 2pt_2D.cpp

2pt_projected.o: 2pt_projected.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c 2pt_projected.cpp

2pt_angular.o: 2pt_angular.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c 2pt_angular.cpp

2pt_multipoles.o: 2pt_multipoles.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c 2pt_multipoles.cpp

power_spectrum_angular.o: power_spectrum_angular.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c power_spectrum_angular.cpp

3pt.o: 3pt.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c 3pt.cpp 

3pt_multipoles.o: 3pt_multipoles.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c 3pt_multipoles.cpp 

model_2pt_monopole_BAO.o: model_2pt_monopole_BAO.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_2pt_monopole_BAO.cpp

model_2pt_monopole_RSD.o: model_2pt_monopole_RSD.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_2pt_monopole_RSD.cpp

model_2pt_projected.o: model_2pt_projected.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_2pt_projected.cpp

model_2pt_2D.o: model_2pt_2D.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_2pt_2D.cpp

model_2pt_multipoles.o: model_2pt_multipoles.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_2pt_multipoles.cpp\

model_power_spectrum_angular.o: model_power_spectrum_angular.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_power_spectrum_angular.cpp

model_3pt.o: model_3pt.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_3pt.cpp 
	
model_power_spectrum_multipoles.o: model_power_spectrum_multipoles.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_power_spectrum_multipoles.cpp 
