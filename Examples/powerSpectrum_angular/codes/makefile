CXX = g++

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.4.0/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2.1/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ1 = power_spectrum_angular.o
OBJ2 = model_power_spectrum_angular.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

all: power_spectrum_angular model_power_spectrum_angular

power_spectrum_angular: $(OBJ1) 
	$(CXX) $(OBJ1) -o power_spectrum_angular $(FLAGS_LIB)
model_power_spectrum_angular: $(OBJ2) 
	$(CXX) $(OBJ2) -o model_power_spectrum_angular $(FLAGS_LIB)

clean:
	rm -f *.o power_spectrum_angular model_power_spectrum_angular *~ \#* temp* core*

power_spectrum_angular.o: power_spectrum_angular.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c power_spectrum_angular.cpp
model_power_spectrum_angular.o: model_power_spectrum_angular.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c model_power_spectrum_angular.cpp
