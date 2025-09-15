CXX = /opt/homebrew/bin/g++-15

FLAGS0 = -std=c++17 -fopenmp
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

dirLib = $(PWD)/../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.4.0/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2.1/

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL
FLAGS_INC = -I/opt/homebrew/include -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

OBJ0 = LCDM.o
OBJ1 = cosmology.o
OBJ2 = distances.o 
OBJ3 = model_cosmology.o
OBJ4 = growth_rate.o
OBJ5 = Pk_BoltzmannSolver.o
OBJ6 = velocities.o
OBJ7 = three_point_correlation.o
OBJ8 = size_function.o
OBJ9 = ng.o 
OBJ10 = bias.o
OBJ11 = mass_function.o
OBJ12 = fsigma8.o
OBJ13 = bao.o
OBJ14 = mass_growth.o
OBJ15 = Pk_DeltaDelta.o
OBJ16 = xi0_Kaiser.o
OBJ17 = sigma.o
OBJ18 = xi_matter.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif


LCDM: $(OBJ0) 
	$(CXX) $(OBJ0) -o LCDM $(FLAGS_LIB)

cosmology: $(OBJ1) 
	$(CXX) $(OBJ1) -o cosmology $(FLAGS_LIB)

distances: $(OBJ2) 
	$(CXX) $(OBJ2) -o distances $(FLAGS_LIB) 

model_cosmology: $(OBJ3) 
	$(CXX) $(OBJ3) -o model_cosmology $(FLAGS_LIB)

growth_rate: $(OBJ4) 
	$(CXX) $(OBJ4) -o growth_rate $(FLAGS_LIB)

Pk_BoltzmannSolver: $(OBJ5) 
	$(CXX) $(OBJ5) -o Pk_BoltzmannSolver $(FLAGS_LIB)

velocities: $(OBJ6) 
	$(CXX) $(OBJ6) -o velocities $(FLAGS_LIB)

three_point_correlation: $(OBJ7) 
	$(CXX) $(OBJ7) -o three_point_correlation $(FLAGS_LIB)

size_function: $(OBJ8) 
	$(CXX) $(OBJ8) -o size_function $(FLAGS_LIB)

ng: $(OBJ9) 
	$(CXX) $(OBJ9) -o ng $(FLAGS_LIB)

bias: $(OBJ10) 
	$(CXX) $(OBJ10) -o bias $(FLAGS_LIB)

mass_function: $(OBJ11) 
	$(CXX) $(OBJ11) -o mass_function $(FLAGS_LIB)

fsigma8: $(OBJ12) 
	$(CXX) $(OBJ12) -o fsigma8 $(FLAGS_LIB)

bao: $(OBJ13)
	$(CXX) $(OBJ13) -o bao $(FLAGS_LIB)

mass_growth: $(OBJ14)
	$(CXX) $(OBJ14) -o mass_growth $(FLAGS_LIB)

Pk_DeltaDelta: $(OBJ15)
	$(CXX) $(OBJ15) -o Pk_DeltaDelta $(FLAGS_LIB)

xi0_Kaiser: $(OBJ16)
	$(CXX) $(OBJ16) -o xi0_Kaiser $(FLAGS_LIB)

sigma: $(OBJ17)
	$(CXX) $(OBJ17) -o sigma $(FLAGS_LIB)

xi_matter: $(OBJ18)
	$(CXX) $(OBJ18) -o xi_matter $(FLAGS_LIB)


clean:
	rm -f *.o LCDM cosmology distances model_cosmology growth_rate Pk_BoltzmannSolver three_point_correlation velocities size_function ng bias mass_function fsigma8 bao mass_growth Pk_DeltaDelta xi0_Kaiser sigma xi_matter *~ \#* temp* core*


LCDM.o: LCDM.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c LCDM.cpp

cosmology.o: cosmology.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c cosmology.cpp

distances.o: distances.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c distances.cpp 

model_cosmology.o: model_cosmology.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c model_cosmology.cpp

growth_rate.o: growth_rate.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c growth_rate.cpp

Pk_BoltzmannSolver.o: Pk_BoltzmannSolver.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c Pk_BoltzmannSolver.cpp

velocities.o: velocities.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c velocities.cpp

three_point_correlation.o: three_point_correlation.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c three_point_correlation.cpp

size_function.o: size_function.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c size_function.cpp

ng.o: ng.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c ng.cpp

bias.o: bias.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c bias.cpp

mass_function.o: mass_function.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c mass_function.cpp

fsigma8.o: fsigma8.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c fsigma8.cpp

bao.o: bao.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c bao.cpp

mass_growth.o: mass_growth.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c mass_growth.cpp

Pk_DeltaDelta.o: Pk_DeltaDelta.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c Pk_DeltaDelta.cpp

xi0_Kaiser.o: xi0_Kaiser.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c xi0_Kaiser.cpp

sigma.o: sigma.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c sigma.cpp

xi_matter.o: xi_matter.cpp makefile $(dirLib)*.$(ES)
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c xi_matter.cpp
