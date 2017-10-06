C = g++

FLAGS = -O3 -unroll -std=c++11 -Wextra #-openmp
#FLAGS = -O0 -g -std=c++11

dir = $(PWD)/
dirLib = $(dir)../../
dirH = $(dir)../../Headers/Lib/
dirHO = $(dir)../../Headers/Objects/	

FLAGS_LIB = -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC
FLAGS_INC = -I$(dir) -I$(dirH) -I$(dirHO)
FLAGS_GSL = -lgsl -lgslcblas -lm

HH = $(dirH)*.h

OBJ = $(dir)Prior.o


#################################################################### 


libPRIOR.so: $(OBJ) $(dir)makefile
	$(C) -shared -o $(dirLib)libPRIOR.so $(OBJ) $(FLAGS_LIB) $(FLAGS_GSL)

clean:
	rm -f *.o *~ temp* \#* core* 

#################################################################### 


$(dir)Prior.o: $(dir)Prior.cpp $(HH) $(dir)makefile
	$(C) $(FLAGS) -c -fPIC $(FLAGS_INC) $(dir)Prior.cpp -o $(dir)Prior.o
