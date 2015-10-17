CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

DIR = $(PWD)/../../../../

dirLib = $(DIR)CosmoBolognaLib/
varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirNum = $(DIR)Numerical/
dirEH = $(dirLib)Cosmology/EH/

FLAGS_LIB = -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lRANDOM -lTWOP -lTHREEP -lGLOB
FLAGS_INC = -I$(dirH) -I$(dirHO) -I$(dirNum) -I$(dirEH) 

OBJ1 = 2pt.o
OBJ2 = 2pt_jackknife.o
OBJ3 = 3pt.o 

2pt: $(OBJ1) 
	$(CC) $(OBJ1) -o 2pt $(FLAGS_LIB)
2pt_jackknife: $(OBJ2) 
	$(CC) $(OBJ2) -o 2pt_jackknife $(FLAGS_LIB) 
3pt: $(OBJ3) 
	$(CC) $(OBJ3) -o 3pt $(FLAGS_LIB) 

clean:
	rm *.o 2pt 2pt_jackknife 3pt *~ \#* temp* core* -f

2pt.o: 2pt.cpp makefile $(dirLib)*.so
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt.cpp

2pt_jackknife.o: 2pt_jackknife.cpp makefile $(dirLib)*.so
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 2pt_jackknife.cpp

3pt.o: 3pt.cpp makefile $(dirLib)*.so
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c 3pt.cpp 
