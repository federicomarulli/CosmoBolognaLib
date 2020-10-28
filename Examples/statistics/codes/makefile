CXX = g++
SWIG = swig
PY = python

FLAGS0 = -std=c++11 -fopenmp 
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror
FLAGST = $(FLAGS0) $(FLAGS)

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.3.7/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM
FLAGS_INC = -I$(HOME)/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits) 

FLAGS_LINK = -shared
ES = so

# Python-related flags
python_version_full := $(wordlist 2,4,$(subst ., ,$(shell $(PY) --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})

ifeq ($(python_version_major),2)
	PYINC = $(shell $(PY) -c 'from distutils import sysconfig; print sysconfig.get_config_var("INCLUDEDIR")')
	SWIG_FLAG = -python -c++
	PYVERSION = $(python_version_major).$(python_version_minor)
endif
ifeq ($(python_version_major),3)
	PYINC = $(shell $(PY) -c 'from distutils import sysconfig; print(sysconfig.get_config_var("INCLUDEDIR"))')
	SWIG_FLAG = -python -c++ -py3
	PYVERSION = $(python_version_major).$(python_version_minor)
endif

PFLAGS = -I$(PYINC)/python$(PYVERSION)


Dvar = -DLINUX

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
	Dvar = -DMAC 
	FLAGS_LINK = -dynamiclib -undefined dynamic_lookup
        ES = dylib
	FLAGS_PY = -L$(PYLIB) -lpython$(PYVERSION) -ldl	
endif

OBJ1 = prior.o
OBJ2 = fit.o
OBJ3 = sampler.o

ES = so

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
        ES = dylib
endif

prior: $(OBJ1) 
	$(CXX) $(OBJ1) -o prior $(FLAGS_LIB) $(FLAGS_GSL)

fit: $(OBJ2) 
	$(CXX) $(OBJ2) -o fit $(FLAGS_LIB) $(FLAGS_GSL)

sampler: $(OBJ3) 
	$(CXX) $(OBJ3) -o sampler $(FLAGS_LIB) $(FLAGS_GSL)

modelpy: ModelPY_wrap.o ModelPY.i
	$(CXX) $(FLAGS_LINK) -o _modelpy.so ModelPY_wrap.o $(FLAGS_LIB)

clean:
	rm -f *.o *.cxx ModelPY_wrap.cxx modelpy.py _modelpy.so *~ \#* temp* core* prior fit sampler
	rm -rf __pycache__

prior.o: prior.cpp makefile $(dirLib)*.$(ES) 
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c prior.cpp

fit.o: fit.cpp makefile $(dirLib)*.$(ES) 
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c fit.cpp

sampler.o: sampler.cpp makefile $(dirLib)*.$(ES) 
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) $(varDIR) -c sampler.cpp

ModelPY_wrap.o: ModelPY_wrap.cxx ModelPY.i ModelPY.h
	$(CXX) $(FLAGST) -Wno-uninitialized -Wno-missing-field-initializers -Wno-unused-parameter $(PFLAGS) -c -fPIC $(FLAGS_INC) ModelPY_wrap.cxx -o ModelPY_wrap.o

ModelPY_wrap.cxx: ModelPY.h ModelPY.i
	$(call colorecho, "\n"Running swig... "\n")
	$(SWIG) $(SWIG_FLAG) -I$(dirH) -I$(dirLib)/Python/ ModelPY.i
