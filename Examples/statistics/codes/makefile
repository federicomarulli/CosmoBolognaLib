CXX = /opt/homebrew/opt/llvm/bin/clang++
SWIG = swig
PY = python

# Set the flags. Note that '-fopenmp' is now added with the required library paths
FLAGS0 = -std=c++17 -fopenmp=libomp 
FLAGS = -O3 -Wall -Wextra -pedantic -Wfatal-errors
FLAGST = $(FLAGS0) $(FLAGS)

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/
dir_Eigen = $(dirLib)External/Eigen/eigen-3.4.0/
dir_CCfits = $(dirLib)External/CCfits/include
dirCUBA = $(dirLib)External/Cuba-4.2.1/

# Set the Linker Flags to find the OpenMP library and others
LDFLAGS = -L/opt/homebrew/opt/libomp/lib
FLAGS_LIB = -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(dirLib) -L$(dirLib) -lCBL $(LDFLAGS)

# Set the Include Flags to find the OpenMP headers and others
FLAGS_INC = -I/opt/homebrew/include -I/opt/homebrew/opt/libomp/include -I/Users/arun/include/ -I$(dirH) -I$(dirCUBA) -I$(dir_Eigen) -I$(dir_CCfits)

FLAGS_LINK = -shared
ES = so

# Python-related flags
python_version_full := $(wordlist 2,4,$(subst ., ,$(shell $(PY) --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})

ifeq ($(python_version_major),2)
	PYINC = $(shell $(PY) -c 'import sysconfig; print sysconfig.get_config_var("INCLUDEDIR")')
	SWIG_FLAG = -python -c++
	PYVERSION = $(python_version_major).$(python_version_minor)
endif
ifeq ($(python_version_major),3)
	PYINC = $(shell $(PY) -c 'import sysconfig; print(sysconfig.get_config_var("INCLUDEDIR"))')
	PYABI = $(shell $(PY) -c 'import sysconfig; print(sysconfig.get_config_var("ABIFLAGS"))')
	SWIG_FLAG = -python -c++
	PYVERSION = $(python_version_major).$(python_version_minor)$(PYABI)
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
	rm -f *.o *.cxx ModelPY_wrap.cxx modelpy.py modelpy.pyc _modelpy.so *~ \#* temp* core* prior fit sampler valgrind*dat
	rm -rf __pycache__

prior.o: prior.cpp makefile $(dirLib)*.$(ES) 
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c prior.cpp

fit.o: fit.cpp makefile $(dirLib)*.$(ES) 
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c fit.cpp

sampler.o: sampler.cpp makefile $(dirLib)*.$(ES) 
	$(CXX) $(FLAGS0) $(FLAGS) $(FLAGS_INC) -c sampler.cpp

ModelPY_wrap.o: ModelPY_wrap.cxx ModelPY.i ModelPY.h
	$(CXX) $(FLAGST) -Wno-uninitialized -Wno-missing-field-initializers -Wno-unused-parameter $(PFLAGS) -c -fPIC $(FLAGS_INC) ModelPY_wrap.cxx -o ModelPY_wrap.o

ModelPY_wrap.cxx: ModelPY.h ModelPY.i
	$(call colorecho, "\n"Running swig... "\n")
	$(SWIG) $(SWIG_FLAG) -I$(dirH) -I$(dirLib)Python/ ModelPY.i
