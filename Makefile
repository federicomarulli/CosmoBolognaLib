# C++ compiler, main compiler, openMP support required
CXX = g++

# C compiler, used to compile CUBA libraries
CC = gcc

# Fortran 90 compiler, used to compile some external libraries
F = gfortran

# Python, used to compile the python wrapper 
PY = python

# swig, used to create the python wrapper
SWIG = swig3.0

# doxygen, used to create the documentation
Doxygen = doxygen

# GSL installation directories
GSL_VERSION = $(shell gsl-config --version)
dir_INC_GSL = $(shell gsl-config --prefix)/include/
dir_LIB_GSL = $(shell gsl-config --prefix)/lib/

# FFTW installation directories
dir_INC_FFTW =
dir_LIB_FFTW =

# cfitsio installation directories
dir_INC_cfitsio =
dir_LIB_cfitsio =


############################################################################
### hopefully, the user would never modify the makefile after this point ###
############################################################################

Dir_H = Headers/
Dir_CCfits = External/CCfits/
Dir_CUBA = External/Cuba-4.2/
Dir_FFTLOG = External/fftlog-f90-master/
Dir_Eigen = External/eigen-3.3.4/
Dir_Recfast = External/Recfast/

dir_H = $(addprefix $(PWD)/,$(Dir_H))
dir_CCfits = $(addprefix $(PWD)/,$(Dir_CCfits))
dir_CUBA = $(addprefix $(PWD)/,$(Dir_CUBA))
dir_FFTLOG = $(addprefix $(PWD)/,$(Dir_FFTLOG))
dir_Eigen = $(addprefix $(PWD)/, $(Dir_Eigen))
dir_Recfast = $(addprefix $(PWD)/, $(Dir_Recfast))

dir_Python = $(PWD)/Python/

HH = $(dir_H)*.h 

###################
### BASIC FLAGS ###
###################

FLAGS0 = -std=c++11 -fopenmp
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror
FLAGST = $(FLAGS0) $(FLAGS)

FLAGS_INC = -I$(dir_Eigen) -I$(dir_CUBA) -I$(dir_CCfits)include/ -I$(dir_Recfast)include/ -I$(dir_H)

FLAGS_LINK = -shared

#################
### GSL FLAGS ###
#################

#Check the GSL version
gsl_version_full := $(wordlist 1,3,$(subst ., ,$(GSL_VERSION)))
gsl_version_major := $(word 1,${gsl_version_full})
gsl_version_minor := $(word 2,${gsl_version_full})

GSL_VERSION_OK = 0
ifeq ($(gsl_version_major),2)
  ifeq ($(shell test $(gsl_version_minor) -gt 4; echo $$?),0)
    GSL_VERSION_OK = 1
  endif
endif
FLAGS0 := $(FLAGS0) -DGSL_VERSION_OK=$(GSL_VERSION_OK)

# add in FLAGS_INC
ifeq ($(dir_INC_GSL),)
  else
    FLAGS_INC := $(FLAGS_INC) -I$(dir_INC_GSL)
endif

# add search path when linking
FLAGS_GSL = -lgsl -lgslcblas -lm
ifeq ($(dir_LIB_GSL),)
  else 
  FLAGS_GSL := -Wl,-rpath,$(dir_LIB_GSL) -L$(dir_LIB_GSL) $(FLAGS_GSL) 
endif

##################
### FFTW FLAGS ###
##################

# add in FLAGS_INC
ifeq ($(dir_INC_FFTW),)
  else
  FLAGS_INC :=  $(FLAGS_INC) -I$(dir_INC_FFTW)
endif

# add search path when linking
FLAGS_FFTW = -lfftw3 #-lfftw3_omp
ifeq ($(dir_LIB_FFTW),)
  else
  FLAGS_FFTW := -Wl,-rpath,$(dir_LIB_FFTW) -L$(dir_LIB_FFTW) $(FLAGS_FFTW) #-lfftw3_omp
endif


####################
### CCFITS FLAGS ###
####################
FLAGS_CCFITS = -Wl,-rpath,$(dir_CCfits)/lib -L$(dir_CCfits)/lib -lCCfits
CCfits_LIB = $(dir_CCfits)/lib/libCCfits.$(ES)

ifeq ($(dir_INC_cfitsio),)
    CCfits_COMPILE = cd $(dir_CCfits) && tar -xzf CCfits-2.5.tar.gz && cd CCfits && sed -i -e "s/bad_cast/bad_cast\&/g" ColumnT.h && ./configure CXX=$(CXX) --prefix=$(dir_CCfits) && make && make install
  else
    CCfits_COMPILE = cd $(dir_CCfits) && tar -xzf CCfits-2.5.tar.gz && cd CCfits &&  sed -i -e "s/bad_cast/bad_cast\&/g" ColumnT.h && ./configure CXX=$(CXX) --with-cfitsio-include=$(dir_INC_cfitsio) --with-cfitsio-libdir=$(dir_LIB_cfitsio) --prefix=$(dir_CCfits) && make && make install 
endif


####################
### FFTLOG FLAGS ###
####################
FLAGS_FFTLOG = -fPIC -w

#####################
### RECFAST FLAGS ###
#####################
FLAGS_Recfast = -Wall -O3 -fPIC -D RECFASTPPPATH=\"$(PWD)/External/Recfast/\"
FLAGST_Recfast = $(FLAGS0) $(FLAGS_Recfast)

##################
### CUBA FLAGS ###
##################
CUBA_LIB = $(dir_CUBA)libcuba.a
CUBA_COMPILE = cd $(dir_CUBA) && ./configure CC=$(CC) CFLAGS=-fPIC && make lib


####################
### PYTHON FLAGS ###
####################
python_version_full := $(wordlist 2,4,$(subst ., ,$(shell $(PY) --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})

ifeq ($(python_version_major),2)
	PYINC = $(shell $(PY) -c 'from distutils import sysconfig; print sysconfig.get_config_var("INCLUDEDIR")')
	PYLIB = $(shell $(PY) -c 'from distutils import sysconfig; print sysconfig.get_config_var("LIBDIR")')    
	SWIG_FLAG = -python -c++ -threads
	PYVERSION = $(python_version_major).$(python_version_minor)
endif
ifeq ($(python_version_major),3)
	PYINC = $(shell $(PY) -c 'from distutils import sysconfig; print(sysconfig.get_config_var("INCLUDEDIR"))')
	PYLIB = $(shell $(PY) -c 'from distutils import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')    
	SWIG_FLAG = -python -c++ -py3 -threads
	PYVERSION = $(python_version_major).$(python_version_minor)m
endif

PFLAGS = -I$(PYINC)/python$(PYVERSION)

ES = so

Dvar = -DLINUX

SYS:=$(shell uname -s)

ifeq ($(SYS),Darwin)
	Dvar = -DMAC 
	FLAGS_LINK = -dynamiclib -undefined dynamic_lookup
        ES = dylib
	FLAGS_PY = -L$(PYLIB) -lpython$(PYVERSION) -ldl	
endif

####################################################################


##### CBL directories #####

Dir_KERNEL = Kernel/
Dir_WRAP = Wrappers/
Dir_FUNCGRID = FuncGrid/
Dir_FFT = FFT/
Dir_RAN = RandomNumbers/
Dir_FUNC = Func/
Dir_DATA = Data/
Dir_FIELD = Field/
Dir_HIST = Histogram/
Dir_DISTR = Distribution/
Dir_STAT = Statistics/
Dir_COSM = Cosmology/Lib/
Dir_CM = ChainMesh/
Dir_CAT = Catalogue/
Dir_LN = LogNormal/
Dir_NC = Measure/NumberCounts/
Dir_TWOP = Measure/TwoPointCorrelation/
Dir_THREEP = Measure/ThreePointCorrelation/
Dir_MODEL_GLOB = Modelling/Global/
Dir_MODEL_COSM = Modelling/Cosmology/
Dir_MODEL_NC = Modelling/NumberCounts/
Dir_MODEL_TWOP = Modelling/TwoPointCorrelation/
Dir_MODEL_THREEP = Modelling/ThreePointCorrelation/
Dir_GLOB = GlobalFunc/
Dir_READP = ReadParameters/

dir_KERNEL = $(addprefix $(PWD)/,$(Dir_KERNEL))
dir_WRAP = $(addprefix $(PWD)/,$(Dir_WRAP))
dir_FUNCGRID = $(addprefix $(PWD)/,$(Dir_FUNCGRID))
dir_FFT = $(addprefix $(PWD)/,$(Dir_FFT))
dir_RAN = $(addprefix $(PWD)/,$(Dir_RAN))
dir_FUNC = $(addprefix $(PWD)/,$(Dir_FUNC))
dir_DATA = $(addprefix $(PWD)/,$(Dir_DATA))
dir_FIELD = $(addprefix $(PWD)/,$(Dir_FIELD))
dir_HIST = $(addprefix $(PWD)/,$(Dir_HIST))
dir_DISTR = $(addprefix $(PWD)/,$(Dir_DISTR))
dir_STAT = $(addprefix $(PWD)/,$(Dir_STAT))
dir_COSM = $(addprefix $(PWD)/,$(Dir_COSM))
dir_CM = $(addprefix $(PWD)/,$(Dir_CM))
dir_CAT = $(addprefix $(PWD)/,$(Dir_CAT))
dir_LN = $(addprefix $(PWD)/,$(Dir_LN))
dir_NC = $(addprefix $(PWD)/,$(Dir_NC))
dir_TWOP = $(addprefix $(PWD)/,$(Dir_TWOP))
dir_THREEP = $(addprefix $(PWD)/,$(Dir_THREEP))
dir_MODEL_GLOB = $(addprefix $(PWD)/,$(Dir_MODEL_GLOB))
dir_MODEL_COSM = $(addprefix $(PWD)/,$(Dir_MODEL_COSM))
dir_MODEL_NC = $(addprefix $(PWD)/,$(Dir_MODEL_NC))
dir_MODEL_TWOP = $(addprefix $(PWD)/,$(Dir_MODEL_TWOP))
dir_MODEL_THREEP = $(addprefix $(PWD)/,$(Dir_MODEL_THREEP))
dir_GLOB = $(addprefix $(PWD)/,$(Dir_GLOB))
dir_READP = $(addprefix $(PWD)/,$(Dir_READP))


##### FFTlog object files #####

OBJ_FFTLOG = $(dir_FFTLOG)drffti.o $(dir_FFTLOG)drfftb.o $(dir_FFTLOG)drfftf.o $(dir_FFTLOG)fftlog.o $(dir_FFTLOG)cdgamma.o


##### RECfast++ object files #####

OBJ_RECfast = $(dir_Recfast)/src/cosmology.Recfast.o \
	   	$(dir_Recfast)/src/evalode.Recfast.o \
	   	$(dir_Recfast)/src/recombination.Recfast.o \
	   	$(dir_Recfast)/src/ODE_solver.Recfast.o \
	   	$(dir_Recfast)/src/DM_annihilation.Recfast.o \
	  	$(dir_Recfast)/src/Rec_corrs_CT.Recfast.o 


##### CBL object files #####

OBJ_KERNEL = $(dir_KERNEL)Kernel.o

OBJ_WRAP = $(dir_WRAP)EigenWrapper.o $(dir_WRAP)GSLwrapper.o $(dir_WRAP)CUBAwrapper.o $(dir_WRAP)FITSwrapper.o 

OBJ_FUNCGRID = $(dir_FUNCGRID)FuncGrid.o $(dir_FUNCGRID)FuncGrid_Bspline.o

OBJ_FFT = $(OBJ_FFTLOG) $(dir_FFT)FFTlog.o 

OBJ_RAN = $(dir_RAN)RandomNumbers.o

OBJ_FUNC = $(dir_FUNC)Func.o $(dir_FUNC)FuncXi.o $(dir_FUNC)FuncMultipoles.o $(dir_FUNC)SphericalHarmonics_Coefficients.o

OBJ_DATA = $(dir_DATA)Data.o $(dir_DATA)Data1D.o $(dir_DATA)Data1D_collection.o $(dir_DATA)Data2D.o $(dir_DATA)Data1D_extra.o $(dir_DATA)Data2D_extra.o $(dir_DATA)CovarianceMatrix.o $(dir_DATA)TaperedCovarianceMatrix.o $(dir_DATA)Table.o

OBJ_FIELD = $(dir_FIELD)Field3D.o

OBJ_HIST = $(dir_HIST)Histogram.o 

OBJ_DISTR = $(dir_DISTR)Distribution.o 

OBJ_STAT =  $(dir_STAT)Prior.o $(dir_STAT)ModelParameters.o $(dir_STAT)LikelihoodParameters.o $(dir_STAT)PosteriorParameters.o $(dir_STAT)Model.o $(dir_STAT)Model1D.o $(dir_STAT)Model2D.o $(dir_STAT)LikelihoodFunction.o $(dir_STAT)Likelihood.o $(dir_STAT)Chi2.o $(dir_STAT)Sampler.o $(dir_STAT)Posterior.o

OBJ_COSM = $(dir_COSM)Cosmology.o $(dir_COSM)Sigma.o $(dir_COSM)PkXi.o $(dir_COSM)PkXizSpace.o $(dir_COSM)PkXiNonLinear.o $(dir_COSM)MassFunction.o $(dir_COSM)Bias.o $(dir_COSM)RSD.o $(dir_COSM)DensityProfile.o $(dir_COSM)Velocities.o $(dir_COSM)MassGrowth.o $(dir_COSM)NG.o $(dir_COSM)BAO.o $(dir_COSM)SizeFunction.o  $(dir_COSM)3PCF.o $(OBJ_RECfast)

OBJ_CM = $(dir_CM)ChainMesh.o

OBJ_CAT = $(dir_CAT)Object.o $(dir_CAT)Catalogue.o $(dir_CAT)RandomCatalogue.o $(dir_CAT)ChainMesh_Catalogue.o $(dir_CAT)RandomCatalogueVIPERS.o $(dir_CAT)VoidCatalogue.o $(dir_CAT)GadgetCatalogue.o $(dir_CAT)FITSCatalogue.o

OBJ_LN = $(dir_LN)LogNormal.o $(dir_LN)LogNormalFull.o

OBJ_NC = $(dir_NC)NumberCounts.o $(dir_NC)NumberCounts1D.o $(dir_NC)NumberCounts2D.o $(dir_NC)NumberCounts1D_Redshift.o $(dir_NC)NumberCounts1D_Mass.o $(dir_NC)NumberCounts2D_RedshiftMass.o $(dir_NC)NumberCounts1D_Size.o

OBJ_TWOP = $(dir_TWOP)Pair.o $(dir_TWOP)Pair1D.o $(dir_TWOP)Pair2D.o $(dir_TWOP)Pair1D_extra.o $(dir_TWOP)Pair2D_extra.o $(dir_TWOP)TwoPointCorrelation.o $(dir_TWOP)TwoPointCorrelation1D.o $(dir_TWOP)TwoPointCorrelation1D_angular.o $(dir_TWOP)TwoPointCorrelation1D_monopole.o $(dir_TWOP)TwoPointCorrelation2D.o $(dir_TWOP)TwoPointCorrelation2D_cartesian.o $(dir_TWOP)TwoPointCorrelation2D_polar.o $(dir_TWOP)TwoPointCorrelation_projected.o $(dir_TWOP)TwoPointCorrelation_deprojected.o $(dir_TWOP)TwoPointCorrelation_multipoles_direct.o $(dir_TWOP)TwoPointCorrelation_multipoles_integrated.o $(dir_TWOP)TwoPointCorrelation_wedges.o $(dir_TWOP)TwoPointCorrelation1D_filtered.o $(dir_TWOP)TwoPointCorrelationCross.o $(dir_TWOP)TwoPointCorrelationCross1D.o $(dir_TWOP)TwoPointCorrelationCross1D_angular.o $(dir_TWOP)TwoPointCorrelationCross1D_monopole.o

OBJ_THREEP = $(dir_THREEP)Triplet.o $(dir_THREEP)ThreePointCorrelation.o $(dir_THREEP)ThreePointCorrelation_angular_connected.o $(dir_THREEP)ThreePointCorrelation_angular_reduced.o $(dir_THREEP)ThreePointCorrelation_comoving_connected.o $(dir_THREEP)ThreePointCorrelation_comoving_reduced.o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles.o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_single.o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_all.o

OBJ_MODEL_GLOB = $(dir_MODEL_GLOB)Modelling.o  $(dir_MODEL_GLOB)Modelling1D.o $(dir_MODEL_GLOB)Modelling2D.o

OBJ_MODEL_COSM = $(dir_MODEL_COSM)ModelFunction_Cosmology.o $(dir_MODEL_COSM)Modelling_Cosmology.o

OBJ_MODEL_NC =  $(dir_MODEL_NC)Modelling_NumberCounts.o $(dir_MODEL_NC)ModelFunction_NumberCounts.o $(dir_MODEL_NC)Modelling_NumberCounts1D.o $(dir_MODEL_NC)Modelling_NumberCounts2D.o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Redshift.o $(dir_MODEL_NC)Modelling_NumberCounts1D_Redshift.o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Mass.o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Size.o $(dir_MODEL_NC)Modelling_NumberCounts1D_Mass.o $(dir_MODEL_NC)ModelFunction_NumberCounts2D_RedshiftMass.o $(dir_MODEL_NC)Modelling_NumberCounts2D_RedshiftMass.o $(dir_MODEL_NC)Modelling_NumberCounts1D_Size.o

OBJ_MODEL_TWOP = $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.o 

OBJ_MODEL_THREEP = $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.o 

OBJ_GLOB = $(dir_GLOB)FuncCosmology.o $(dir_GLOB)Func.o $(dir_GLOB)SubSample.o $(dir_GLOB)Reconstruction.o $(dir_GLOB)Forecast.o

OBJ_READP = $(dir_READP)ReadParameters.o

OBJ_CBL = $(OBJ_KERNEL) $(OBJ_WRAP) $(OBJ_FUNCGRID) $(OBJ_FFT) $(OBJ_RAN) $(OBJ_FUNC) $(OBJ_DATA) $(OBJ_FIELD) $(OBJ_HIST) $(OBJ_DISTR) $(OBJ_STAT) $(OBJ_COSM) $(OBJ_CM) $(OBJ_CAT) $(OBJ_LN) $(OBJ_NC) $(OBJ_TWOP) $(OBJ_THREEP) $(OBJ_MODEL_GLOB) $(OBJ_MODEL_COSM) $(OBJ_MODEL_NC) $(OBJ_MODEL_TWOP) $(OBJ_MODEL_THREEP) $(OBJ_GLOB) $(OBJ_READP)

OBJ_ALL = $(OBJ_CBL) $(PWD)/External/CAMB/*.o $(PWD)/External/CLASS/*.o $(PWD)/External/mangle/*.o $(PWD)/External/MPTbreeze-v1/*.o $(OBJ_CBL) $(PWD)/External/CPT_Library/*.o


# objects for python compilation -> if OBJ_PYTHON=OBJ_CBL then all the CBL will be converted in python modules

OBJ_PYTHON = $(OBJ_KERNEL) $(OBJ_WRAP) $(OBJ_FUNCGRID) $(OBJ_FFT) $(OBJ_RAN) $(OBJ_FUNC) $(OBJ_DATA) $(OBJ_FIELD) $(OBJ_HIST) $(OBJ_DISTR) $(OBJ_STAT) $(OBJ_COSM) $(OBJ_CM) $(OBJ_CAT) $(OBJ_LN) $(OBJ_NC) $(OBJ_TWOP) $(OBJ_THREEP) $(OBJ_MODEL_GLOB) $(OBJ_MODEL_COSM) $(OBJ_MODEL_NC) $(OBJ_MODEL_TWOP) $(OBJ_MODEL_THREEP) $(OBJ_GLOB) $(OBJ_READP)


##### CBL source files #####

SOURCE_CBL = $(subst .o,.cpp,$(OBJ_CBL))



#################################################################### 

define colorecho
      @tput setaf 3
      @echo $1
      @tput sgr0
endef

ALL:
	make CUBA  
	make CCfits
	make CAMB
	make CLASS
	make MPTbreeze
	make mangle
	make CPT_Library
	$(call colorecho, "\n"Compiling the library: libKERNEL... "\n")
	make -j3 libKERNEL
	$(call colorecho, "\n"Compiling the library: libWRAP... "\n")
	make -j3 libWRAP
	$(call colorecho, "\n"Compiling the library: libFUNCGRID... "\n")
	make -j3 libFUNCGRID
	$(call colorecho, "\n"Compiling the library: libFFT... "\n")
	make -j3 libFFT
	$(call colorecho, "\n"Compiling the library: libRAN... "\n")
	make -j3 libRAN
	$(call colorecho, "\n"Compiling the library: libFUNC... "\n")
	make -j3 libFUNC
	$(call colorecho, "\n"Compiling the library: libDATA... "\n")
	make -j3 libDATA
	$(call colorecho, "\n"Compiling the library: libFIELD... "\n")
	make -j3 libFIELD
	$(call colorecho, "\n"Compiling the library: libHIST... "\n")
	make -j3 libHIST
	$(call colorecho, "\n"Compiling the library: libDISTR... "\n")
	make -j3 libDISTR
	$(call colorecho, "\n"Compiling the library: libSTAT... "\n")
	make -j3 libSTAT
	$(call colorecho, "\n"Compiling the library: libCOSM... "\n")
	make -j3 libCOSM
	$(call colorecho, "\n"Compiling the library: libCM... "\n")
	make -j3 libCM
	$(call colorecho, "\n"Compiling the library: libCAT... "\n")
	make -j3 libCAT
	$(call colorecho, "\n"Compiling the library: libLN... "\n")
	make -j3 libLN
	$(call colorecho, "\n"Compiling the library: libNC... "\n")
	make -j3 libNC
	$(call colorecho, "\n"Compiling the library: libTWOP... "\n")
	make -j3 libTWOP
	$(call colorecho, "\n"Compiling the library: libTHREEP... "\n")
	make -j3 libTHREEP
	$(call colorecho, "\n"Compiling the library: libMODEL_GLOB... "\n")
	make -j3 libMODEL_GLOB
	$(call colorecho, "\n"Compiling the library: libMODEL_COSM... "\n")
	make -j3 libMODEL_COSM
	$(call colorecho, "\n"Compiling the library: libMODEL_NC... "\n")
	make -j3 libMODEL_NC
	$(call colorecho, "\n"Compiling the library: libMODEL_TWOP... "\n")
	make -j3 libMODEL_TWOP
	$(call colorecho, "\n"Compiling the library: libMODEL_THREEP... "\n")
	make -j3 libMODEL_THREEP
	$(call colorecho, "\n"Compiling the library: libGLOB... "\n")
	make -j3 libGLOB
	$(call colorecho, "\n"Compiling the library: libREADP... "\n")
	make -j3 libREADP
	$(call colorecho, "\n"Compiling the full library: libCBL... "\n")
	make -j3 libCBL

libKERNEL: $(OBJ_KERNEL) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libKERNEL.$(ES) $(OBJ_KERNEL) -lgfortran

libWRAP: $(OBJ_WRAP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libWRAP.$(ES) $(OBJ_WRAP) $(CUBA_LIB) $(FLAGS_CCFITS) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL

libFUNCGRID: $(OBJ_FUNCGRID) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libFUNCGRID.$(ES) $(OBJ_FUNCGRID) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP

libFFT: $(OBJ_FFT) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libFFT.$(ES) $(OBJ_FFT) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lgfortran

libRAN: $(OBJ_RAN) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libRAN.$(ES) $(OBJ_RAN) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT

libFUNC: $(OBJ_FUNC) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libFUNC.$(ES) $(OBJ_FUNC) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN

libDATA: $(OBJ_DATA) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libDATA.$(ES) $(OBJ_DATA) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC

libFIELD: $(OBJ_FIELD) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libFIELD.$(ES) $(OBJ_FIELD) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA

libHIST: $(OBJ_HIST) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libHIST.$(ES) $(OBJ_HIST) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD 

libDISTR: $(OBJ_DISTR) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libDISTR.$(ES) $(OBJ_DISTR) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST

libSTAT: $(OBJ_STAT) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libSTAT.$(ES) $(OBJ_STAT) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR 

libCOSM: $(OBJ_COSM) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libCOSM.$(ES) $(OBJ_COSM) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT

libCM: $(OBJ_CM) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libCM.$(ES) $(OBJ_CM) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM

libCAT: $(OBJ_CAT) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libCAT.$(ES) $(OBJ_CAT) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM

libLN: $(OBJ_LN) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libLN.$(ES) $(OBJ_LN) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT 

libNC: $(OBJ_NC) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libNC.$(ES) $(OBJ_NC) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN

libTWOP: $(OBJ_TWOP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libTWOP.$(ES) $(OBJ_TWOP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC

libTHREEP: $(OBJ_THREEP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libTHREEP.$(ES) $(OBJ_THREEP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC -lTWOP

libMODEL_GLOB: $(OBJ_MODEL_GLOB) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_GLOB.$(ES) $(OBJ_MODEL_GLOB) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC -lTWOP -lTHREEP

libMODEL_COSM: $(OBJ_MODEL_COSM) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_COSM.$(ES) $(OBJ_MODEL_COSM) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC -lTWOP -lTHREEP -lMODEL_GLOB

libMODEL_NC: $(OBJ_MODEL_NC) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_NC.$(ES) $(OBJ_MODEL_NC) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC -lTWOP -lTHREEP -lMODEL_GLOB -lMODEL_COSM

libMODEL_TWOP: $(OBJ_MODEL_TWOP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_TWOP.$(ES) $(OBJ_MODEL_TWOP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC -lTWOP -lTHREEP -lMODEL_GLOB -lMODEL_COSM -lMODEL_NC

libMODEL_THREEP: $(OBJ_MODEL_THREEP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_THREEP.$(ES) $(OBJ_MODEL_THREEP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC -lTWOP -lTHREEP -lMODEL_GLOB -lMODEL_COSM -lMODEL_NC -lMODEL_TWOP

libGLOB: $(OBJ_GLOB) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libGLOB.$(ES) $(OBJ_GLOB) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC -lTWOP -lTHREEP -lMODEL_GLOB -lMODEL_COSM -lMODEL_NC -lMODEL_TWOP -lMODEL_THREEP

libREADP: $(OBJ_READP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libREADP.$(ES) $(OBJ_READP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lSTAT -lCOSM -lCM -lCAT -lLN -lNC -lTWOP -lTHREEP -lMODEL_GLOB -lMODEL_COSM -lMODEL_NC -lMODEL_TWOP -lMODEL_THREEP -lGLOB

libCBL: $(OBJ_CBL) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libCBL.$(ES) $(OBJ_CBL) $(CUBA_LIB) $(FLAGS_CCFITS) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -lgfortran 

CUBA: $(CUBA_LIB)

CCfits: $(CCfits_LIB)

CAMB: $(PWD)/External/CAMB/camb

CLASS: $(PWD)/External/CLASS/class

MPTbreeze: $(PWD)/External/MPTbreeze-v1/mptbreeze

mangle: $(PWD)/External/mangle/bin/ransack

venice: $(PWD)/External/VIPERS/venice3.9/venice

CPT_Library: 
	cd $(PWD)/External/CPT_Library/ ; make all

allExamples:
	$(call colorecho, "\n"Compiling the example code: vector.cpp ... "\n")
	cd $(PWD)/Examples/vectors ; make CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: eigen.cpp ... "\n")
	cd $(PWD)/Examples/eigen ; make CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: randomNumbers.cpp ... "\n")
	cd $(PWD)/Examples/randomNumbers ; make CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: randomNumbers_custom.cpp ... "\n")
	cd $(PWD)/Examples/randomNumbers ; make randomNumbers_custom CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)' 
	$(call colorecho, "\n"Compiling the example code: correlated_samples.cpp ... "\n")
	cd $(PWD)/Examples/randomNumbers ; make correlated_samples CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)' 
	$(call colorecho, "\n"Compiling the example code: histogram.cpp ... "\n")
	cd $(PWD)/Examples/histogram ; make histogram CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)' 	
	$(call colorecho, "\n"Compiling the example code: integration_gsl.cpp ... "\n")
	cd $(PWD)/Examples/wrappers ; make integration_gsl CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: minimisation_gsl.cpp ... "\n")
	cd $(PWD)/Examples/wrappers ; make minimisation_gsl CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: integration_cuba.cpp ... "\n")
	cd $(PWD)/Examples/wrappers ; make integration_cuba CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: fits.cpp ... "\n")
	cd $(PWD)/Examples/wrappers ; make fits CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: covsample.cpp ... "\n")
	cd $(PWD)/Examples/covsample ; make CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)' 
	$(call colorecho, "\n"Compiling the example code: cosmology.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make cosmology CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: fsigma8.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make fsigma8 CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: distances.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make distances CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: model_cosmology.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make model_cosmology CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: data1D.cpp ... "\n")
	cd $(PWD)/Examples/data ; make data1D CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: prior.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make prior CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: fit.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make fit CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the wrapper for the example code: fit.py ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make modelpy CXX=$(CXX) SWIG=$(SWIG) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: sampler.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make sampler CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: catalogue.cpp ... "\n")
	cd $(PWD)/Examples/catalogue ; make catalogue CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)' 
	$(call colorecho, "\n"Compiling the example code: numberCounts.cpp ... "\n")
	cd $(PWD)/Examples/numberCounts/codes ; make numberCounts CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: numberCounts_errors.cpp ... "\n")
	cd $(PWD)/Examples/numberCounts/codes ; make numberCounts_errors CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: 2pt_monopole.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_monopole CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: 2pt_monopole_errors.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_monopole_errors CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: 2pt_multipoles.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_multipoles CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: 2pt_2D.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_2D CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: 2pt_projected.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_projected CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: 2pt_angular.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_angular CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: 3pt.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 3pt CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: 3pt_multipoles.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 3pt_multipoles CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: model_2pt_monopole_BAO.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_monopole_BAO CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: model_2pt_monopole_RSD.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_monopole_RSD CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: model_2pt_projected.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_projected CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: model_2pt_multipoles.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_multipoles CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: model_2pt_2D.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_2D CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: model_3pt.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_3pt CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: sizeFunction.cpp ... "\n")
	cd $(PWD)/Examples/cosmicVoids/codes ; make sizeFunction CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: cleanVoidCatalogue.cpp ... "\n")
	cd $(PWD)/Examples/cosmicVoids/codes ; make cleanVoidCatalogue CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: modelling_VoidAbundances ... "\n")
	cd $(PWD)/Examples/cosmicVoids/codes ; make modelling_VoidAbundances CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: readParameterFile.cpp ... "\n")
	cd $(PWD)/Examples/readParameterFile/ ; make CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'

python: $(dir_Python)CBL_wrap.o libCBL $(dir_Python)CBL.i $(PWD)/Makefile
	make ALL
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_CosmoBolognaLib.so  $(dir_Python)CBL_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

documentation:
	rm -rf Doc/html/* Doc/xml/*
	$(Doxygen) Doc/dconfig
	rm -f Doc/doxygen_sqlite3.db 
#	python ../bin/doxy2swig2.py Doc/xml/index.xml Doc/documentation.i
#	python ../doxy2swig/doxy2swig.py Doc/xml/index.xml Doc/documentation.i

doct:
	rm -rf Doc/html/* Doc/xml/*
	doxygen Doc/dconfigT
	rm -f Doc/doxygen_sqlite3.db

cleanExamples:
	cd $(PWD)/Examples/vectors ; make clean && cd ../..
	cd $(PWD)/Examples/eigen ; make clean && cd ../..
	cd $(PWD)/Examples/randomNumbers ; make clean && cd ../..
	cd $(PWD)/Examples/histogram ; make clean && cd ../..
	cd $(PWD)/Examples/wrappers ; make clean && cd ../..
	cd $(PWD)/Examples/covsample ; make clean && cd ../..
	cd $(PWD)/Examples/cosmology ; make clean && cd ../..
	cd $(PWD)/Examples/data ; make clean && cd ../..
	cd $(PWD)/Examples/statistics/codes ; make clean && cd ../..
	cd $(PWD)/Examples/catalogue ; make clean && cd ../..
	cd $(PWD)/Examples/numberCounts/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/clustering/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/cosmicVoids/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/readParameterFile ; make clean && cd ../..
	rm -rf $(PWD)/Examples/cosmology/results* $(PWD)/Examples/histogram/*.dat $(PWD)/Examples/statistics/output/* $(PWD)/Examples/numberCounts/output/* $(PWD)/Examples/clustering/output/* $(PWD)/Examples/cosmicVoids/output/*


cleanpy:
	rm -f $(dir_Python)*~ $(dir_Python)CBL_wrap.o $(dir_Python)CBL_wrap.cxx $(dir_Python)CosmoBolognaLib.py*
	rm -f $(dir_Python)Lib/*~ $(dir_Python)Lib/*.o $(dir_Python)Lib/*.cxx $(dir_Python)Lib/*.py
	rm -f $(dir_Python)CosmoBolognaLib/*CosmoBolognaLib* $(dir_Python)CosmoBolognaLib/*~ $(dir_Python)CosmoBolognaLib/*.pyc
	rm -f $(dir_Python)_CosmoBolognaLib.so $(dir_Python)CosmoBolognaLib.py
	rm -rf $(dir_Python)dist $(dir_Python)build $(dir_Python)CosmoBolognaLib.egg-info $(dir_Python)__pycache__


cleanTEMP:
	rm -f $(OBJ_ALL) core* $(PWD)/*~ $(dir_KERNEL)*~ $(dir_WRAP)*~ $(dir_FUNCGRID)*~ $(dir_FFT)*~ $(dir_RAN)*~ $(dir_FUNC)*~ $(dir_DATA)*~ $(dir_FIELD)*~ $(dir_HIST)*~ $(dir_DISTR)*~ $(dir_STAT)*~ $(dir_COSM)*~ $(dir_CM)*~ $(dir_CAT)*~ $(dir_LN)*~ $(dir_NC)*~ $(dir_TWOP)*~ $(dir_MODEL_GLOB)*~ $(dir_MODEL_COSM)*~ $(dir_MODEL_NC)*~ $(dir_MODEL_TWOP)*~ $(dir_MODEL_THREEP)*~ $(dir_THREEP)*~ $(dir_GLOB)*~ $(dir_READP)*~ $(dir_H)*~ $(PWD)/\#* $(dir_KERNEL)\#* $(dir_WRAP)\#* $(dir_FUNCGRID)\#* $(dir_FFT)\#* $(dir_RAN)\#* $(dir_FUNC)\#* $(dir_DATA)\#* $(dir_FIELD)\#* $(dir_HIST)\#*  $(dir_DISTR)\#* $(dir_STAT)\#* $(dir_COSM)\#* $(dir_CM)\#* $(dir_CAT)\#* $(dir_LN)\#* $(dir_TWOP)\#* $(dir_THREEP)\#* $(dir_MODEL_GLOB)\#* $(dir_MODEL_COSM)\#* $(dir_MODEL_NC)\#* $(dir_MODEL_TWOP)\#* $(dir_MODEL_THREEP)\#* $(dir_GLOB)\#* $(dir_READP)\#* $(dir_H)\#* $(PWD)/Doc/WARNING_LOGFILE* $(PWD)/Doc/*~

clean:
	make cleanExamples
	make cleanTEMP

purge:
	make clean
	rm -f *.$(ES) temp*

purgeALL:
	make purge
	make cleanpy
	rm -rf Doc/html/* Doc/xml/* 
	rm -rf Cosmology/Tables/*
	rm -rf External/EisensteinHu/output_linear/*
	cd External/CAMB ; make clean 
	rm -rf External/CAMB/camb
	rm -rf External/CAMB/output_linear/*
	rm -rf External/CAMB/output_nonlinear/*
	rm -rf External/CAMB/test_*
	rm -rf External/CAMB/NULL*
	rm -rf External/VIPERS/venice3.9/venice
	rm -rf External/mangle/bin
	cd External/mangle/src; make cleaner ; rm -f Makefile libmangle.a; true
	cd External/CLASS/ ; make clean ; rm -rf class libclass.a python/build/* ; true
	rm -rf External/CLASS/output_linear/*
	rm -rf External/CLASS/output_nonlinear/*
	cd External/fftlog-f90-master/ ; make clean ; rm -f fftlog-f90 ; true
	cd External/MPTbreeze-v1/Cuba-1.4/ ; rm -rf config.h config.log config.status demo-fortran.dSYM/ libcuba.a makefile *~ ; true
	rm -rf External/MPTbreeze-v1/mptbreeze
	rm -rf External/MPTbreeze-v1/*~
	rm -rf External/MPTbreeze-v1/output_linear/*
	rm -rf External/MPTbreeze-v1/output_nonlinear/*
	cd External/Cuba-4.2 ; make distclean ; rm -rf config.h config.log config.status demo-fortran.dSYM/ libcuba.a makefile *~ ; true
	cd External/CCfits ; rm -rf bin lib include CCfits *.fit .deps .libs; true
	cd External/Recfast; make tidy
	cd External/CPT_Library ; make clean 

#################################################################### 


$(CUBA_LIB):
	$(CUBA_COMPILE)

$(CCfits_LIB):
	$(CCfits_COMPILE)


#################################################################### 


$(dir_KERNEL)Kernel.o: $(dir_KERNEL)Kernel.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_KERNEL)Kernel.cpp -o $(dir_KERNEL)Kernel.o 


####################################################################

$(dir_WRAP)EigenWrapper.o: $(dir_WRAP)EigenWrapper.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP)EigenWrapper.cpp -o $(dir_WRAP)EigenWrapper.o 

$(dir_WRAP)GSLwrapper.o: $(dir_WRAP)GSLwrapper.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP)GSLwrapper.cpp -o $(dir_WRAP)GSLwrapper.o 

$(dir_WRAP)CUBAwrapper.o: $(dir_WRAP)CUBAwrapper.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP)CUBAwrapper.cpp -o $(dir_WRAP)CUBAwrapper.o

$(dir_WRAP)FITSwrapper.o: $(dir_WRAP)FITSwrapper.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP)FITSwrapper.cpp -o $(dir_WRAP)FITSwrapper.o


#################################################################### 


$(dir_FUNCGRID)FuncGrid.o: $(dir_FUNCGRID)FuncGrid.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNCGRID)FuncGrid.cpp -o $(dir_FUNCGRID)FuncGrid.o

$(dir_FUNCGRID)FuncGrid_Bspline.o: $(dir_FUNCGRID)FuncGrid_Bspline.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNCGRID)FuncGrid_Bspline.cpp -o $(dir_FUNCGRID)FuncGrid_Bspline.o


####################################################################


$(dir_FFT)FFTlog.o: $(OBJ_FFTLOG) $(dir_FFT)FFTlog.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_FFT)FFTlog.cpp -o $(dir_FFT)FFTlog.o 


#################################################################### 


$(dir_FUNC)Func.o: $(dir_FUNC)Func.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Func.cpp -o $(dir_FUNC)Func.o

$(dir_FUNC)FuncXi.o: $(dir_FUNC)FuncXi.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)FuncXi.cpp -o $(dir_FUNC)FuncXi.o

$(dir_FUNC)FuncMultipoles.o: $(dir_FUNC)FuncMultipoles.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)FuncMultipoles.cpp -o $(dir_FUNC)FuncMultipoles.o

$(dir_FUNC)SphericalHarmonics_Coefficients.o: $(dir_FUNC)SphericalHarmonics_Coefficients.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)SphericalHarmonics_Coefficients.cpp -o $(dir_FUNC)SphericalHarmonics_Coefficients.o


#################################################################### 


$(dir_DATA)Data.o: $(dir_DATA)Data.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)Data.cpp -o $(dir_DATA)Data.o

$(dir_DATA)Data1D.o: $(dir_DATA)Data1D.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)Data1D.cpp -o $(dir_DATA)Data1D.o

$(dir_DATA)Data1D_collection.o: $(dir_DATA)Data1D_collection.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)Data1D_collection.cpp -o $(dir_DATA)Data1D_collection.o

$(dir_DATA)Data2D.o: $(dir_DATA)Data2D.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)Data2D.cpp -o $(dir_DATA)Data2D.o

$(dir_DATA)Data1D_extra.o: $(dir_DATA)Data1D_extra.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)Data1D_extra.cpp -o $(dir_DATA)Data1D_extra.o

$(dir_DATA)Data2D_extra.o: $(dir_DATA)Data2D_extra.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)Data2D_extra.cpp -o $(dir_DATA)Data2D_extra.o

$(dir_DATA)CovarianceMatrix.o: $(dir_DATA)CovarianceMatrix.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)CovarianceMatrix.cpp -o $(dir_DATA)CovarianceMatrix.o

$(dir_DATA)TaperedCovarianceMatrix.o: $(dir_DATA)TaperedCovarianceMatrix.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)TaperedCovarianceMatrix.cpp -o $(dir_DATA)TaperedCovarianceMatrix.o

$(dir_DATA)Table.o: $(dir_DATA)Table.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_DATA)Table.cpp -o $(dir_DATA)Table.o

#################################################################### 


$(dir_FIELD)Field3D.o: $(dir_FIELD)Field3D.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FIELD)Field3D.cpp -o $(dir_FIELD)Field3D.o


#################################################################### 


$(dir_HIST)Histogram.o: $(dir_HIST)Histogram.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_HIST)Histogram.cpp -o $(dir_HIST)Histogram.o


#################################################################### 


$(dir_RAN)RandomNumbers.o: $(dir_RAN)RandomNumbers.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_RAN)RandomNumbers.cpp -o $(dir_RAN)RandomNumbers.o 


#################################################################### 


$(dir_DISTR)Distribution.o: $(dir_DISTR)Distribution.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_DISTR)Distribution.cpp -o $(dir_DISTR)Distribution.o 


#################################################################### 


$(dir_STAT)Prior.o: $(dir_STAT)Prior.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Prior.cpp -o $(dir_STAT)Prior.o

$(dir_STAT)ModelParameters.o: $(dir_STAT)ModelParameters.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)ModelParameters.cpp -o $(dir_STAT)ModelParameters.o

$(dir_STAT)LikelihoodParameters.o: $(dir_STAT)LikelihoodParameters.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)LikelihoodParameters.cpp -o $(dir_STAT)LikelihoodParameters.o

$(dir_STAT)PosteriorParameters.o: $(dir_STAT)PosteriorParameters.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)PosteriorParameters.cpp -o $(dir_STAT)PosteriorParameters.o

$(dir_STAT)Model.o: $(dir_STAT)Model.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Model.cpp -o $(dir_STAT)Model.o

$(dir_STAT)Model1D.o: $(dir_STAT)Model1D.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Model1D.cpp -o $(dir_STAT)Model1D.o

$(dir_STAT)Model2D.o: $(dir_STAT)Model2D.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Model2D.cpp -o $(dir_STAT)Model2D.o

$(dir_STAT)LikelihoodFunction.o: $(dir_STAT)LikelihoodFunction.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)LikelihoodFunction.cpp -o $(dir_STAT)LikelihoodFunction.o

$(dir_STAT)Likelihood.o: $(dir_STAT)Likelihood.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Likelihood.cpp -o $(dir_STAT)Likelihood.o

$(dir_STAT)Chi2.o: $(dir_STAT)Chi2.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Chi2.cpp -o $(dir_STAT)Chi2.o

$(dir_STAT)Sampler.o: $(dir_STAT)Sampler.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Sampler.cpp -o $(dir_STAT)Sampler.o

$(dir_STAT)Posterior.o: $(dir_STAT)Posterior.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Posterior.cpp -o $(dir_STAT)Posterior.o


#################################################################### 


$(dir_COSM)Cosmology.o: $(dir_COSM)Cosmology.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Cosmology.cpp -o $(dir_COSM)Cosmology.o

$(dir_COSM)MassFunction.o: $(dir_COSM)MassFunction.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)MassFunction.cpp -o $(dir_COSM)MassFunction.o

$(dir_COSM)SizeFunction.o: $(dir_COSM)SizeFunction.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)SizeFunction.cpp -o $(dir_COSM)SizeFunction.o

$(dir_COSM)PkXi.o: $(dir_COSM)PkXi.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)PkXi.cpp -o $(dir_COSM)PkXi.o

$(dir_COSM)PkXizSpace.o: $(dir_COSM)PkXizSpace.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)PkXizSpace.cpp -o $(dir_COSM)PkXizSpace.o

$(dir_COSM)PkXiNonLinear.o: $(dir_COSM)PkXiNonLinear.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)PkXiNonLinear.cpp -o $(dir_COSM)PkXiNonLinear.o

$(dir_COSM)Bias.o: $(dir_COSM)Bias.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Bias.cpp -o $(dir_COSM)Bias.o

$(dir_COSM)RSD.o: $(dir_COSM)RSD.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)RSD.cpp -o $(dir_COSM)RSD.o

$(dir_COSM)DensityProfile.o: $(dir_COSM)DensityProfile.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)DensityProfile.cpp -o $(dir_COSM)DensityProfile.o

$(dir_COSM)Sigma.o: $(dir_COSM)Sigma.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Sigma.cpp -o $(dir_COSM)Sigma.o

$(dir_COSM)Velocities.o: $(dir_COSM)Velocities.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Velocities.cpp -o $(dir_COSM)Velocities.o

$(dir_COSM)MassGrowth.o: $(dir_COSM)MassGrowth.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)MassGrowth.cpp -o $(dir_COSM)MassGrowth.o

$(dir_COSM)NG.o: $(dir_COSM)NG.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)NG.cpp -o $(dir_COSM)NG.o

$(dir_COSM)BAO.o: $(dir_COSM)BAO.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)BAO.cpp -o $(dir_COSM)BAO.o

$(dir_COSM)3PCF.o: $(dir_COSM)3PCF.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)3PCF.cpp -o $(dir_COSM)3PCF.o


#################################################################### 


$(dir_CM)ChainMesh.o: $(dir_CM)ChainMesh.cpp $(HH) $(PWD)/Makefile  
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CM)ChainMesh.cpp -o $(dir_CM)ChainMesh.o


#################################################################### 


$(dir_CAT)Object.o: $(dir_CAT)Object.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)Object.cpp -o $(dir_CAT)Object.o

$(dir_CAT)Catalogue.o: $(dir_CAT)Catalogue.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)Catalogue.cpp -o $(dir_CAT)Catalogue.o

$(dir_CAT)RandomCatalogue.o: $(dir_CAT)RandomCatalogue.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)RandomCatalogue.cpp -o $(dir_CAT)RandomCatalogue.o

$(dir_CAT)ChainMesh_Catalogue.o: $(dir_CAT)ChainMesh_Catalogue.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)ChainMesh_Catalogue.cpp -o $(dir_CAT)ChainMesh_Catalogue.o

$(dir_CAT)RandomCatalogueVIPERS.o: $(dir_CAT)RandomCatalogueVIPERS.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)RandomCatalogueVIPERS.cpp -o $(dir_CAT)RandomCatalogueVIPERS.o

$(dir_CAT)VoidCatalogue.o: $(dir_CAT)VoidCatalogue.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)VoidCatalogue.cpp -o $(dir_CAT)VoidCatalogue.o

$(dir_CAT)GadgetCatalogue.o: $(dir_CAT)GadgetCatalogue.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)GadgetCatalogue.cpp -o $(dir_CAT)GadgetCatalogue.o

$(dir_CAT)FITSCatalogue.o: $(dir_CAT)FITSCatalogue.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)FITSCatalogue.cpp -o $(dir_CAT)FITSCatalogue.o


#################################################################### 


$(dir_LN)LogNormal.o: $(dir_LN)LogNormal.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_LN)LogNormal.cpp -o $(dir_LN)LogNormal.o

$(dir_LN)LogNormalFull.o: $(dir_LN)LogNormalFull.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_LN)LogNormalFull.cpp -o $(dir_LN)LogNormalFull.o


#################################################################### 


$(dir_NC)NumberCounts.o: $(dir_NC)NumberCounts.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts.cpp -o $(dir_NC)NumberCounts.o

$(dir_NC)NumberCounts1D.o: $(dir_NC)NumberCounts1D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts1D.cpp -o $(dir_NC)NumberCounts1D.o

$(dir_NC)NumberCounts2D.o: $(dir_NC)NumberCounts2D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts2D.cpp -o $(dir_NC)NumberCounts2D.o

$(dir_NC)NumberCounts1D_Redshift.o: $(dir_NC)NumberCounts1D_Redshift.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts1D_Redshift.cpp -o $(dir_NC)NumberCounts1D_Redshift.o

$(dir_NC)NumberCounts1D_Mass.o: $(dir_NC)NumberCounts1D_Mass.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts1D_Mass.cpp -o $(dir_NC)NumberCounts1D_Mass.o

$(dir_NC)NumberCounts2D_RedshiftMass.o: $(dir_NC)NumberCounts2D_RedshiftMass.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts2D_RedshiftMass.cpp -o $(dir_NC)NumberCounts2D_RedshiftMass.o

$(dir_NC)NumberCounts1D_Size.o: $(dir_NC)NumberCounts1D_Size.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts1D_Size.cpp -o $(dir_NC)NumberCounts1D_Size.o

#################################################################### 


$(dir_TWOP)Pair.o: $(dir_TWOP)Pair.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair.cpp -o $(dir_TWOP)Pair.o

$(dir_TWOP)Pair1D.o: $(dir_TWOP)Pair1D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair1D.cpp -o $(dir_TWOP)Pair1D.o

$(dir_TWOP)Pair2D.o: $(dir_TWOP)Pair2D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair2D.cpp -o $(dir_TWOP)Pair2D.o

$(dir_TWOP)Pair1D_extra.o: $(dir_TWOP)Pair1D_extra.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair1D_extra.cpp -o $(dir_TWOP)Pair1D_extra.o

$(dir_TWOP)Pair2D_extra.o: $(dir_TWOP)Pair2D_extra.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair2D_extra.cpp -o $(dir_TWOP)Pair2D_extra.o

$(dir_TWOP)TwoPointCorrelation.o: $(dir_TWOP)TwoPointCorrelation.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation.cpp -o $(dir_TWOP)TwoPointCorrelation.o

$(dir_TWOP)TwoPointCorrelation1D.o: $(dir_TWOP)TwoPointCorrelation1D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation1D.cpp -o $(dir_TWOP)TwoPointCorrelation1D.o

$(dir_TWOP)TwoPointCorrelation2D.o: $(dir_TWOP)TwoPointCorrelation2D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation2D.cpp -o $(dir_TWOP)TwoPointCorrelation2D.o

$(dir_TWOP)TwoPointCorrelation1D_monopole.o: $(dir_TWOP)TwoPointCorrelation1D_monopole.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation1D_monopole.cpp -o $(dir_TWOP)TwoPointCorrelation1D_monopole.o

$(dir_TWOP)TwoPointCorrelation1D_angular.o: $(dir_TWOP)TwoPointCorrelation1D_angular.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation1D_angular.cpp -o $(dir_TWOP)TwoPointCorrelation1D_angular.o

$(dir_TWOP)TwoPointCorrelation2D_cartesian.o: $(dir_TWOP)TwoPointCorrelation2D_cartesian.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation2D_cartesian.cpp -o $(dir_TWOP)TwoPointCorrelation2D_cartesian.o

$(dir_TWOP)TwoPointCorrelation2D_polar.o: $(dir_TWOP)TwoPointCorrelation2D_polar.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation2D_polar.cpp -o $(dir_TWOP)TwoPointCorrelation2D_polar.o

$(dir_TWOP)TwoPointCorrelation_projected.o: $(dir_TWOP)TwoPointCorrelation_projected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_projected.cpp -o $(dir_TWOP)TwoPointCorrelation_projected.o

$(dir_TWOP)TwoPointCorrelation_deprojected.o: $(dir_TWOP)TwoPointCorrelation_deprojected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_deprojected.cpp -o $(dir_TWOP)TwoPointCorrelation_deprojected.o

$(dir_TWOP)TwoPointCorrelation_multipoles_direct.o: $(dir_TWOP)TwoPointCorrelation_multipoles_direct.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_multipoles_direct.cpp -o $(dir_TWOP)TwoPointCorrelation_multipoles_direct.o

$(dir_TWOP)TwoPointCorrelation_multipoles_integrated.o: $(dir_TWOP)TwoPointCorrelation_multipoles_integrated.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_multipoles_integrated.cpp -o $(dir_TWOP)TwoPointCorrelation_multipoles_integrated.o

$(dir_TWOP)TwoPointCorrelation_wedges.o: $(dir_TWOP)TwoPointCorrelation_wedges.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_wedges.cpp -o $(dir_TWOP)TwoPointCorrelation_wedges.o

$(dir_TWOP)TwoPointCorrelation1D_filtered.o: $(dir_TWOP)TwoPointCorrelation1D_filtered.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation1D_filtered.cpp -o $(dir_TWOP)TwoPointCorrelation1D_filtered.o

$(dir_TWOP)TwoPointCorrelationCross.o: $(dir_TWOP)TwoPointCorrelationCross.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelationCross.cpp -o $(dir_TWOP)TwoPointCorrelationCross.o

$(dir_TWOP)TwoPointCorrelationCross1D.o: $(dir_TWOP)TwoPointCorrelationCross1D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelationCross1D.cpp -o $(dir_TWOP)TwoPointCorrelationCross1D.o

$(dir_TWOP)TwoPointCorrelaation2D.o: $(dir_TWOP)TwoPointCorrelation2D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation2D.cpp -o $(dir_TWOP)TwoPointCorrelation2D.o

$(dir_TWOP)TwoPointCorrelationCross1D_monopole.o: $(dir_TWOP)TwoPointCorrelationCross1D_monopole.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelationCross1D_monopole.cpp -o $(dir_TWOP)TwoPointCorrelationCross1D_monopole.o

$(dir_TWOP)TwoPointCorrelationCross1D_angular.o: $(dir_TWOP)TwoPointCorrelationCross1D_angular.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelationCross1D_angular.cpp -o $(dir_TWOP)TwoPointCorrelationCross1D_angular.o


#################################################################### 


$(dir_THREEP)Triplet.o: $(dir_THREEP)Triplet.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)Triplet.cpp -o $(dir_THREEP)Triplet.o

$(dir_THREEP)ThreePointCorrelation.o: $(dir_THREEP)ThreePointCorrelation.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation.cpp -o $(dir_THREEP)ThreePointCorrelation.o

$(dir_THREEP)ThreePointCorrelation_angular_connected.o: $(dir_THREEP)ThreePointCorrelation_angular_connected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_angular_connected.cpp -o $(dir_THREEP)ThreePointCorrelation_angular_connected.o

$(dir_THREEP)ThreePointCorrelation_angular_reduced.o: $(dir_THREEP)ThreePointCorrelation_angular_reduced.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_angular_reduced.cpp -o $(dir_THREEP)ThreePointCorrelation_angular_reduced.o

$(dir_THREEP)ThreePointCorrelation_comoving_connected.o: $(dir_THREEP)ThreePointCorrelation_comoving_connected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_comoving_connected.cpp -o $(dir_THREEP)ThreePointCorrelation_comoving_connected.o

$(dir_THREEP)ThreePointCorrelation_comoving_reduced.o: $(dir_THREEP)ThreePointCorrelation_comoving_reduced.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_comoving_reduced.cpp -o $(dir_THREEP)ThreePointCorrelation_comoving_reduced.o

$(dir_THREEP)ThreePointCorrelation_comoving_multipoles.o: $(dir_THREEP)ThreePointCorrelation_comoving_multipoles.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_comoving_multipoles.cpp -o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles.o

$(dir_THREEP)ThreePointCorrelation_comoving_multipoles_single.o: $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_single.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_single.cpp -o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_single.o

$(dir_THREEP)ThreePointCorrelation_comoving_multipoles_all.o: $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_all.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_all.cpp -o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_all.o


#################################################################### 


$(dir_MODEL_GLOB)Modelling.o: $(dir_MODEL_GLOB)Modelling.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)Modelling.cpp -o $(dir_MODEL_GLOB)Modelling.o

$(dir_MODEL_GLOB)Modelling1D.o: $(dir_MODEL_GLOB)Modelling1D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)Modelling1D.cpp -o $(dir_MODEL_GLOB)Modelling1D.o

$(dir_MODEL_GLOB)Modelling2D.o: $(dir_MODEL_GLOB)Modelling2D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)Modelling2D.cpp -o $(dir_MODEL_GLOB)Modelling2D.o


#################################################################### 


$(dir_MODEL_COSM)Modelling_Cosmology.o: $(dir_MODEL_COSM)Modelling_Cosmology.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_COSM)Modelling_Cosmology.cpp -o $(dir_MODEL_COSM)Modelling_Cosmology.o

$(dir_MODEL_COSM)ModelFunction_Cosmology.o: $(dir_MODEL_COSM)ModelFunction_Cosmology.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_COSM)ModelFunction_Cosmology.cpp -o $(dir_MODEL_COSM)ModelFunction_Cosmology.o


#################################################################### 


$(dir_MODEL_NC)Modelling_NumberCounts.o: $(dir_MODEL_NC)Modelling_NumberCounts.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)Modelling_NumberCounts.cpp -o $(dir_MODEL_NC)Modelling_NumberCounts.o

$(dir_MODEL_NC)ModelFunction_NumberCounts.o: $(dir_MODEL_NC)ModelFunction_NumberCounts.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)ModelFunction_NumberCounts.cpp -o $(dir_MODEL_NC)ModelFunction_NumberCounts.o

$(dir_MODEL_NC)Modelling_NumberCounts1D.o: $(dir_MODEL_NC)Modelling_NumberCounts1D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)Modelling_NumberCounts1D.cpp -o $(dir_MODEL_NC)Modelling_NumberCounts1D.o

$(dir_MODEL_NC)Modelling_NumberCounts2D.o: $(dir_MODEL_NC)Modelling_NumberCounts2D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)Modelling_NumberCounts2D.cpp -o $(dir_MODEL_NC)Modelling_NumberCounts2D.o

$(dir_MODEL_NC)Modelling_NumberCounts1D_Redshift.o: $(dir_MODEL_NC)Modelling_NumberCounts1D_Redshift.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)Modelling_NumberCounts1D_Redshift.cpp -o $(dir_MODEL_NC)Modelling_NumberCounts1D_Redshift.o

$(dir_MODEL_NC)ModelFunction_NumberCounts1D_Redshift.o: $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Redshift.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Redshift.cpp -o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Redshift.o

$(dir_MODEL_NC)Modelling_NumberCounts1D_Mass.o: $(dir_MODEL_NC)Modelling_NumberCounts1D_Mass.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)Modelling_NumberCounts1D_Mass.cpp -o $(dir_MODEL_NC)Modelling_NumberCounts1D_Mass.o

$(dir_MODEL_NC)ModelFunction_NumberCounts1D_Mass.o: $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Mass.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Mass.cpp -o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Mass.o

$(dir_MODEL_NC)ModelFunction_NumberCounts1D_Size.o: $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Size.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Size.cpp -o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Size.o

$(dir_MODEL_NC)ModelFunction_NumberCounts2D_RedshiftMass.o: $(dir_MODEL_NC)ModelFunction_NumberCounts2D_RedshiftMass.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)ModelFunction_NumberCounts2D_RedshiftMass.cpp -o $(dir_MODEL_NC)ModelFunction_NumberCounts2D_RedshiftMass.o

$(dir_MODEL_NC)Modelling_NumberCounts2D_RedshiftMass.o: $(dir_MODEL_NC)Modelling_NumberCounts2D_RedshiftMass.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)Modelling_NumberCounts2D_RedshiftMass.cpp -o $(dir_MODEL_NC)Modelling_NumberCounts2D_RedshiftMass.o

$(dir_MODEL_NC)Modelling_NumberCounts1D_Size.o: $(dir_MODEL_NC)Modelling_NumberCounts1D_Size.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)Modelling_NumberCounts1D_Size.cpp -o $(dir_MODEL_NC)Modelling_NumberCounts1D_Size.o

#################################################################### 


$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.o


#################################################################### 


$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.o

$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.o

$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.o

$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.o

$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.o


#################################################################### 


$(dir_GLOB)FuncCosmology.o: $(dir_GLOB)FuncCosmology.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)FuncCosmology.cpp -o $(dir_GLOB)FuncCosmology.o

$(dir_GLOB)Func.o: $(dir_GLOB)Func.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)Func.cpp -o $(dir_GLOB)Func.o

$(dir_GLOB)SubSample.o: $(dir_GLOB)SubSample.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)SubSample.cpp -o $(dir_GLOB)SubSample.o

$(dir_GLOB)Reconstruction.o: $(dir_GLOB)Reconstruction.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)Reconstruction.cpp -o $(dir_GLOB)Reconstruction.o

$(dir_GLOB)Forecast.o: $(dir_GLOB)Forecast.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)Forecast.cpp -o $(dir_GLOB)Forecast.o


#################################################################### 


$(dir_READP)ReadParameters.o: $(dir_READP)ReadParameters.cpp $(HH) $(PWD)/Makefile 
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_READP)ReadParameters.cpp -o $(dir_READP)ReadParameters.o


#################################################################### 


$(dir_FFTLOG)drffti.o: $(dir_FFTLOG)drffti.f
	$(F) $(FLAGS_FFTLOG) -c $(dir_FFTLOG)drffti.f -o $(dir_FFTLOG)drffti.o

$(dir_FFTLOG)drfftb.o: $(dir_FFTLOG)drfftb.f
	$(F) $(FLAGS_FFTLOG) -c $(dir_FFTLOG)drfftb.f -o $(dir_FFTLOG)drfftb.o

$(dir_FFTLOG)drfftf.o: $(dir_FFTLOG)drfftf.f
	$(F) $(FLAGS_FFTLOG) -c $(dir_FFTLOG)drfftf.f -o $(dir_FFTLOG)drfftf.o

$(dir_FFTLOG)fftlog.o: $(dir_FFTLOG)fftlog.f
	$(F) $(FLAGS_FFTLOG) -c $(dir_FFTLOG)fftlog.f -o $(dir_FFTLOG)fftlog.o

$(dir_FFTLOG)cdgamma.o: $(dir_FFTLOG)cdgamma.f
	$(F) $(FLAGS_FFTLOG) -c $(dir_FFTLOG)cdgamma.f -o $(dir_FFTLOG)cdgamma.o


#################################################################### 


$(dir_Recfast)/src/cosmology.Recfast.o: $(dir_Recfast)/src/cosmology.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)/src/cosmology.Recfast.cpp -o $(dir_Recfast)/src/cosmology.Recfast.o

$(dir_Recfast)/src/evalode.Recfast.o: $(dir_Recfast)/src/evalode.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)/src/evalode.Recfast.cpp -o $(dir_Recfast)/src/evalode.Recfast.o

$(dir_Recfast)/src/recombination.Recfast.o: $(dir_Recfast)/src/recombination.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)/src/recombination.Recfast.cpp -o $(dir_Recfast)/src/recombination.Recfast.o

$(dir_Recfast)/src/ODE_solver.Recfast.o: $(dir_Recfast)/src/ODE_solver.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)/src/ODE_solver.Recfast.cpp -o $(dir_Recfast)/src/ODE_solver.Recfast.o

$(dir_Recfast)/src/DM_annihilation.Recfast.o: $(dir_Recfast)/src/DM_annihilation.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)/src/DM_annihilation.Recfast.cpp -o $(dir_Recfast)/src/DM_annihilation.Recfast.o

$(dir_Recfast)/src/Rec_corrs_CT.Recfast.o: $(dir_Recfast)/src/Rec_corrs_CT.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)/src/Rec_corrs_CT.Recfast.cpp -o $(dir_Recfast)/src/Rec_corrs_CT.Recfast.o


#################################################################### 


$(dir_Python)CBL_wrap.o: $(dir_Python)CBL_wrap.cxx $(dir_Python)CBL.i $(HH) $(PWD)/Makefile 
	$(call colorecho, "\n"Compiling the python wrapper. It may take a few minutes ... "\n")
	$(CXX) $(FLAGST) -Wno-stringop-overflow -Wno-uninitialized $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)CBL_wrap.cxx -o $(dir_Python)CBL_wrap.o

$(dir_Python)CBL_wrap.cxx: $(dir_Python)CBL.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Running swig. It may take a few minutes ... "\n")
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -I$(dir_EH) $(dir_Python)CBL.i


#################################################################### 


$(PWD)/External/CAMB/camb:
	cd $(PWD)/External/CAMB ; make clean && make F90C=$(F) && make clean && cd ../..

$(PWD)/External/CLASS/class:
	cd $(PWD)/External/CLASS ; make clean && make CC=$(CC) OPTFLAG=-O3 && make clean && cd ../..

$(PWD)/External/MPTbreeze-v1/mptbreeze:
	cd External/MPTbreeze-v1/Cuba-1.4/ ; ./configure CC=$(CC) F77=$(F) && make lib && cd ../../../
	cd $(PWD)/External/MPTbreeze-v1 ; ./compile.sh $(F) && cd ../..

$(PWD)/External/mangle/bin/ransack:
	cd $(PWD)/External/mangle/src && mkdir -p ../bin && chmod +x configure && ./configure && make CC=$(CC) F77=$(F) && make clean && cd -

$(PWD)/External/VIPERS/venice3.9/venice:
	cd $(PWD)/External/VIPERS/venice3.9 ; make clean && make CC=$(CC) && make && make clean && cd -

$(PWD)/External/CPT_Library:
	cd $(PWD)/External/CPT_Library ; make clean && make F90C=$(F) && make clean && cd ../..

