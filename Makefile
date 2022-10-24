# C++ compiler, main compiler, openMP support required
CXX = g++

# old C++ compiler for CCfits
CXX_OLD = g++-9

# C compiler, used to compile CUBA libraries
CC = gcc

# Fortran 90 compiler, used to compile some external libraries
F = gfortran
FC = mpif90

# Python, used to compile the python wrapper
PY = python

# swig, used to create the python wrapper
SWIG = swig

# doxygen, used to create the documentation
Doxygen = doxygen

# the path to CosmoBolognaLib
dirCOSMO = -DDIRCOSMO=\"$(PWD)\"

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

# boost installation directories
dir_INC_BOOST = $(shell $(PY) -c 'from distutils import sysconfig; print(sysconfig.get_config_var("INCLUDEDIR"))')
dir_LIB_BOOST =


############################################################################
### hopefully, the user would never modify the makefile after this point ###
############################################################################

EIGEN_VERSION = 3.4.0

Dir_Path = Path/
Dir_H = Headers/
Dir_CCfits = External/CCfits/
Dir_CUBA = External/Cuba-4.2.1/
Dir_FFTLOG = External/fftlog-f90-master/
Dir_CAMB = External/CAMB/fortran/
Dir_CAMButils = External/CAMB/forutils/
Dir_Eigen = External/Eigen/eigen-$(EIGEN_VERSION)/
Dir_Recfast = External/Recfast/

dir_H = $(addprefix $(PWD)/,$(Dir_H))
dir_CCfits = $(addprefix $(PWD)/,$(Dir_CCfits))
dir_CUBA = $(addprefix $(PWD)/,$(Dir_CUBA))
dir_FFTLOG = $(addprefix $(PWD)/,$(Dir_FFTLOG))
dir_CAMB = $(addprefix $(PWD)/,$(Dir_CAMB))
dir_CAMButils = $(addprefix $(PWD)/,$(Dir_CAMButils))
dir_Eigen = $(addprefix $(PWD)/, $(Dir_Eigen))
dir_Recfast = $(addprefix $(PWD)/, $(Dir_Recfast))

dir_Python = $(PWD)/Python/

HH = $(dir_H)*.h


###################
### BASIC FLAGS ###
###################

FLAGS0 = -std=c++14 -fopenmp
FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror 
FLAGST = $(FLAGS0) $(FLAGS) $(dirCOSMO)

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

FLAGS_CCFITS = -Wl,-rpath,$(dir_CCfits)lib -L$(dir_CCfits)lib -lCCfits
CCfits_LIB = $(dir_CCfits)lib/libCCfits.$(ES)

ifeq ($(dir_INC_cfitsio),)
    CCfits_COMPILE = cd $(dir_CCfits) && tar -xzf CCfits-2.6.tar.gz && cd CCfits && ./configure CXX=$(CXX_OLD) --prefix=$(dir_CCfits) CXXFLAGS="-w" && make && make install
  else
    CCfits_COMPILE = cd $(dir_CCfits) && tar -xzf CCfits-2.6.tar.gz && cd CCfits && ./configure CXX=$(CXX_OLD) --with-cfitsio-include=$(dir_INC_cfitsio) --with-cfitsio-libdir=$(dir_LIB_cfitsio) --prefix=$(dir_CCfits) CXXFLAGS="-w" && make && make install
  FLAGS_INC := $(FLAGS_INC) -I$(dir_INC_cfitsio)
endif


###################
### BOOST FLAGS ###
###################

# add in FLAGS_INC
ifeq ($(dir_INC_BOOST),)
  else
  FLAGS_INC :=  $(FLAGS_INC) -I$(dir_INC_BOOST)
endif


####################
### FFTLOG FLAGS ###
####################

GCCVERSION := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 10)
FLAGS_FFTLOG = -fPIC -w
ifeq "$(GCCVERSION)" "1"
    FLAGS_FFTLOG += -fallow-argument-mismatch
endif


##################
### CAMB FLAGS ###
##################

FLAGS_CAMB = -cpp -w -fPIC -ffast-math -MMD -fopenmp -ffree-line-length-none -Ofast -O3 -march=native -lstdc++
ifeq "$(GCCVERSION)" "1"
    FLAGS_CAMB += -fallow-argument-mismatch
endif


#####################
### RECFAST FLAGS ###
#####################

FLAGS_Recfast = -Wall -O3 -fPIC -D RECFASTPPPATH=\"$(PWD)/External/Recfast/\"
FLAGST_Recfast = $(FLAGS0) $(FLAGS_Recfast)


##################
### CUBA FLAGS ###
##################

CUBA_LIB = $(dir_CUBA)libcuba.a
CUBA_COMPILE = cd $(dir_CUBA) && ./configure CC=$(CC) CFLAGS=-fPIC && make lib CC=$(CC)" -w" FC=$(F)" -w"


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
	ABIFLAGS = $(shell $(PY) -c 'from distutils import sysconfig; print(sysconfig.get_config_var("ABIFLAGS"))')
	SWIG_FLAG = -python -c++ -py3 -threads
	PYVERSION = $(python_version_major).$(python_version_minor)$(ABIFLAGS)
endif

SWIG_FLAG_ADD = -Wno-stringop-overflow -Wno-uninitialized -Wno-missing-field-initializers -Wno-unused-parameter -Wno-deprecated-declarations

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

Dir_PATH = Path/
Dir_KERNEL = Kernel/
Dir_WRAP_CAMB = Wrappers/CAMBwrapper/
Dir_WRAP_LIB = Wrappers/Libraries/
Dir_FUNCGRID = FuncGrid/
Dir_FFT = FFT/
Dir_RAN = RandomNumbers/
Dir_CM = ChainMesh/
Dir_FUNC = Func/
Dir_DATA = Data/
Dir_FIELD = Field/
Dir_HIST = Histogram/
Dir_DISTR = Distribution/
Dir_COSM = Cosmology/Lib/
Dir_STAT = Statistics/
Dir_CAT = Catalogue/
Dir_LN = LogNormal/
Dir_NC = Measure/NumberCounts/
Dir_STACKPROFILE = Measure/StackedDensityProfile/
Dir_TWOP = Measure/TwoPointCorrelation/
Dir_THREEP = Measure/ThreePointCorrelation/
Dir_ANGPOW = Measure/AngularPowerSpectrum/
Dir_MODEL_GLOB = Modelling/Global/
Dir_MODEL_COSM = Modelling/Cosmology/
Dir_MODEL_MASSOBSERVABLEREL = Modelling/MassObservableRelation/
Dir_MODEL_DENSITYPROFILE = Modelling/DensityProfile/
Dir_MODEL_NC = Modelling/NumberCounts/
Dir_MODEL_TWOP = Modelling/TwoPointCorrelation/
Dir_MODEL_THREEP = Modelling/ThreePointCorrelation/
Dir_MODEL_ANGPOW = Modelling/AngularPowerSpectrum/
Dir_GLOB = GlobalFunc/
Dir_PARF = ParameterFile/

dir_PATH = $(addprefix $(PWD)/,$(Dir_PATH))
dir_KERNEL = $(addprefix $(PWD)/,$(Dir_KERNEL))
dir_WRAP_CAMB = $(addprefix $(PWD)/,$(Dir_WRAP_CAMB))
dir_WRAP_LIB = $(addprefix $(PWD)/,$(Dir_WRAP_LIB))
dir_FUNCGRID = $(addprefix $(PWD)/,$(Dir_FUNCGRID))
dir_FFT = $(addprefix $(PWD)/,$(Dir_FFT))
dir_RAN = $(addprefix $(PWD)/,$(Dir_RAN))
dir_CM = $(addprefix $(PWD)/,$(Dir_CM))
dir_FUNC = $(addprefix $(PWD)/,$(Dir_FUNC))
dir_DATA = $(addprefix $(PWD)/,$(Dir_DATA))
dir_FIELD = $(addprefix $(PWD)/,$(Dir_FIELD))
dir_HIST = $(addprefix $(PWD)/,$(Dir_HIST))
dir_DISTR = $(addprefix $(PWD)/,$(Dir_DISTR))
dir_COSM = $(addprefix $(PWD)/,$(Dir_COSM))
dir_STAT = $(addprefix $(PWD)/,$(Dir_STAT))
dir_CAT = $(addprefix $(PWD)/,$(Dir_CAT))
dir_LN = $(addprefix $(PWD)/,$(Dir_LN))
dir_NC = $(addprefix $(PWD)/,$(Dir_NC))
dir_STACKPROFILE = $(addprefix $(PWD)/,$(Dir_STACKPROFILE))
dir_TWOP = $(addprefix $(PWD)/,$(Dir_TWOP))
dir_THREEP = $(addprefix $(PWD)/,$(Dir_THREEP))
dir_ANGPOW = $(addprefix $(PWD)/,$(Dir_ANGPOW))
dir_MODEL_GLOB = $(addprefix $(PWD)/,$(Dir_MODEL_GLOB))
dir_MODEL_COSM = $(addprefix $(PWD)/,$(Dir_MODEL_COSM))
dir_MODEL_MASSOBSERVABLEREL = $(addprefix $(PWD)/,$(Dir_MODEL_MASSOBSERVABLEREL))
dir_MODEL_DENSITYPROFILE = $(addprefix $(PWD)/,$(Dir_MODEL_DENSITYPROFILE))
dir_MODEL_NC = $(addprefix $(PWD)/,$(Dir_MODEL_NC))
dir_MODEL_TWOP = $(addprefix $(PWD)/,$(Dir_MODEL_TWOP))
dir_MODEL_THREEP = $(addprefix $(PWD)/,$(Dir_MODEL_THREEP))
dir_MODEL_ANGPOW = $(addprefix $(PWD)/,$(Dir_MODEL_ANGPOW))
dir_GLOB = $(addprefix $(PWD)/,$(Dir_GLOB))
dir_PARF = $(addprefix $(PWD)/,$(Dir_PARF))


##### FFTlog object files #####

OBJ_FFTLOG = $(dir_FFTLOG)drffti.o $(dir_FFTLOG)drfftb.o $(dir_FFTLOG)drfftf.o $(dir_FFTLOG)fftlog.o $(dir_FFTLOG)cdgamma.o


##### CAMB object files #####

OBJ_CAMButils = $(dir_CAMButils)MiscUtils.o $(dir_CAMButils)MpiUtils.o $(dir_CAMButils)ArrayUtils.o $(dir_CAMButils)StringUtils.o $(dir_CAMButils)FileUtils.o $(dir_CAMButils)IniObjects.o $(dir_CAMButils)ObjectLists.o $(dir_CAMButils)Interpolation.o $(dir_CAMButils)MatrixUtils.o $(dir_CAMButils)RandUtils.o $(dir_CAMButils)RangeUtils.o

OBJ_CAMB = $(dir_CAMB)constants.o $(dir_CAMB)config.o $(dir_CAMB)MathUtils.o $(dir_CAMB)classes.o $(dir_CAMB)SourceWindows.o $(dir_CAMB)DarkEnergyInterface.o $(dir_CAMB)massive_neutrinos.o $(dir_CAMB)model.o $(dir_CAMB)results.o $(dir_CAMB)DarkEnergyFluid.o $(dir_CAMB)DarkEnergyPPF.o $(dir_CAMB)bessels.o $(dir_CAMB)DarkAge21cm.o $(dir_CAMB)recfast.o $(dir_CAMB)equations.o $(dir_CAMB)InitialPower.o $(dir_CAMB)halofit.o $(dir_CAMB)lensing.o $(dir_CAMB)SeparableBispectrum.o $(dir_CAMB)cmbmain.o $(dir_CAMB)PowellMinimize.o $(dir_CAMB)DarkEnergyQuintessence.o $(dir_CAMB)reionization.o $(dir_CAMB)SecondOrderPK.o $(dir_CAMB)subroutines.o $(dir_CAMB)camb.o $(dir_CAMB)inidriver.o


##### RECfast++ object files #####

OBJ_RECfast = $(dir_Recfast)src/cosmology.Recfast.o \
	   	$(dir_Recfast)src/evalode.Recfast.o \
	   	$(dir_Recfast)src/recombination.Recfast.o \
	   	$(dir_Recfast)src/ODE_solver.Recfast.o \
	   	$(dir_Recfast)src/DM_annihilation.Recfast.o \
	  	$(dir_Recfast)src/Rec_corrs_CT.Recfast.o


##### CBL object files #####

OBJ_PATH = $(dir_PATH)Path.o

OBJ_KERNEL = $(dir_KERNEL)Kernel.o

OBJ_WRAP_CAMB = $(OBJ_CAMButils) $(OBJ_CAMB) $(dir_WRAP_CAMB)CAMBinterface.o $(dir_WRAP_CAMB)CAMB.o

OBJ_WRAP_LIB = $(dir_WRAP_LIB)EigenWrapper.o $(dir_WRAP_LIB)GSLwrapper.o $(dir_WRAP_LIB)CUBAwrapper.o $(dir_WRAP_LIB)FITSwrapper.o

OBJ_FUNCGRID = $(dir_FUNCGRID)FuncGrid.o $(dir_FUNCGRID)FuncGrid_Bspline.o

OBJ_FFT = $(OBJ_FFTLOG) $(dir_FFT)FFTlog.o

OBJ_RAN = $(dir_RAN)RandomNumbers.o

OBJ_CM = $(dir_CM)ChainMesh.o

OBJ_FUNC = $(dir_FUNC)Func.o $(dir_FUNC)FuncXi.o $(dir_FUNC)FuncMultipoles.o $(dir_FUNC)LegendrePolynomials.o $(dir_FUNC)SphericalHarmonics_Coefficients.o

OBJ_DATA = $(dir_DATA)Data.o $(dir_DATA)Data1D.o $(dir_DATA)Data1D_collection.o $(dir_DATA)Data2D.o $(dir_DATA)Data1D_extra.o $(dir_DATA)Data2D_extra.o $(dir_DATA)CovarianceMatrix.o $(dir_DATA)TaperedCovarianceMatrix.o $(dir_DATA)Table.o

OBJ_FIELD = $(dir_FIELD)Field3D.o

OBJ_HIST = $(dir_HIST)Histogram.o

OBJ_DISTR = $(dir_DISTR)Distribution.o $(dir_DISTR)CombinedDistribution.o

OBJ_COSM = $(dir_COSM)Cosmology.o $(dir_COSM)Sigma.o $(dir_COSM)PkXi.o $(dir_COSM)PkXizSpace.o $(dir_COSM)PkXiNonLinear.o $(dir_COSM)MassFunction.o $(dir_COSM)Bias.o $(dir_COSM)RSD.o $(dir_COSM)Velocities.o $(dir_COSM)MassGrowth.o $(dir_COSM)NG.o $(dir_COSM)BAO.o $(dir_COSM)SizeFunction.o  $(dir_COSM)3PCF.o $(OBJ_RECfast) $(dir_COSM)MassFunction_vector.o $(dir_COSM)Bias_vector.o $(dir_COSM)HaloProfile.o $(dir_COSM)SuperSampleCovariance.o

OBJ_STAT =  $(dir_STAT)Prior.o $(dir_STAT)ModelParameters.o $(dir_STAT)LikelihoodParameters.o $(dir_STAT)PosteriorParameters.o $(dir_STAT)Model.o $(dir_STAT)Model1D.o $(dir_STAT)Model2D.o $(dir_STAT)LikelihoodFunction.o $(dir_STAT)Likelihood.o $(dir_STAT)Chi2.o $(dir_STAT)Sampler.o $(dir_STAT)Posterior.o $(dir_STAT)CombinedPosterior.o

OBJ_CAT = $(dir_CAT)Object.o $(dir_CAT)Catalogue.o $(dir_CAT)RandomCatalogue.o $(dir_CAT)CatalogueChainMesh.o $(dir_CAT)ChainMesh_Catalogue.o $(dir_CAT)RandomCatalogueVIPERS.o $(dir_CAT)VoidCatalogue.o $(dir_CAT)GadgetCatalogue.o $(dir_CAT)FITSCatalogue.o $(dir_CAT)HODCatalogue.o

OBJ_LN = $(dir_LN)LogNormal.o $(dir_LN)LogNormalFull.o

OBJ_NC = $(dir_NC)NumberCounts.o $(dir_NC)NumberCounts1D.o $(dir_NC)NumberCounts2D.o $(dir_NC)NumberCounts1D_Redshift.o $(dir_NC)NumberCounts1D_Mass.o $(dir_NC)NumberCounts1D_MassProxy.o $(dir_NC)NumberCounts2D_RedshiftMass.o $(dir_NC)NumberCounts1D_Size.o

OBJ_STACKPROFILE = $(dir_STACKPROFILE)StackedDensityProfile.o

OBJ_TWOP = $(dir_TWOP)Pair.o $(dir_TWOP)Pair1D.o $(dir_TWOP)Pair2D.o $(dir_TWOP)Pair1D_extra.o $(dir_TWOP)Pair2D_extra.o $(dir_TWOP)TwoPointCorrelation.o $(dir_TWOP)TwoPointCorrelation1D.o $(dir_TWOP)TwoPointCorrelation1D_angular.o $(dir_TWOP)TwoPointCorrelation1D_monopole.o $(dir_TWOP)TwoPointCorrelation2D.o $(dir_TWOP)TwoPointCorrelation2D_cartesian.o $(dir_TWOP)TwoPointCorrelation2D_polar.o $(dir_TWOP)TwoPointCorrelation_projected.o $(dir_TWOP)TwoPointCorrelation_deprojected.o $(dir_TWOP)TwoPointCorrelation_multipoles_direct.o $(dir_TWOP)TwoPointCorrelation_multipoles_integrated.o $(dir_TWOP)TwoPointCorrelation_wedges.o $(dir_TWOP)TwoPointCorrelation1D_filtered.o $(dir_TWOP)TwoPointCorrelationCross.o $(dir_TWOP)TwoPointCorrelationCross1D.o $(dir_TWOP)TwoPointCorrelationCross1D_angular.o $(dir_TWOP)TwoPointCorrelationCross1D_monopole.o

OBJ_THREEP = $(dir_THREEP)Triplet.o $(dir_THREEP)ThreePointCorrelation.o $(dir_THREEP)ThreePointCorrelation_angular_connected.o $(dir_THREEP)ThreePointCorrelation_angular_reduced.o $(dir_THREEP)ThreePointCorrelation_comoving_connected.o $(dir_THREEP)ThreePointCorrelation_comoving_reduced.o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles.o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_single.o $(dir_THREEP)ThreePointCorrelation_comoving_multipoles_all.o

OBJ_ANGPOW = $(dir_ANGPOW)AngularPowerSpectrum.o

OBJ_MODEL_GLOB = $(dir_MODEL_GLOB)Modelling.o  $(dir_MODEL_GLOB)Modelling1D.o $(dir_MODEL_GLOB)Modelling2D.o $(dir_MODEL_GLOB)CombinedModelling.o $(dir_MODEL_GLOB)Modelling_Distribution.o  

OBJ_MODEL_COSM = $(dir_MODEL_COSM)ModelFunction_Cosmology.o $(dir_MODEL_COSM)Modelling_Cosmology.o

OBJ_MODEL_MASSOBSERVABLEREL = $(dir_MODEL_MASSOBSERVABLEREL)Modelling_MassObservableRelation.o

OBJ_MODEL_DENSITYPROFILE = $(dir_MODEL_DENSITYPROFILE)Modelling_DensityProfile.o

OBJ_MODEL_NC =  $(dir_MODEL_NC)Modelling_NumberCounts.o $(dir_MODEL_NC)ModelFunction_NumberCounts.o $(dir_MODEL_NC)Modelling_NumberCounts1D.o $(dir_MODEL_NC)Modelling_NumberCounts2D.o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Redshift.o $(dir_MODEL_NC)Modelling_NumberCounts1D_Redshift.o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Mass.o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_MassProxy.o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_Size.o $(dir_MODEL_NC)Modelling_NumberCounts1D_Mass.o $(dir_MODEL_NC)Modelling_NumberCounts1D_MassProxy.o $(dir_MODEL_NC)ModelFunction_NumberCounts2D_RedshiftMass.o $(dir_MODEL_NC)Modelling_NumberCounts2D_RedshiftMass.o $(dir_MODEL_NC)Modelling_NumberCounts1D_Size.o

OBJ_MODEL_TWOP = $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.o

OBJ_MODEL_THREEP = $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.o

OBJ_MODEL_ANGPOW = $(dir_MODEL_ANGPOW)Modelling_PowerSpectrum_Angular.o $(dir_MODEL_ANGPOW)ModelFunction_PowerSpectrum_Angular.o 

OBJ_GLOB = $(dir_GLOB)FuncCosmology.o $(dir_GLOB)Func.o $(dir_GLOB)SubSample.o $(dir_GLOB)Reconstruction.o $(dir_GLOB)Forecast.o

OBJ_PARF = $(dir_PARF)ReadParameters.o $(dir_PARF)ParameterFile.o

OBJ_CBL = $(OBJ_PATH) $(OBJ_KERNEL) $(OBJ_WRAP_CAMB) $(OBJ_WRAP_LIB) $(OBJ_FUNCGRID) $(OBJ_FFT) $(OBJ_RAN) $(OBJ_FUNC) $(OBJ_DATA) $(OBJ_FIELD) $(OBJ_HIST) $(OBJ_DISTR) $(OBJ_STAT) $(OBJ_COSM) $(OBJ_CM) $(OBJ_CAT) $(OBJ_LN) $(OBJ_NC) $(OBJ_STACKPROFILE) $(OBJ_TWOP) $(OBJ_THREEP) $(OBJ_ANGPOW) $(OBJ_MODEL_GLOB) $(OBJ_MODEL_COSM) $(OBJ_MODEL_MASSOBSERVABLEREL) $(OBJ_MODEL_DENSITYPROFILE) $(OBJ_MODEL_NC) $(OBJ_MODEL_TWOP) $(OBJ_MODEL_THREEP) $(OBJ_MODEL_ANGPOW) $(OBJ_GLOB) $(OBJ_PARF)

OBJ_ALL = $(OBJ_CBL) $(PWD)/External/CAMB/fortran/Release/*.o $(PWD)/External/CLASS/*.o $(PWD)/External/mangle/*.o $(PWD)/External/MPTbreeze-v1/*.o $(OBJ_CBL) $(PWD)/External/CPT_Library/*.o $(PWD)/External/CAMB_SPT_private/*.o $(PWD)/External/MGCAMB/*.o


# objects for python compilation -> if OBJ_PYTHON=OBJ_CBL then all the CBL will be converted in python modules

OBJ_PYTHON = $(OBJ_PATH) $(OBJ_KERNEL) $(OBJ_WRAP_CAMB) $(OBJ_WRAP_LIB) $(OBJ_FUNCGRID) $(OBJ_FFT) $(OBJ_RAN) $(OBJ_FUNC) $(OBJ_DATA) $(OBJ_FIELD) $(OBJ_HIST) $(OBJ_DISTR) $(OBJ_STAT) $(OBJ_COSM) $(OBJ_CM) $(OBJ_CAT) $(OBJ_LN) $(OBJ_NC) $(OBJ_STACKPROFILE) $(OBJ_TWOP) $(OBJ_THREEP) $(OBJ_MODEL_GLOB) $(OBJ_MODEL_COSM) $(OBJ_MODEL_MASSOBSERVABLEREL) $(OBJ_MODEL_DENSITYPROFILE) $(OBJ_MODEL_NC) $(OBJ_MODEL_TWOP) $(OBJ_MODEL_THREEP) $(OBJ_MODEL_ANGPOW) $(OBJ_GLOB) $(OBJ_PARF)


##### CBL source files #####

SOURCE_CBL = $(subst .o,.cpp,$(OBJ_CBL))



####################################################################

define colorecho
      @tput setaf 3
      @echo $1
      @tput sgr0
endef

define insertLine
	grep -qxF $(1) $(2) || echo $(1) >> $(2)
endef

ALL:
	make Eigen
	make CUBA
	make CCfits
	make CAMB
	make CLASS
	make MPTbreeze
	make mangle
	make CPT_Library
	#make CAMB_SPT_private
	make MGCAMB
	$(call colorecho, "\n"Compiling the library: libPATH... "\n")
	make -j3 libPATH
	$(call colorecho, "\n"Compiling the library: libKERNEL... "\n")
	make -j3 libKERNEL
	$(call colorecho, "\n"Compiling the library: libWRAP_CAMB... "\n")
	make libWRAP_CAMB
	$(call colorecho, "\n"Compiling the library: libWRAP_LIB... "\n")
	make -j3 libWRAP_LIB
	$(call colorecho, "\n"Compiling the library: libFUNCGRID... "\n")
	make -j3 libFUNCGRID
	$(call colorecho, "\n"Compiling the library: libFFT... "\n")
	make -j3 libFFT
	$(call colorecho, "\n"Compiling the library: libRAN... "\n")
	make -j3 libRAN
	$(call colorecho, "\n"Compiling the library: libCM... "\n")
	make -j3 libCM
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
	$(call colorecho, "\n"Compiling the library: libCOSM... "\n")
	make -j3 libCOSM
	$(call colorecho, "\n"Compiling the library: libSTAT... "\n")
	make -j3 libSTAT
	$(call colorecho, "\n"Compiling the library: libCAT... "\n")
	make -j3 libCAT
	$(call colorecho, "\n"Compiling the library: libLN... "\n")
	make -j3 libLN
	$(call colorecho, "\n"Compiling the library: libNC... "\n")
	make -j3 libNC
	$(call colorecho, "\n"Compiling the library: libSTACKPROFILE... "\n")
	make -j3 libSTACKPROFILE
	$(call colorecho, "\n"Compiling the library: libTWOP... "\n")
	make -j3 libTWOP
	$(call colorecho, "\n"Compiling the library: libTHREEP... "\n")
	make -j3 libTHREEP
	$(call colorecho, "\n"Compiling the library: libANGPOW... "\n")
	make -j3 libANGPOW
	$(call colorecho, "\n"Compiling the library: libMODEL_GLOB... "\n")
	make -j3 libMODEL_GLOB
	$(call colorecho, "\n"Compiling the library: libMODEL_COSM... "\n")
	make -j3 libMODEL_COSM
	$(call colorecho, "\n"Compiling the library: libMODEL_DENSITYPROFILE... "\n")
	make -j3 libMODEL_DENSITYPROFILE
	$(call colorecho, "\n"Compiling the library: libMODEL_MASSOBSERVABLEREL... "\n")
	make -j3 libMODEL_MASSOBSERVABLEREL
	$(call colorecho, "\n"Compiling the library: libMODEL_NC... "\n")
	make -j3 libMODEL_NC
	$(call colorecho, "\n"Compiling the library: libMODEL_TWOP... "\n")
	make -j3 libMODEL_TWOP
	$(call colorecho, "\n"Compiling the library: libMODEL_THREEP... "\n")
	make -j3 libMODEL_THREEP
	$(call colorecho, "\n"Compiling the library: libMODEL_ANGPOW... "\n")
	make -j3 libMODEL_ANGPOW
	$(call colorecho, "\n"Compiling the library: libGLOB... "\n")
	make -j3 libGLOB
	$(call colorecho, "\n"Compiling the library: libPARF... "\n")
	make -j3 libPARF
	$(call colorecho, "\n"Compiling the full library: libCBL... "\n")
	make -j3 libCBL
	make logo

libPATH: $(OBJ_PATH) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libPATH.$(ES) $(OBJ_PATH) -lgfortran -lgomp

libKERNEL: $(OBJ_KERNEL) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libKERNEL.$(ES) $(OBJ_KERNEL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH

libWRAP_CAMB: $(OBJ_WRAP_CAMB) $(PWD)/Makefile
	$(CXX) $(FLAGS_CAMB) $(FLAGS_LINK) $(OBJ_WRAP_CAMB) -o $(PWD)/libWRAP_CAMB.$(ES) -llapack -lblas -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lgfortran

libWRAP_LIB: $(OBJ_WRAP_LIB) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libWRAP_LIB.$(ES) $(OBJ_WRAP_LIB) $(CUBA_LIB) $(FLAGS_CCFITS) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB

libFUNCGRID: $(OBJ_FUNCGRID) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libFUNCGRID.$(ES) $(OBJ_FUNCGRID) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB

libFFT: $(OBJ_FFT) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libFFT.$(ES) $(OBJ_FFT) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lgfortran

libRAN: $(OBJ_RAN) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libRAN.$(ES) $(OBJ_RAN) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT

libCM: $(OBJ_CM) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libCM.$(ES) $(OBJ_CM) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN

libFUNC: $(OBJ_FUNC) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libFUNC.$(ES) $(OBJ_FUNC) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM

libDATA: $(OBJ_DATA) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libDATA.$(ES) $(OBJ_DATA) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC

libFIELD: $(OBJ_FIELD) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libFIELD.$(ES) $(OBJ_FIELD) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA

libHIST: $(OBJ_HIST) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libHIST.$(ES) $(OBJ_HIST) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD

libDISTR: $(OBJ_DISTR) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libDISTR.$(ES) $(OBJ_DISTR) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST

libCOSM: $(OBJ_COSM) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libCOSM.$(ES) $(OBJ_COSM) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR

libSTAT: $(OBJ_STAT) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libSTAT.$(ES) $(OBJ_STAT) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM

libCAT: $(OBJ_CAT) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libCAT.$(ES) $(OBJ_CAT) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT

libLN: $(OBJ_LN) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libLN.$(ES) $(OBJ_LN) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT

libNC: $(OBJ_NC) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libNC.$(ES) $(OBJ_NC) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN

libSTACKPROFILE: $(OBJ_STACKPROFILE) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libSTACKPROFILE.$(ES) $(OBJ_STACKPROFILE) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC

libTWOP: $(OBJ_TWOP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libTWOP.$(ES) $(OBJ_TWOP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE

libTHREEP: $(OBJ_THREEP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libTHREEP.$(ES) $(OBJ_THREEP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP

libANGPOW: $(OBJ_ANGPOW) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libANGPOW.$(ES) $(OBJ_ANGPOW) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP

libMODEL_GLOB: $(OBJ_MODEL_GLOB) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_GLOB.$(ES) $(OBJ_MODEL_GLOB) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCM -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW

libMODEL_COSM: $(OBJ_MODEL_COSM) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_COSM.$(ES) $(OBJ_MODEL_COSM) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB

libMODEL_DENSITYPROFILE: $(OBJ_MODEL_DENSITYPROFILE) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_DENSITYPROFILE.$(ES) $(OBJ_MODEL_DENSITYPROFILE) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB -lMODEL_COSM

libMODEL_MASSOBSERVABLEREL: $(OBJ_MODEL_MASSOBSERVABLEREL) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_MASSOBSERVABLEREL.$(ES) $(OBJ_MODEL_MASSOBSERVABLEREL) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB -lMODEL_COSM -lMODEL_DENSITYPROFILE

libMODEL_NC: $(OBJ_MODEL_NC) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_NC.$(ES) $(OBJ_MODEL_NC) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB -lMODEL_COSM -lMODEL_DENSITYPROFILE

libMODEL_TWOP: $(OBJ_MODEL_TWOP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_TWOP.$(ES) $(OBJ_MODEL_TWOP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB -lMODEL_COSM -lMODEL_DENSITYPROFILE -lMODEL_NC

libMODEL_THREEP: $(OBJ_MODEL_THREEP) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_THREEP.$(ES) $(OBJ_MODEL_THREEP) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lCM -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB -lMODEL_COSM -lMODEL_DENSITYPROFILE -lMODEL_NC -lMODEL_TWOP

libMODEL_ANGPOW: $(OBJ_MODEL_ANGPOW) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libMODEL_ANGPOW.$(ES) $(OBJ_MODEL_ANGPOW) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lFUNC -lDATA -lFIELD -lCM -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB -lMODEL_COSM -lMODEL_DENSITYPROFILE -lMODEL_NC -lMODEL_TWOP -lMODEL_THREEP

libGLOB: $(OBJ_GLOB) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libGLOB.$(ES) $(OBJ_GLOB) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB -lMODEL_COSM -lMODEL_DENSITYPROFILE -lMODEL_NC -lMODEL_TWOP -lMODEL_THREEP -lMODEL_ANGPOW

libPARF: $(OBJ_PARF) $(PWD)/Makefile
	$(CXX) $(FLAGS_LINK) -o $(PWD)/libPARF.$(ES) $(OBJ_PARF) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lPATH -lKERNEL -lWRAP_CAMB -lWRAP_LIB -lFUNCGRID -lFFT -lRAN -lCM -lFUNC -lDATA -lFIELD -lHIST -lDISTR -lCOSM -lSTAT -lCAT -lLN -lNC -lSTACKPROFILE -lTWOP -lTHREEP -lANGPOW -lMODEL_GLOB -lMODEL_COSM -lMODEL_DENSITYPROFILE -lMODEL_NC -lMODEL_TWOP -lMODEL_THREEP -lMODEL_ANGPOW -lGLOB

libCBL: $(OBJ_CBL) $(PWD)/Makefile
	$(CXX) $(FLAGS_CAMB) $(FLAGS_LINK) $(OBJ_CBL) -o $(PWD)/libCBL.$(ES) -llapack -lblas $(CUBA_LIB) $(FLAGS_CCFITS) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -lgfortran

Eigen: $(PWD)/External/Eigen/eigen-$(EIGEN_VERSION)/Eigen/Dense

CUBA: $(CUBA_LIB)

CCfits: $(CCfits_LIB)

CAMB: $(PWD)/External/CAMB/fortran/camb

CLASS: $(PWD)/External/CLASS/class

MPTbreeze: $(PWD)/External/MPTbreeze-v1/mptbreeze

mangle: $(PWD)/External/mangle/bin/ransack

venice: $(PWD)/External/VIPERS/venice3.9/venice

CPT_Library: $(PWD)/External/CPT_Library/read_pk_sum

CAMB_SPT_private: $(PWD)/External/CAMB_SPT_private/camb

MGCAMB: $(PWD)/External/MGCAMB/camb

logo: $(PWD)/Logo/logo

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
	$(call colorecho, "\n"Compiling the example code: Pk_dynamical_DE.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make Pk_dynamical_DE CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: Pk_BoltzmannSolver.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make Pk_BoltzmannSolver CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: data1D.cpp ... "\n")
	cd $(PWD)/Examples/data ; make data1D CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: prior.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make prior CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: fit.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make fit CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the wrapper for the example code: fit.py ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make modelpy CXX=$(CXX) SWIG=$(SWIG) PY=python3 FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: sampler.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make sampler CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: catalogue.cpp ... "\n")
	cd $(PWD)/Examples/catalogue ; make catalogue CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: catalogueHOD.cpp ... "\n")
	cd $(PWD)/Examples/HOD/codes ; make catalogueHOD CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: lognormal.cpp ... "\n")
	cd $(PWD)/Examples/lognormal/codes ; make lognormal CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
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
	$(call colorecho, "\n"Compiling the example code: power_spectrum_angular.cpp ... "\n")
	cd $(PWD)/Examples/powerSpectrum_angular/codes ; make power_spectrum_angular CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: model_power_spectrum_angular.cpp ... "\n")
	cd $(PWD)/Examples/powerSpectrum_angular/codes ; make model_power_spectrum_angular CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: sizeFunction.cpp ... "\n")
	cd $(PWD)/Examples/cosmicVoids/codes ; make sizeFunction CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: cleanVoidCatalogue.cpp ... "\n")
	cd $(PWD)/Examples/cosmicVoids/codes ; make cleanVoidCatalogue CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: modelling_VoidAbundances ... "\n")
	cd $(PWD)/Examples/cosmicVoids/codes ; make modelling_VoidAbundances CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'
	$(call colorecho, "\n"Compiling the example code: parameterFile.cpp ... "\n")
	cd $(PWD)/Examples/parameterFile/ ; make CXX=$(CXX) FLAGS_INC='$(FLAGS_INC)'


python:	
	make ALL
	make pythonPath
	make pythonKernel
	make pythonWrappers
	make pythonFuncGrid
	make pythonFFT
	make pythonCAMB
	make pythonRandom
	make pythonFunc
	make pythonData
	make pythonField
	make pythonHistogram
	make pythonDistribution
	make pythonStat
	make pythonCosmology
	make pythonChainMesh
	make pythonCatalogue
	make pythonLogNormal
	make pythonMeasure
	make pythonNumberCounts
	make pythonStackedDensityProfile
	make pythonTwoPointCorrelation
	make pythonThreePointCorrelation
	make pythonPowerSpectrumAngular
	make pythonModelling
	make pythonModelling_Cosmology
	make pythonModelling_DensityProfile
	make pythonModelling_MassObservableRelation
	make pythonModelling_NumberCounts
	make pythonModelling_TwoPointCorrelation
	make pythonModelling_ThreePointCorrelation
	make pythonModelling_PowerSpectrumAngular
	make pythonGlobalFunc
	make pythonReadParameters
	$(call colorecho, "\n"The CosmoBolognaLib have been fully wrapped in python!"\n")

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
	cd $(PWD)/Examples/lognormal/codes ; make clean && cd ../..
	cd $(PWD)/Examples/catalogue ; make clean && cd ../..
	cd $(PWD)/Examples/HOD/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/numberCounts/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/clustering/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/powerSpectrum_angular/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/cosmicVoids/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/parameterFile ; make clean && cd ../..
	rm -rf $(PWD)/Examples/cosmology/results* $(PWD)/Examples/histogram/*.dat $(PWD)/Examples/statistics/output/* $(PWD)/Examples/lognormal/output/* $(PWD)/Examples/numberCounts/output/* $(PWD)/Examples/clustering/output/* $(PWD)/Examples/cosmicVoids/output/*


cleanpy:
	rm -rf $(dir_Python)CosmoBolognaLib $(dir_Python)Build $(dir_Python)__pycache__
	rm -f $(dir_Python)*.py $(dir_Python)*.so

cleanTEMP:
	rm -f $(OBJ_ALL) core* $(PWD)/*~ $(dir_PATH)*~ $(dir_KERNEL)*~ $(dir_WRAP_CAMB)*~ $(dir_WRAP_LIB)*~ $(dir_FUNCGRID)*~ $(dir_FFT)*~ $(dir_RAN)*~ $(dir_FUNC)*~ $(dir_DATA)*~ $(dir_FIELD)*~ $(dir_HIST)*~ $(dir_DISTR)*~ $(dir_STAT)*~ $(dir_COSM)*~ $(dir_CM)*~ $(dir_CAT)*~ $(dir_LN)*~ $(dir_NC)*~ $(dir_STACKPROFILE)*~ $(dir_TWOP)*~ $(dir_ANGPOW)*~ $(dir_MODEL_GLOB)*~ $(dir_MODEL_COSM)*~ $(dir_MODEL_MASSOBSERVABLEREL)*~ $(dir_MODEL_DENSITYPROFILE)*~ $(dir_MODEL_NC)*~ $(dir_MODEL_TWOP)*~ $(dir_MODEL_THREEP)*~ $(dir_THREEP)*~ $(dir_GLOB)*~ $(dir_PARF)*~ $(dir_H)*~ $(PWD)/\#* $(dir_KERNEL)\#* $(dir_WRAP_CAMB)\#* $(dir_WRAP_LIB)\#* $(dir_FUNCGRID)\#* $(dir_FFT)\#* $(dir_RAN)\#* $(dir_FUNC)\#* $(dir_DATA)\#* $(dir_FIELD)\#* $(dir_HIST)\#*  $(dir_DISTR)\#* $(dir_STAT)\#* $(dir_COSM)\#* $(dir_CM)\#* $(dir_CAT)\#* $(dir_LN)\#* $(dir_TWOP)\#* $(dir_THREEP)\#* $(dir_MODEL_GLOB)\#* $(dir_MODEL_COSM)\#* $(dir_MODEL_MASSOBSERVABLEREL)\#* $(dir_MODEL_DENSITYPROFILE)\#* $(dir_MODEL_MASSOBSERVABLEREL)\#* $(dir_MODEL_DENSITYPROFILE)\#* $(dir_MODEL_NC)\#* $(dir_MODEL_TWOP)\#* $(dir_MODEL_THREEP)\#* $(dir_MODEL_ANGPOW)\#* $(dir_GLOB)\#* $(dir_PARF)\#* $(dir_H)\#* $(PWD)/Doc/WARNING_LOGFILE* $(PWD)/Doc/*~

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
	rm -rf External/Eigen/eigen-$(EIGEN_VERSION)/
	cd External/CAMB/forutils ; make clean F90C=$(F) ; rm -f *.mod
	cd External/CAMB/fortran ; make clean F90C=$(F) COMPILER=$(F)
	rm -rf External/CAMB/fortran/camb
	rm -rf External/CAMB/fortran/test_*
	rm -rf External/CAMB/output_linear/*
	rm -rf External/CAMB/output_nonlinear/*
	rm -rf External/CAMB/inifiles/test_*
	rm -rf External/CAMB/inifiles/NULL*
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
	rm -rf External/MPTbreeze-v1/output/*
	cd External/Cuba-4.2.1 ; make distclean ; rm -rf config.h config.log config.status demo-fortran.dSYM/ libcuba.a makefile *~ ; true
	cd External/CCfits ; rm -rf bin lib include CCfits *.fit .deps .libs; true
	cd External/Recfast; make tidy
	cd External/CPT_Library ; make clean
	cd External/CAMB_SPT_private ; make clean
	rm -rf External/CAMB_SPT_private/camb
	cd External/MGCAMB ; make clean F90C=$(F) COMPILER=$(F)
	rm -rf External/MGCAMB/camb
	rm -rf External/MGCAMB/test_*
	rm -rf External/MGCAMB/output_linear/*
	rm -rf External/MGCAMB/output_nonlinear/*
	rm -rf External/MGCAMB/test_*
	rm -rf External/MGCAMB/NULL*
	rm -rf Logo/logo
	rm -f Wrapper/CAMBwrapper/*.d Wrapper/CAMBwrapper/*.mod

####################################################################


$(dir_PATH)Path.o: $(dir_PATH)Path.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_PATH)Path.cpp -o $(dir_PATH)Path.o


####################################################################


$(dir_KERNEL)Kernel.o: $(dir_KERNEL)Kernel.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_KERNEL)Kernel.cpp -o $(dir_KERNEL)Kernel.o


####################################################################


$(dir_WRAP_CAMB)CAMB.o: $(dir_WRAP_CAMB)CAMB.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP_CAMB)CAMB.cpp -o $(dir_WRAP_CAMB)CAMB.o


####################################################################


$(dir_WRAP_LIB)EigenWrapper.o: $(dir_WRAP_LIB)EigenWrapper.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP_LIB)EigenWrapper.cpp -o $(dir_WRAP_LIB)EigenWrapper.o

$(dir_WRAP_LIB)GSLwrapper.o: $(dir_WRAP_LIB)GSLwrapper.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP_LIB)GSLwrapper.cpp -o $(dir_WRAP_LIB)GSLwrapper.o

$(dir_WRAP_LIB)CUBAwrapper.o: $(dir_WRAP_LIB)CUBAwrapper.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP_LIB)CUBAwrapper.cpp -o $(dir_WRAP_LIB)CUBAwrapper.o

$(dir_WRAP_LIB)FITSwrapper.o: $(dir_WRAP_LIB)FITSwrapper.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_WRAP_LIB)FITSwrapper.cpp -o $(dir_WRAP_LIB)FITSwrapper.o


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

$(dir_FUNC)LegendrePolynomials.o: $(dir_FUNC)LegendrePolynomials.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)LegendrePolynomials.cpp -o $(dir_FUNC)LegendrePolynomials.o

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


$(dir_CM)ChainMesh.o: $(dir_CM)ChainMesh.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CM)ChainMesh.cpp -o $(dir_CM)ChainMesh.o


####################################################################


$(dir_DISTR)Distribution.o: $(dir_DISTR)Distribution.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_DISTR)Distribution.cpp -o $(dir_DISTR)Distribution.o

$(dir_DISTR)CombinedDistribution.o: $(dir_DISTR)CombinedDistribution.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_DISTR)CombinedDistribution.cpp -o $(dir_DISTR)CombinedDistribution.o


####################################################################


$(dir_COSM)Cosmology.o: $(dir_COSM)Cosmology.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Cosmology.cpp -o $(dir_COSM)Cosmology.o

$(dir_COSM)MassFunction.o: $(dir_COSM)MassFunction.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)MassFunction.cpp -o $(dir_COSM)MassFunction.o

$(dir_COSM)MassFunction_vector.o: $(dir_COSM)MassFunction_vector.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)MassFunction_vector.cpp -o $(dir_COSM)MassFunction_vector.o

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

$(dir_COSM)Bias_vector.o: $(dir_COSM)Bias_vector.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Bias_vector.cpp -o $(dir_COSM)Bias_vector.o

$(dir_COSM)RSD.o: $(dir_COSM)RSD.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)RSD.cpp -o $(dir_COSM)RSD.o

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

$(dir_COSM)HaloProfile.o: $(dir_COSM)HaloProfile.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)HaloProfile.cpp -o $(dir_COSM)HaloProfile.o

$(dir_COSM)SuperSampleCovariance.o: $(dir_COSM)SuperSampleCovariance.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)SuperSampleCovariance.cpp -o $(dir_COSM)SuperSampleCovariance.o


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

$(dir_STAT)CombinedPosterior.o: $(dir_STAT)CombinedPosterior.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)CombinedPosterior.cpp -o $(dir_STAT)CombinedPosterior.o


####################################################################


$(dir_CAT)Object.o: $(dir_CAT)Object.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)Object.cpp -o $(dir_CAT)Object.o

$(dir_CAT)Catalogue.o: $(dir_CAT)Catalogue.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)Catalogue.cpp -o $(dir_CAT)Catalogue.o

$(dir_CAT)RandomCatalogue.o: $(dir_CAT)RandomCatalogue.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)RandomCatalogue.cpp -o $(dir_CAT)RandomCatalogue.o

$(dir_CAT)CatalogueChainMesh.o: $(dir_CAT)CatalogueChainMesh.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)CatalogueChainMesh.cpp -o $(dir_CAT)CatalogueChainMesh.o

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

$(dir_CAT)HODCatalogue.o: $(dir_CAT)HODCatalogue.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)HODCatalogue.cpp -o $(dir_CAT)HODCatalogue.o


####################################################################################################


$(dir_LN)LogNormal.o: $(dir_LN)LogNormal.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_LN)LogNormal.cpp -o $(dir_LN)LogNormal.o

$(dir_LN)LogNormalFull.o: $(dir_LN)LogNormalFull.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_LN)LogNormalFull.cpp -o $(dir_LN)LogNormalFull.o


##################################################################################################


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

$(dir_NC)NumberCounts1D_MassProxy.o: $(dir_NC)NumberCounts1D_MassProxy.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts1D_MassProxy.cpp -o $(dir_NC)NumberCounts1D_MassProxy.o

$(dir_NC)NumberCounts2D_RedshiftMass.o: $(dir_NC)NumberCounts2D_RedshiftMass.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts2D_RedshiftMass.cpp -o $(dir_NC)NumberCounts2D_RedshiftMass.o

$(dir_NC)NumberCounts1D_Size.o: $(dir_NC)NumberCounts1D_Size.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_NC)NumberCounts1D_Size.cpp -o $(dir_NC)NumberCounts1D_Size.o



####################################################################


$(dir_STACKPROFILE)StackedDensityProfile.o: $(dir_STACKPROFILE)StackedDensityProfile.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STACKPROFILE)StackedDensityProfile.cpp -o $(dir_STACKPROFILE)StackedDensityProfile.o


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


$(dir_ANGPOW)AngularPowerSpectrum.o: $(dir_ANGPOW)AngularPowerSpectrum.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_ANGPOW)AngularPowerSpectrum.cpp -o $(dir_ANGPOW)AngularPowerSpectrum.o


####################################################################


$(dir_MODEL_GLOB)Modelling.o: $(dir_MODEL_GLOB)Modelling.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)Modelling.cpp -o $(dir_MODEL_GLOB)Modelling.o

$(dir_MODEL_GLOB)Modelling1D.o: $(dir_MODEL_GLOB)Modelling1D.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)Modelling1D.cpp -o $(dir_MODEL_GLOB)Modelling1D.o

$(dir_MODEL_GLOB)Modelling2D.o: $(dir_MODEL_GLOB)Modelling2D.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)Modelling2D.cpp -o $(dir_MODEL_GLOB)Modelling2D.o

$(dir_MODEL_GLOB)CombinedModelling.o: $(dir_MODEL_GLOB)CombinedModelling.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)CombinedModelling.cpp -o $(dir_MODEL_GLOB)CombinedModelling.o

$(dir_MODEL_GLOB)Modelling_Distribution.o: $(dir_MODEL_GLOB)Modelling_Distribution.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)Modelling_Distribution.cpp -o $(dir_MODEL_GLOB)Modelling_Distribution.o


####################################################################


$(dir_MODEL_COSM)Modelling_Cosmology.o: $(dir_MODEL_COSM)Modelling_Cosmology.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_COSM)Modelling_Cosmology.cpp -o $(dir_MODEL_COSM)Modelling_Cosmology.o

$(dir_MODEL_COSM)ModelFunction_Cosmology.o: $(dir_MODEL_COSM)ModelFunction_Cosmology.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_COSM)ModelFunction_Cosmology.cpp -o $(dir_MODEL_COSM)ModelFunction_Cosmology.o


####################################################################


$(dir_MODEL_MASSOBSERVABLEREL)Modelling_MassObservableRelation.o: $(dir_MODEL_MASSOBSERVABLEREL)Modelling_MassObservableRelation.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_MASSOBSERVABLEREL)Modelling_MassObservableRelation.cpp -o $(dir_MODEL_MASSOBSERVABLEREL)Modelling_MassObservableRelation.o


####################################################################


$(dir_MODEL_DENSITYPROFILE)Modelling_DensityProfile.o: $(dir_MODEL_DENSITYPROFILE)Modelling_DensityProfile.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_DENSITYPROFILE)Modelling_DensityProfile.cpp -o $(dir_MODEL_DENSITYPROFILE)Modelling_DensityProfile.o


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

$(dir_MODEL_NC)Modelling_NumberCounts1D_MassProxy.o: $(dir_MODEL_NC)Modelling_NumberCounts1D_MassProxy.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)Modelling_NumberCounts1D_MassProxy.cpp -o $(dir_MODEL_NC)Modelling_NumberCounts1D_MassProxy.o

$(dir_MODEL_NC)ModelFunction_NumberCounts1D_MassProxy.o: $(dir_MODEL_NC)ModelFunction_NumberCounts1D_MassProxy.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_NC)ModelFunction_NumberCounts1D_MassProxy.cpp -o $(dir_MODEL_NC)ModelFunction_NumberCounts1D_MassProxy.o

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


$(dir_MODEL_ANGPOW)Modelling_PowerSpectrum_Angular.o: $(dir_MODEL_ANGPOW)Modelling_PowerSpectrum_Angular.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_ANGPOW)Modelling_PowerSpectrum_Angular.cpp -o $(dir_MODEL_ANGPOW)Modelling_PowerSpectrum_Angular.o

$(dir_MODEL_ANGPOW)ModelFunction_PowerSpectrum_Angular.o: $(dir_MODEL_ANGPOW)ModelFunction_PowerSpectrum_Angular.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_ANGPOW)ModelFunction_PowerSpectrum_Angular.cpp -o $(dir_MODEL_ANGPOW)ModelFunction_PowerSpectrum_Angular.o


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


$(dir_PARF)ReadParameters.o: $(dir_PARF)ReadParameters.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_PARF)ReadParameters.cpp -o $(dir_PARF)ReadParameters.o

$(dir_PARF)ParameterFile.o: $(dir_PARF)ParameterFile.cpp $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_PARF)ParameterFile.cpp -o $(dir_PARF)ParameterFile.o


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


$(dir_WRAP_CAMB)CAMBinterface.o: $(dir_WRAP_CAMB)CAMBinterface.f90
	$(F) -I$(dir_CAMB) -J$(dir_WRAP_CAMB) $(FLAGS_CAMB) -c $(dir_WRAP_CAMB)CAMBinterface.f90 -o $(dir_WRAP_CAMB)CAMBinterface.o

$(dir_CAMButils)ArrayUtils.o: $(dir_CAMButils)ArrayUtils.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)ArrayUtils.f90 -o $(dir_CAMButils)ArrayUtils.o

$(dir_CAMButils)FileUtils.o: $(dir_CAMButils)FileUtils.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)FileUtils.f90 -o $(dir_CAMButils)FileUtils.o

$(dir_CAMButils)IniObjects.o: $(dir_CAMButils)IniObjects.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)IniObjects.f90 -o $(dir_CAMButils)IniObjects.o

$(dir_CAMButils)Interpolation.o: $(dir_CAMButils)Interpolation.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)Interpolation.f90 -o $(dir_CAMButils)Interpolation.o

$(dir_CAMButils)MatrixUtils.o: $(dir_CAMButils)MatrixUtils.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)MatrixUtils.f90 -o $(dir_CAMButils)MatrixUtils.o

$(dir_CAMButils)MiscUtils.o: $(dir_CAMButils)MiscUtils.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)MiscUtils.f90 -o $(dir_CAMButils)MiscUtils.o

$(dir_CAMButils)MpiUtils.o: $(dir_CAMButils)MpiUtils.f90
	$(FC) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)MpiUtils.f90 -o $(dir_CAMButils)MpiUtils.o

$(dir_CAMButils)ObjectLists.o: $(dir_CAMButils)ObjectLists.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)ObjectLists.f90 -o $(dir_CAMButils)ObjectLists.o

$(dir_CAMButils)RandUtils.o: $(dir_CAMButils)RandUtils.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)RandUtils.f90 -o $(dir_CAMButils)RandUtils.o

$(dir_CAMButils)RangeUtils.o: $(dir_CAMButils)RangeUtils.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)RangeUtils.f90 -o $(dir_CAMButils)RangeUtils.o

$(dir_CAMButils)StringUtils.o: $(dir_CAMButils)StringUtils.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMButils)StringUtils.f90 -o $(dir_CAMButils)StringUtils.o

$(dir_CAMB)bessels.o: $(dir_CAMB)bessels.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMButils) $(FLAGS_CAMB) -c $(dir_CAMB)bessels.f90 -o $(dir_CAMB)bessels.o

$(dir_CAMB)classes.o: $(dir_CAMB)classes.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)classes.f90 -o $(dir_CAMB)classes.o

$(dir_CAMB)cmbmain.o: $(dir_CAMB)cmbmain.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)cmbmain.f90 -o $(dir_CAMB)cmbmain.o

$(dir_CAMB)config.o: $(dir_CAMB)config.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)config.f90 -o $(dir_CAMB)config.o

$(dir_CAMB)constants.o: $(dir_CAMB)constants.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)constants.f90 -o $(dir_CAMB)constants.o

$(dir_CAMB)DarkAge21cm.o: $(dir_CAMB)DarkAge21cm.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)DarkAge21cm.f90 -o $(dir_CAMB)DarkAge21cm.o

$(dir_CAMB)DarkEnergyFluid.o: $(dir_CAMB)DarkEnergyFluid.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)DarkEnergyFluid.f90 -o $(dir_CAMB)DarkEnergyFluid.o

$(dir_CAMB)DarkEnergyInterface.o: $(dir_CAMB)DarkEnergyInterface.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)DarkEnergyInterface.f90 -o $(dir_CAMB)DarkEnergyInterface.o

$(dir_CAMB)DarkEnergyPPF.o: $(dir_CAMB)DarkEnergyPPF.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)DarkEnergyPPF.f90 -o $(dir_CAMB)DarkEnergyPPF.o

$(dir_CAMB)DarkEnergyQuintessence.o: $(dir_CAMB)DarkEnergyQuintessence.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)DarkEnergyQuintessence.f90 -o $(dir_CAMB)DarkEnergyQuintessence.o

$(dir_CAMB)equations.o: $(dir_CAMB)equations.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)equations.f90 -o $(dir_CAMB)equations.o

$(dir_CAMB)halofit.o: $(dir_CAMB)halofit.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)halofit.f90 -o $(dir_CAMB)halofit.o

$(dir_CAMB)inidriver.o: $(dir_CAMB)inidriver.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)inidriver.f90 -o $(dir_CAMB)inidriver.o

$(dir_CAMB)InitialPower.o: $(dir_CAMB)InitialPower.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)InitialPower.f90 -o $(dir_CAMB)InitialPower.o

$(dir_CAMB)lensing.o: $(dir_CAMB)lensing.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)lensing.f90 -o $(dir_CAMB)lensing.o

$(dir_CAMB)massive_neutrinos.o: $(dir_CAMB)massive_neutrinos.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)massive_neutrinos.f90 -o $(dir_CAMB)massive_neutrinos.o

$(dir_CAMB)MathUtils.o: $(dir_CAMB)MathUtils.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)MathUtils.f90 -o $(dir_CAMB)MathUtils.o

$(dir_CAMB)model.o: $(dir_CAMB)model.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)model.f90 -o $(dir_CAMB)model.o

$(dir_CAMB)PowellMinimize.o: $(dir_CAMB)PowellMinimize.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)PowellMinimize.f90 -o $(dir_CAMB)PowellMinimize.o

$(dir_CAMB)recfast.o: $(dir_CAMB)recfast.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)recfast.f90 -o $(dir_CAMB)recfast.o

$(dir_CAMB)reionization.o: $(dir_CAMB)reionization.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)reionization.f90 -o $(dir_CAMB)reionization.o

$(dir_CAMB)results.o: $(dir_CAMB)results.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)results.f90 -o $(dir_CAMB)results.o

$(dir_CAMB)SecondOrderPK.o: $(dir_CAMB)SecondOrderPK.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)SecondOrderPK.f90 -o $(dir_CAMB)SecondOrderPK.o

$(dir_CAMB)SeparableBispectrum.o: $(dir_CAMB)SeparableBispectrum.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)SeparableBispectrum.f90 -o $(dir_CAMB)SeparableBispectrum.o

$(dir_CAMB)SourceWindows.o: $(dir_CAMB)SourceWindows.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)SourceWindows.f90 -o $(dir_CAMB)SourceWindows.o

$(dir_CAMB)subroutines.o: $(dir_CAMB)subroutines.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)subroutines.f90 -o $(dir_CAMB)subroutines.o

$(dir_CAMB)camb.o: $(dir_CAMB)camb.f90
	$(F) -I$(dir_CAMButils) -J$(dir_CAMB) $(FLAGS_CAMB) -c $(dir_CAMB)camb.f90 -o $(dir_CAMB)camb.o


####################################################################


$(dir_Recfast)src/cosmology.Recfast.o: $(dir_Recfast)src/cosmology.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)src/cosmology.Recfast.cpp -o $(dir_Recfast)src/cosmology.Recfast.o

$(dir_Recfast)src/evalode.Recfast.o: $(dir_Recfast)src/evalode.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)src/evalode.Recfast.cpp -o $(dir_Recfast)src/evalode.Recfast.o

$(dir_Recfast)src/recombination.Recfast.o: $(dir_Recfast)src/recombination.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)src/recombination.Recfast.cpp -o $(dir_Recfast)src/recombination.Recfast.o

$(dir_Recfast)src/ODE_solver.Recfast.o: $(dir_Recfast)src/ODE_solver.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)src/ODE_solver.Recfast.cpp -o $(dir_Recfast)src/ODE_solver.Recfast.o

$(dir_Recfast)src/DM_annihilation.Recfast.o: $(dir_Recfast)src/DM_annihilation.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)src/DM_annihilation.Recfast.cpp -o $(dir_Recfast)src/DM_annihilation.Recfast.o

$(dir_Recfast)src/Rec_corrs_CT.Recfast.o: $(dir_Recfast)src/Rec_corrs_CT.Recfast.cpp
	$(CXX) $(FLAGST_Recfast) -c -I$(dir_Recfast)include/ $(dir_Recfast)src/Rec_corrs_CT.Recfast.cpp -o $(dir_Recfast)src/Rec_corrs_CT.Recfast.o


####################################################################


pythonPath: $(dir_Python)/Build $(dir_Python)Build/Path_wrap.o $(dir_Python)Lib/Path.i
	$(call insertLine, "from cblPath import *",  $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Path_wrap.cxx: $(dir_Python)Lib/Path.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Path with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Path_wrap.cxx $(dir_Python)Lib/Path.i

$(dir_Python)Build/Path_wrap.o: $(dir_Python)Build/Path_wrap.cxx $(dir_Python)Lib/Path.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Path_wrap.cxx -o $(dir_Python)Build/Path_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblPath.so $(dir_Python)Build/Path_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonKernel: $(dir_Python)/Build $(dir_Python)Build/Kernel_wrap.o $(dir_Python)Lib/Kernel.i
	$(call insertLine, "from cblKernel import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Kernel_wrap.cxx: $(dir_Python)Lib/Kernel.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Kernel with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Kernel_wrap.cxx $(dir_Python)Lib/Kernel.i

$(dir_Python)Build/Kernel_wrap.o: $(dir_Python)Build/Kernel_wrap.cxx $(dir_Python)Lib/Kernel.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Kernel_wrap.cxx -o $(dir_Python)Build/Kernel_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblKernel.so $(dir_Python)Build/Kernel_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonWrappers: $(dir_Python)/Build $(dir_Python)Build/Wrappers_wrap.o $(dir_Python)Lib/Wrappers.i
	$(call insertLine, "from cblWrappers import *",  $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Wrappers_wrap.cxx: $(dir_Python)Lib/Wrappers.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Wrappers with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Wrappers_wrap.cxx $(dir_Python)Lib/Wrappers.i

$(dir_Python)Build/Wrappers_wrap.o: $(dir_Python)Build/Wrappers_wrap.cxx $(dir_Python)Lib/Wrappers.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Wrappers_wrap.cxx -o $(dir_Python)Build/Wrappers_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblWrappers.so $(dir_Python)Build/Wrappers_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonFuncGrid: $(dir_Python)/Build $(dir_Python)Build/FuncGrid_wrap.o $(dir_Python)Lib/FuncGrid.i
	$(call insertLine, "from cblFuncGrid import *", $(dir_Python)CosmoBolognaLib/__init__.py)	

$(dir_Python)Build/FuncGrid_wrap.cxx: $(dir_Python)Lib/FuncGrid.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping FuncGrid with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/FuncGrid_wrap.cxx $(dir_Python)Lib/FuncGrid.i

$(dir_Python)Build/FuncGrid_wrap.o: $(dir_Python)Build/FuncGrid_wrap.cxx $(dir_Python)Lib/FuncGrid.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/FuncGrid_wrap.cxx -o $(dir_Python)Build/FuncGrid_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblFuncGrid.so $(dir_Python)Build/FuncGrid_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonFFT: $(dir_Python)/Build $(dir_Python)Build/FFT_wrap.o $(dir_Python)Lib/FFT.i
	$(call insertLine, "from cblFFT import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/FFT_wrap.cxx: $(dir_Python)Lib/FFT.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping FFT with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/FFT_wrap.cxx $(dir_Python)Lib/FFT.i

$(dir_Python)Build/FFT_wrap.o: $(dir_Python)Build/FFT_wrap.cxx $(dir_Python)Lib/FFT.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/FFT_wrap.cxx -o $(dir_Python)Build/FFT_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblFFT.so $(dir_Python)Build/FFT_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonCAMB: $(dir_Python)/Build $(dir_Python)Build/CAMB_wrap.o $(dir_Python)Lib/CAMB.i
	$(call insertLine, "from cblCAMB import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/CAMB_wrap.cxx: $(dir_Python)Lib/CAMB.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping CAMB with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/CAMB_wrap.cxx $(dir_Python)Lib/CAMB.i

$(dir_Python)Build/CAMB_wrap.o: $(dir_Python)Build/CAMB_wrap.cxx $(dir_Python)Lib/CAMB.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/CAMB_wrap.cxx -o $(dir_Python)Build/CAMB_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblCAMB.so $(dir_Python)Build/CAMB_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonRandom: $(dir_Python)/Build $(dir_Python)Build/Random_wrap.o $(dir_Python)Lib/Random.i
	$(call insertLine, "from cblRandom import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Random_wrap.cxx: $(dir_Python)Lib/Random.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Random with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Random_wrap.cxx $(dir_Python)Lib/Random.i

$(dir_Python)Build/Random_wrap.o: $(dir_Python)Build/Random_wrap.cxx $(dir_Python)Lib/Random.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Random_wrap.cxx -o $(dir_Python)Build/Random_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblRandom.so $(dir_Python)Build/Random_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonFunc: $(dir_Python)/Build $(dir_Python)Build/Func_wrap.o $(dir_Python)Lib/Func.i
	$(call insertLine, "from cblFunc import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Func_wrap.cxx: $(dir_Python)Lib/Func.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Func with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Func_wrap.cxx $(dir_Python)Lib/Func.i

$(dir_Python)Build/Func_wrap.o: $(dir_Python)Build/Func_wrap.cxx $(dir_Python)Lib/Func.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Func_wrap.cxx -o $(dir_Python)Build/Func_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblFunc.so $(dir_Python)Build/Func_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonData: $(dir_Python)/Build $(dir_Python)Build/Data_wrap.o $(dir_Python)Lib/Data.i
	$(call insertLine, "from cblData import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Data_wrap.cxx: $(dir_Python)Lib/Data.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Data with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Data_wrap.cxx $(dir_Python)Lib/Data.i

$(dir_Python)Build/Data_wrap.o: $(dir_Python)Build/Data_wrap.cxx $(dir_Python)Lib/Data.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Data_wrap.cxx -o $(dir_Python)Build/Data_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblData.so $(dir_Python)Build/Data_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonField: $(dir_Python)/Build $(dir_Python)Build/Field_wrap.o $(dir_Python)Lib/Field.i
	$(call insertLine, "from cblField import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Field_wrap.cxx: $(dir_Python)Lib/Field.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Field with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Field_wrap.cxx $(dir_Python)Lib/Field.i

$(dir_Python)Build/Field_wrap.o: $(dir_Python)Build/Field_wrap.cxx $(dir_Python)Lib/Field.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Field_wrap.cxx -o $(dir_Python)Build/Field_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblField.so $(dir_Python)Build/Field_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonHistogram: $(dir_Python)/Build $(dir_Python)Build/Histogram_wrap.o $(dir_Python)Lib/Histogram.i
	 $(call insertLine, "from cblHistogram import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Histogram_wrap.cxx: $(dir_Python)Lib/Histogram.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Histogram with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Histogram_wrap.cxx $(dir_Python)Lib/Histogram.i

$(dir_Python)Build/Histogram_wrap.o: $(dir_Python)Build/Histogram_wrap.cxx $(dir_Python)Lib/Histogram.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Histogram_wrap.cxx -o $(dir_Python)Build/Histogram_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblHistogram.so $(dir_Python)Build/Histogram_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonDistribution: $(dir_Python)/Build $(dir_Python)Build/Distribution_wrap.o $(dir_Python)Lib/Distribution.i
	$(call insertLine, "from cblDistribution import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Distribution_wrap.cxx: $(dir_Python)Lib/Distribution.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Distribution with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Distribution_wrap.cxx $(dir_Python)Lib/Distribution.i

$(dir_Python)Build/Distribution_wrap.o: $(dir_Python)Build/Distribution_wrap.cxx $(dir_Python)Lib/Distribution.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Distribution_wrap.cxx -o $(dir_Python)Build/Distribution_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblDistribution.so $(dir_Python)Build/Distribution_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonStat: $(dir_Python)/Build $(dir_Python)Build/Stat_wrap.o $(dir_Python)Lib/Stat.i
	$(call insertLine, "from cblStat import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Stat_wrap.cxx: $(dir_Python)Lib/Stat.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Stat with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Stat_wrap.cxx $(dir_Python)Lib/Stat.i

$(dir_Python)Build/Stat_wrap.o: $(dir_Python)Build/Stat_wrap.cxx $(dir_Python)Lib/Stat.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Stat_wrap.cxx -o $(dir_Python)Build/Stat_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblStat.so $(dir_Python)Build/Stat_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonCosmology: $(dir_Python)/Build $(dir_Python)Build/Cosmology_wrap.o $(dir_Python)Lib/Cosmology.i
	$(call insertLine, "from cblCosmology import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Cosmology_wrap.cxx: $(dir_Python)Lib/Cosmology.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Cosmology with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Cosmology_wrap.cxx $(dir_Python)Lib/Cosmology.i

$(dir_Python)Build/Cosmology_wrap.o: $(dir_Python)Build/Cosmology_wrap.cxx $(dir_Python)Lib/Cosmology.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Cosmology_wrap.cxx -o $(dir_Python)Build/Cosmology_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblCosmology.so $(dir_Python)Build/Cosmology_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonChainMesh: $(dir_Python)/Build $(dir_Python)Build/ChainMesh_wrap.o $(dir_Python)Lib/ChainMesh.i
	$(call insertLine, "from cblChainMesh import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/ChainMesh_wrap.cxx: $(dir_Python)Lib/ChainMesh.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping ChainMesh with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/ChainMesh_wrap.cxx $(dir_Python)Lib/ChainMesh.i

$(dir_Python)Build/ChainMesh_wrap.o: $(dir_Python)Build/ChainMesh_wrap.cxx $(dir_Python)Lib/ChainMesh.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/ChainMesh_wrap.cxx -o $(dir_Python)Build/ChainMesh_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblChainMesh.so $(dir_Python)Build/ChainMesh_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonCatalogue: $(dir_Python)/Build $(dir_Python)Build/Catalogue_wrap.o $(dir_Python)Lib/Catalogue.i
	$(call insertLine, "from cblCatalogue import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Catalogue_wrap.cxx: $(dir_Python)Lib/Catalogue.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Catalogue with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Catalogue_wrap.cxx $(dir_Python)Lib/Catalogue.i

$(dir_Python)Build/Catalogue_wrap.o: $(dir_Python)Build/Catalogue_wrap.cxx $(dir_Python)Lib/Catalogue.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Catalogue_wrap.cxx -o $(dir_Python)Build/Catalogue_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblCatalogue.so $(dir_Python)Build/Catalogue_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonLogNormal: $(dir_Python)/Build $(dir_Python)Build/LogNormal_wrap.o $(dir_Python)Lib/LogNormal.i
	$(call insertLine, "from cblLogNormal import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/LogNormal_wrap.cxx: $(dir_Python)Lib/LogNormal.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping LogNormal with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/LogNormal_wrap.cxx $(dir_Python)Lib/LogNormal.i

$(dir_Python)Build/LogNormal_wrap.o: $(dir_Python)Build/LogNormal_wrap.cxx $(dir_Python)Lib/LogNormal.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/LogNormal_wrap.cxx -o $(dir_Python)Build/LogNormal_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblLogNormal.so $(dir_Python)Build/LogNormal_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonMeasure: $(dir_Python)/Build $(dir_Python)Build/Measure_wrap.o $(dir_Python)Lib/Measure.i
	$(call insertLine, "from cblMeasure import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Measure_wrap.cxx: $(dir_Python)Lib/Measure.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Measure with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Measure_wrap.cxx $(dir_Python)Lib/Measure.i

$(dir_Python)Build/Measure_wrap.o: $(dir_Python)Build/Measure_wrap.cxx $(dir_Python)Lib/Measure.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Measure_wrap.cxx -o $(dir_Python)Build/Measure_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblMeasure.so $(dir_Python)Build/Measure_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonNumberCounts: $(dir_Python)/Build $(dir_Python)Build/NumberCounts_wrap.o $(dir_Python)Lib/NumberCounts.i
	$(call insertLine, "from cblNumberCounts import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/NumberCounts_wrap.cxx: $(dir_Python)Lib/NumberCounts.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping NumberCounts with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/NumberCounts_wrap.cxx $(dir_Python)Lib/NumberCounts.i

$(dir_Python)Build/NumberCounts_wrap.o: $(dir_Python)Build/NumberCounts_wrap.cxx $(dir_Python)Lib/NumberCounts.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/NumberCounts_wrap.cxx -o $(dir_Python)Build/NumberCounts_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblNumberCounts.so $(dir_Python)Build/NumberCounts_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonStackedDensityProfile: $(dir_Python)/Build $(dir_Python)Build/StackedDensityProfile_wrap.o $(dir_Python)Lib/StackedDensityProfile.i
	$(call insertLine, "from cblStackedDensityProfile import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/StackedDensityProfile_wrap.cxx: $(dir_Python)Lib/StackedDensityProfile.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping StackedDensityProfile with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/StackedDensityProfile_wrap.cxx $(dir_Python)Lib/StackedDensityProfile.i

$(dir_Python)Build/StackedDensityProfile_wrap.o: $(dir_Python)Build/StackedDensityProfile_wrap.cxx $(dir_Python)Lib/StackedDensityProfile.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/StackedDensityProfile_wrap.cxx -o $(dir_Python)Build/StackedDensityProfile_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblStackedDensityProfile.so $(dir_Python)Build/StackedDensityProfile_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonMassObservableRelation: $(dir_Python)/Build $(dir_Python)Build/MassObservableRelation_wrap.o $(dir_Python)Lib/MassObservableRelation.i
	$(call insertLine, "from cblMassObservableRelation import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/MassObservableRelation_wrap.cxx: $(dir_Python)Lib/MassObservableRelation.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping MassObservableRelation with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/MassObservableRelation_wrap.cxx $(dir_Python)Lib/MassObservableRelation.i

$(dir_Python)Build/MassObservableRelation_wrap.o: $(dir_Python)Build/MassObservableRelation_wrap.cxx $(dir_Python)Lib/MassObservableRelation.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/MassObservableRelation_wrap.cxx -o $(dir_Python)Build/MassObservableRelation_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblMassObservableRelation.so $(dir_Python)Build/MassObservableRelation_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonTwoPointCorrelation: $(dir_Python)/Build $(dir_Python)Build/TwoPointCorrelation_wrap.o $(dir_Python)Lib/TwoPointCorrelation.i
	$(call insertLine, "from cblTwoPointCorrelation import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/TwoPointCorrelation_wrap.cxx: $(dir_Python)Lib/TwoPointCorrelation.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping TwoPointCorrelation with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/TwoPointCorrelation_wrap.cxx $(dir_Python)Lib/TwoPointCorrelation.i

$(dir_Python)Build/TwoPointCorrelation_wrap.o: $(dir_Python)Build/TwoPointCorrelation_wrap.cxx $(dir_Python)Lib/TwoPointCorrelation.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/TwoPointCorrelation_wrap.cxx -o $(dir_Python)Build/TwoPointCorrelation_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblTwoPointCorrelation.so $(dir_Python)Build/TwoPointCorrelation_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonThreePointCorrelation: $(dir_Python)/Build $(dir_Python)Build/ThreePointCorrelation_wrap.o $(dir_Python)Lib/ThreePointCorrelation.i
	$(call insertLine, "from cblThreePointCorrelation import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/ThreePointCorrelation_wrap.cxx: $(dir_Python)Lib/ThreePointCorrelation.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping ThreePointCorrelation with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/ThreePointCorrelation_wrap.cxx $(dir_Python)Lib/ThreePointCorrelation.i

$(dir_Python)Build/ThreePointCorrelation_wrap.o: $(dir_Python)Build/ThreePointCorrelation_wrap.cxx $(dir_Python)Lib/ThreePointCorrelation.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/ThreePointCorrelation_wrap.cxx -o $(dir_Python)Build/ThreePointCorrelation_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblThreePointCorrelation.so $(dir_Python)Build/ThreePointCorrelation_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonPowerSpectrumAngular: $(dir_Python)/Build $(dir_Python)Build/PowerSpectrumAngular_wrap.o $(dir_Python)Lib/PowerSpectrumAngular.i
	$(call insertLine, "from cblPowerSpectrumAngular import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/PowerSpectrumAngular_wrap.cxx: $(dir_Python)Lib/PowerSpectrumAngular.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping PowerSpectrumAngular with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/PowerSpectrumAngular_wrap.cxx $(dir_Python)Lib/PowerSpectrumAngular.i

$(dir_Python)Build/PowerSpectrumAngular_wrap.o: $(dir_Python)Build/PowerSpectrumAngular_wrap.cxx $(dir_Python)Lib/PowerSpectrumAngular.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/PowerSpectrumAngular_wrap.cxx -o $(dir_Python)Build/PowerSpectrumAngular_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblPowerSpectrumAngular.so $(dir_Python)Build/PowerSpectrumAngular_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonModelling: $(dir_Python)/Build $(dir_Python)Build/Modelling_wrap.o $(dir_Python)Lib/Modelling.i
	$(call insertLine, "from cblModelling import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Modelling_wrap.cxx: $(dir_Python)Lib/Modelling.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Modelling with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Modelling_wrap.cxx $(dir_Python)Lib/Modelling.i

$(dir_Python)Build/Modelling_wrap.o: $(dir_Python)Build/Modelling_wrap.cxx $(dir_Python)Lib/Modelling.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Modelling_wrap.cxx -o $(dir_Python)Build/Modelling_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblModelling.so $(dir_Python)Build/Modelling_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonModelling_Cosmology: $(dir_Python)/Build $(dir_Python)Build/Modelling_Cosmology_wrap.o $(dir_Python)Lib/Modelling_Cosmology.i
	$(call insertLine, "from cblModelling_Cosmology import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Modelling_Cosmology_wrap.cxx: $(dir_Python)Lib/Modelling_Cosmology.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Modelling_Cosmology with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Modelling_Cosmology_wrap.cxx $(dir_Python)Lib/Modelling_Cosmology.i

$(dir_Python)Build/Modelling_Cosmology_wrap.o: $(dir_Python)Build/Modelling_Cosmology_wrap.cxx $(dir_Python)Lib/Modelling_Cosmology.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Modelling_Cosmology_wrap.cxx -o $(dir_Python)Build/Modelling_Cosmology_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblModelling_Cosmology.so $(dir_Python)Build/Modelling_Cosmology_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonModelling_DensityProfile: $(dir_Python)/Build $(dir_Python)Build/Modelling_DensityProfile_wrap.o $(dir_Python)Lib/Modelling_DensityProfile.i
	$(call insertLine, "from cblModelling_DensityProfile import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Modelling_DensityProfile_wrap.cxx: $(dir_Python)Lib/Modelling_DensityProfile.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Modelling_DensityProfile with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Modelling_DensityProfile_wrap.cxx $(dir_Python)Lib/Modelling_DensityProfile.i

$(dir_Python)Build/Modelling_DensityProfile_wrap.o: $(dir_Python)Build/Modelling_DensityProfile_wrap.cxx $(dir_Python)Lib/Modelling_DensityProfile.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Modelling_DensityProfile_wrap.cxx -o $(dir_Python)Build/Modelling_DensityProfile_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblModelling_DensityProfile.so $(dir_Python)Build/Modelling_DensityProfile_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonModelling_MassObservableRelation: $(dir_Python)/Build $(dir_Python)Build/Modelling_MassObservableRelation_wrap.o $(dir_Python)Lib/Modelling_MassObservableRelation.i
	$(call insertLine, "from cblModelling_MassObservableRelation import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Modelling_MassObservableRelation_wrap.cxx: $(dir_Python)Lib/Modelling_MassObservableRelation.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Modelling_MassObservableRelation with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Modelling_MassObservableRelation_wrap.cxx $(dir_Python)Lib/Modelling_MassObservableRelation.i

$(dir_Python)Build/Modelling_MassObservableRelation_wrap.o: $(dir_Python)Build/Modelling_MassObservableRelation_wrap.cxx $(dir_Python)Lib/Modelling_MassObservableRelation.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Modelling_MassObservableRelation_wrap.cxx -o $(dir_Python)Build/Modelling_MassObservableRelation_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblModelling_MassObservableRelation.so $(dir_Python)Build/Modelling_MassObservableRelation_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonModelling_NumberCounts: $(dir_Python)/Build $(dir_Python)Build/Modelling_NumberCounts_wrap.o $(dir_Python)Lib/Modelling_NumberCounts.i
	$(call insertLine, "from cblModelling_NumberCounts import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Modelling_NumberCounts_wrap.cxx: $(dir_Python)Lib/Modelling_NumberCounts.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Modelling_NumberCounts with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Modelling_NumberCounts_wrap.cxx $(dir_Python)Lib/Modelling_NumberCounts.i

$(dir_Python)Build/Modelling_NumberCounts_wrap.o: $(dir_Python)Build/Modelling_NumberCounts_wrap.cxx $(dir_Python)Lib/Modelling_NumberCounts.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Modelling_NumberCounts_wrap.cxx -o $(dir_Python)Build/Modelling_NumberCounts_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblModelling_NumberCounts.so $(dir_Python)Build/Modelling_NumberCounts_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonModelling_TwoPointCorrelation: $(dir_Python)/Build $(dir_Python)Build/Modelling_TwoPointCorrelation_wrap.o $(dir_Python)Lib/Modelling_TwoPointCorrelation.i
	$(call insertLine, "from cblModelling_TwoPointCorrelation import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Modelling_TwoPointCorrelation_wrap.cxx: $(dir_Python)Lib/Modelling_TwoPointCorrelation.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Modelling_TwoPointCorrelation with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Modelling_TwoPointCorrelation_wrap.cxx $(dir_Python)Lib/Modelling_TwoPointCorrelation.i

$(dir_Python)Build/Modelling_TwoPointCorrelation_wrap.o: $(dir_Python)Build/Modelling_TwoPointCorrelation_wrap.cxx $(dir_Python)Lib/Modelling_TwoPointCorrelation.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Modelling_TwoPointCorrelation_wrap.cxx -o $(dir_Python)Build/Modelling_TwoPointCorrelation_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblModelling_TwoPointCorrelation.so $(dir_Python)Build/Modelling_TwoPointCorrelation_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonModelling_ThreePointCorrelation: $(dir_Python)/Build $(dir_Python)Build/Modelling_ThreePointCorrelation_wrap.o $(dir_Python)Lib/Modelling_ThreePointCorrelation.i
	$(call insertLine, "from cblModelling_ThreePointCorrelation import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Modelling_ThreePointCorrelation_wrap.cxx: $(dir_Python)Lib/Modelling_ThreePointCorrelation.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Modelling_ThreePointCorrelation with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Modelling_ThreePointCorrelation_wrap.cxx $(dir_Python)Lib/Modelling_ThreePointCorrelation.i

$(dir_Python)Build/Modelling_ThreePointCorrelation_wrap.o: $(dir_Python)Build/Modelling_ThreePointCorrelation_wrap.cxx $(dir_Python)Lib/Modelling_ThreePointCorrelation.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Modelling_ThreePointCorrelation_wrap.cxx -o $(dir_Python)Build/Modelling_ThreePointCorrelation_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblModelling_ThreePointCorrelation.so $(dir_Python)Build/Modelling_ThreePointCorrelation_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonModelling_PowerSpectrumAngular: $(dir_Python)/Build $(dir_Python)Build/Modelling_PowerSpectrumAngular_wrap.o $(dir_Python)Lib/Modelling_PowerSpectrumAngular.i
	$(call insertLine, "from cblModelling_PowerSpectrumAngular import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/Modelling_PowerSpectrumAngular_wrap.cxx: $(dir_Python)Lib/Modelling_PowerSpectrumAngular.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping Modelling_PowerSpectrumAngular with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/Modelling_PowerSpectrumAngular_wrap.cxx $(dir_Python)Lib/Modelling_PowerSpectrumAngular.i

$(dir_Python)Build/Modelling_PowerSpectrumAngular_wrap.o: $(dir_Python)Build/Modelling_PowerSpectrumAngular_wrap.cxx $(dir_Python)Lib/Modelling_PowerSpectrumAngular.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/Modelling_PowerSpectrumAngular_wrap.cxx -o $(dir_Python)Build/Modelling_PowerSpectrumAngular_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblModelling_PowerSpectrumAngular.so $(dir_Python)Build/Modelling_PowerSpectrumAngular_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonGlobalFunc: $(dir_Python)/Build $(dir_Python)Build/GlobalFunc_wrap.o $(dir_Python)Lib/GlobalFunc.i
	$(call insertLine, "from cblGlobalFunc import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/GlobalFunc_wrap.cxx: $(dir_Python)Lib/GlobalFunc.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping GlobalFunc with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/GlobalFunc_wrap.cxx $(dir_Python)Lib/GlobalFunc.i

$(dir_Python)Build/GlobalFunc_wrap.o: $(dir_Python)Build/GlobalFunc_wrap.cxx $(dir_Python)Lib/GlobalFunc.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/GlobalFunc_wrap.cxx -o $(dir_Python)Build/GlobalFunc_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblGlobalFunc.so $(dir_Python)Build/GlobalFunc_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL

###

pythonReadParameters: $(dir_Python)/Build $(dir_Python)Build/ReadParameters_wrap.o $(dir_Python)Lib/ReadParameters.i
	$(call insertLine, "from cblReadParameters import *", $(dir_Python)CosmoBolognaLib/__init__.py)

$(dir_Python)Build/ReadParameters_wrap.cxx: $(dir_Python)Lib/ReadParameters.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Wrapping ReadParameters with swig... "\n")	
	$(SWIG) $(SWIG_FLAG) -I$(dir_H) -outdir $(dir_Python) -o $(dir_Python)Build/ReadParameters_wrap.cxx $(dir_Python)Lib/ReadParameters.i

$(dir_Python)Build/ReadParameters_wrap.o: $(dir_Python)Build/ReadParameters_wrap.cxx $(dir_Python)Lib/ReadParameters.i $(HH) $(PWD)/Makefile
	$(CXX) $(FLAGST) $(SWIG_FLAG_ADD) $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)Build/ReadParameters_wrap.cxx -o $(dir_Python)Build/ReadParameters_wrap.o
	$(CXX) $(FLAGS_LINK) -o $(dir_Python)_cblReadParameters.so $(dir_Python)Build/ReadParameters_wrap.o -Wl,-rpath,$(PWD)/ -L$(PWD)/ -lCBL


####################################################################


$(PWD)/External/Eigen/eigen-$(EIGEN_VERSION)/Eigen/Dense:
	cd $(PWD)/External/Eigen/ && tar -xzf eigen-$(EIGEN_VERSION).tar.gz

$(CUBA_LIB):
	$(CUBA_COMPILE)

$(CCfits_LIB):
	$(CCfits_COMPILE)

$(PWD)/External/CAMB/fortran/camb:
	cd $(PWD)/External/CAMB/forutils ; make clean F90C=$(F)
	cd $(PWD)/External/CAMB/fortran ; make clean F90C=$(F) COMPILER=$(F) && make F90C=$(F) COMPILER=$(F)" -w" && make clean F90C=$(F) COMPILER=$(F) && cd ../../../

$(PWD)/External/CLASS/class:
	cd $(PWD)/External/CLASS ; mkdir -p output && make clean && make CC=$(CC) PYTHON=$(PY) OPTFLAG="-O3 -w" && make clean && cd ../..
        #cd $(PWD)/External/CLASS ; make clean && make CC=$(CC) PYTHON=$(PY) OPTFLAG="-O3 -w" && make clean && cd ../..

$(PWD)/External/MPTbreeze-v1/mptbreeze:
	cd External/MPTbreeze-v1/Cuba-1.4/ ; ./configure CC=$(CC)" -w" F77=$(F)" -w" && make lib && cd ../../../
	cd $(PWD)/External/MPTbreeze-v1 ; ./compile.sh $(F) && cd ../..

$(PWD)/External/mangle/bin/ransack:
	cd $(PWD)/External/mangle/src && mkdir -p ../bin && chmod +x configure && ./configure && make CC=$(CC)" -w" F77=$(F)" -w" && make clean && cd -

$(PWD)/External/VIPERS/venice3.9/venice:
	cd $(PWD)/External/VIPERS/venice3.9 ; make clean && make CC=$(CC) && make && make clean && cd -

$(PWD)/External/CPT_Library/read_pk_sum:
	cd $(PWD)/External/CPT_Library ; make F90C=$(F) switch="-O3 -w" -Wargument-mismatch && cd ../..

$(PWD)/External/CAMB_SPT_private/camb:
	cd $(PWD)/External/CAMB_SPT_private ; make clean && make F90C=$(F) && make clean && cd ../..

$(PWD)/External/MGCAMB/camb:
	cd $(PWD)/External/MGCAMB ; make clean F90C=$(F) COMPILER=$(F) && make F90C=$(F) COMPILER=$(F)" -w" && make clean F90C=$(F) COMPILER=$(F) && cd ../../

$(PWD)/Logo/logo:
	$(CXX) $(PWD)/Logo/Logo.cpp -Wl,-rpath,$(HOME)/lib/ -Wl,-rpath,$(PWD) -L$(PWD) $(FLAGS_INC) -lPATH -lKERNEL -o $(PWD)/Logo/logo ; $(PWD)/Logo/logo

$(dir_Python)/Build:
	mkdir -p $(dir_Python)Build
	mkdir -p $(dir_Python)CosmoBolognaLib
	cd $(dir_Python)CosmoBolognaLib && touch -a __init__.py
