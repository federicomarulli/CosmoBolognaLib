C = g++
F = gfortran
SWIG = swig3.0

FLAGS0 = -std=c++11 -fopenmp

FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wfatal-errors -Werror

FLAGS_FFTLOG = -fPIC -w

Dir_H = Headers/Lib/
Dir_O = Headers/Objects/
Dir_CUBA = External/Cuba-4.2/
Dir_FFTLOG = External/fftlog-f90-master/

dir_H = $(addprefix $(PWD)/,$(Dir_H))
dir_O = $(addprefix $(PWD)/,$(Dir_O))
dir_CUBA = $(addprefix $(PWD)/,$(Dir_CUBA))
dir_FFTLOG = $(addprefix $(PWD)/,$(Dir_FFTLOG))

dir_Python = $(PWD)/Python/

HH = $(dir_H)*.h $(dir_O)*.h

FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dir_CUBA) -I$(dir_H) -I$(dir_O)
FLAGS_FFTW = -lfftw3 #-lfftw3_omp
FLAGS_GSL = -lgsl -lgslcblas -lm -L$(HOME)/lib

CUBA_LIB = $(dir_CUBA)libcuba.so
FLAGS_CUBA = -Wl,-rpath,$(dir_CUBA) -L$(dir_CUBA) -lcuba
CUBA_COMPILE = cd $(dir_CUBA) && ./makeshared.sh 

FLAGS_LINK = -shared

PFLAGS = -I/usr/include/python2.7 

ES = so

Dvar = -DLINUX


FLAGS_PY:=$(shell python -c 'from distutils import sysconfig; print sysconfig.get_config_var("LIBDIR")')

ifeq ($(SYS),MAC)
	Dvar = -DMAC
	FLAGS0 = -std=c++11 -fopenmp
	FLAGS_FFTW = -lfftw3 
	FLAGS_LINK = -dynamiclib -undefined suppress -flat_namespace
        ES = dylib
	FLAGS_PY = -L$(shell python -c 'from distutils import sysconfig; print sysconfig.get_config_var("LIBDIR")') -lpython2.7 -ldl	
	CUBA_LIB = $(dir_CUBA)libcuba.a
	FLAGS_CUBA = $(dir_CUBA)libcuba.a
	CUBA_COMPILE = cd $(dir_CUBA) && ./configure && make lib
endif

FLAGST = $(FLAGS0) $(FLAGS)


####################################################################


##### CBL directories #####

Dir_FUNC = Func/
Dir_STAT = Statistics/
Dir_COSM = Cosmology/Lib/
Dir_CM = ChainMesh/
Dir_CAT = Catalogue/
Dir_LN = LogNormal/
Dir_TWOP = Measure/TwoPointCorrelation/
Dir_THREEP = Measure/ThreePointCorrelation/
Dir_MODEL_GLOB = Modelling/Global/
Dir_MODEL_COSM = Modelling/Cosmology/
Dir_MODEL_TWOP = Modelling/TwoPointCorrelation/
Dir_MODEL_THREEP = Modelling/ThreePointCorrelation/
Dir_GLOB = GlobalFunc/
Dir_READP = ReadParameters/

dir_FUNC = $(addprefix $(PWD)/,$(Dir_FUNC))
dir_STAT = $(addprefix $(PWD)/,$(Dir_STAT))
dir_COSM = $(addprefix $(PWD)/,$(Dir_COSM))
dir_CM = $(addprefix $(PWD)/,$(Dir_CM))
dir_CAT = $(addprefix $(PWD)/,$(Dir_CAT))
dir_LN = $(addprefix $(PWD)/,$(Dir_LN))
dir_TWOP = $(addprefix $(PWD)/,$(Dir_TWOP))
dir_THREEP = $(addprefix $(PWD)/,$(Dir_THREEP))
dir_MODEL_GLOB = $(addprefix $(PWD)/,$(Dir_MODEL_GLOB))
dir_MODEL_COSM = $(addprefix $(PWD)/,$(Dir_MODEL_COSM))
dir_MODEL_TWOP = $(addprefix $(PWD)/,$(Dir_MODEL_TWOP))
dir_MODEL_THREEP = $(addprefix $(PWD)/,$(Dir_MODEL_THREEP))
dir_GLOB = $(addprefix $(PWD)/,$(Dir_GLOB))
dir_READP = $(addprefix $(PWD)/,$(Dir_READP))

##### FFTlog object files #####

OBJ_FFTLOG = $(dir_FFTLOG)drffti.o $(dir_FFTLOG)drfftb.o $(dir_FFTLOG)drfftf.o $(dir_FFTLOG)fftlog.o $(dir_FFTLOG)cdgamma.o

##### CBL object files #####

OBJ_FUNC = $(dir_FUNC)Func.o $(dir_FUNC)FuncXi.o $(dir_FUNC)FuncMultipoles.o $(dir_FUNC)GSLfunction.o  $(dir_FUNC)Data.o $(dir_FUNC)Data1D.o $(dir_FUNC)Data1D_collection.o $(dir_FUNC)Data2D.o $(dir_FUNC)Data1D_extra.o $(dir_FUNC)Data2D_extra.o $(dir_FUNC)Field3D.o $(dir_FUNC)FuncGrid.o $(dir_FUNC)GSLwrapper.o $(dir_FUNC)CUBAwrapper.o $(dir_FUNC)RandomNumbers.o $(dir_FUNC)Distribution.o $(OBJ_FFTLOG) $(dir_FUNC)FFTlog.o

OBJ_STAT = $(dir_STAT)Model.o $(dir_STAT)Model1D.o $(dir_STAT)Model2D.o $(dir_STAT)Chain.o $(dir_STAT)Parameter.o $(dir_STAT)BaseParameter.o $(dir_STAT)DerivedParameter.o $(dir_STAT)LikelihoodParameters.o $(dir_STAT)LikelihoodFunction.o $(dir_STAT)Sampler.o  $(dir_STAT)Likelihood.o

OBJ_COSM = $(dir_COSM)Cosmology.o $(dir_COSM)Sigma.o $(dir_COSM)PkXi.o $(dir_COSM)PkXizSpace.o $(dir_COSM)MassFunction.o $(dir_COSM)Bias.o $(dir_COSM)RSD.o $(dir_COSM)DensityProfile.o $(dir_COSM)Velocities.o $(dir_COSM)MassGrowth.o $(dir_COSM)NG.o $(dir_COSM)BAO.o $(dir_COSM)SizeFunction.o  $(dir_COSM)3PCF.o

OBJ_CM = $(dir_CM)ChainMesh.o

OBJ_CAT = $(dir_CAT)Object.o $(dir_CAT)Catalogue.o $(dir_CAT)RandomCatalogue.o $(dir_CAT)ChainMesh_Catalogue.o $(dir_CAT)RandomCatalogueVIPERS.o $(dir_CAT)VoidCatalogue.o $(dir_CAT)GadgetCatalogue.o

OBJ_LN = $(dir_LN)LogNormal.o $(dir_LN)LogNormalFull.o

OBJ_TWOP = $(dir_TWOP)Pair.o $(dir_TWOP)Pair1D.o $(dir_TWOP)Pair2D.o $(dir_TWOP)Pair1D_extra.o $(dir_TWOP)Pair2D_extra.o $(dir_TWOP)TwoPointCorrelation.o $(dir_TWOP)TwoPointCorrelation1D.o $(dir_TWOP)TwoPointCorrelation1D_angular.o $(dir_TWOP)TwoPointCorrelation1D_monopole.o $(dir_TWOP)TwoPointCorrelation2D.o $(dir_TWOP)TwoPointCorrelation2D_cartesian.o $(dir_TWOP)TwoPointCorrelation2D_polar.o $(dir_TWOP)TwoPointCorrelation_projected.o $(dir_TWOP)TwoPointCorrelation_deprojected.o $(dir_TWOP)TwoPointCorrelation_multipoles.o $(dir_TWOP)TwoPointCorrelation_wedges.o $(dir_TWOP)TwoPointCorrelation1D_filtered.o $(dir_TWOP)TwoPointCorrelationCross.o $(dir_TWOP)TwoPointCorrelationCross1D.o $(dir_TWOP)TwoPointCorrelationCross1D_angular.o $(dir_TWOP)TwoPointCorrelationCross1D_monopole.o

OBJ_THREEP = $(dir_THREEP)Triplet.o $(dir_THREEP)ThreePointCorrelation.o $(dir_THREEP)ThreePointCorrelation_angular_connected.o $(dir_THREEP)ThreePointCorrelation_angular_reduced.o $(dir_THREEP)ThreePointCorrelation_comoving_connected.o $(dir_THREEP)ThreePointCorrelation_comoving_reduced.o 

OBJ_MODEL_GLOB = $(dir_MODEL_GLOB)Modelling.o

OBJ_MODEL_COSM = $(dir_MODEL_COSM)ModelFunction_Cosmology.o $(dir_MODEL_COSM)Modelling_Cosmology.o

OBJ_MODEL_TWOP = $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.o 

OBJ_MODEL_THREEP = $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.o 

OBJ_GLOB = $(dir_GLOB)FuncCosmology.o $(dir_GLOB)Func.o $(dir_GLOB)SubSample.o $(dir_GLOB)Reconstruction.o $(dir_GLOB)Forecast.o

OBJ_READP = $(dir_READP)ReadParameters.o


OBJ_CBL = $(OBJ_FUNC) $(OBJ_STAT) $(OBJ_COSM) $(OBJ_CM) $(OBJ_CAT) $(OBJ_LN) $(OBJ_TWOP) $(OBJ_THREEP) $(OBJ_MODEL_GLOB) $(OBJ_MODEL_COSM) $(OBJ_MODEL_TWOP) $(OBJ_MODEL_THREEP) $(OBJ_GLOB) $(OBJ_READP)

OBJ_ALL = $(OBJ_CBL) $(dir_FUNC)conv.o $(PWD)/External/CAMB/*.o $(PWD)/External/classgal_v1/*.o $(PWD)/External/mangle/*.o $(PWD)/External/MPTbreeze-v1/*.o 


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
	$(call colorecho, "\n"Compiling the library: libFUNC... "\n")
	make -j3 libFUNC
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
	$(call colorecho, "\n"Compiling the library: libTWOP... "\n")
	make -j3 libTWOP
	$(call colorecho, "\n"Compiling the library: libTHREEP... "\n")
	make -j3 libTHREEP
	$(call colorecho, "\n"Compiling the library: libMODEL_GLOB... "\n")
	make -j3 libMODEL_GLOB
	$(call colorecho, "\n"Compiling the library: libMODEL_COSM... "\n")
	make -j3 libMODEL_COSM
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

libFUNC: $(OBJ_FUNC) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libFUNC.$(ES) $(OBJ_FUNC) $(FLAGS_CUBA) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -lgfortran

libSTAT: $(OBJ_STAT) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libSTAT.$(ES) $(OBJ_STAT) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC

libCOSM: $(OBJ_COSM) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libCOSM.$(ES) $(OBJ_COSM) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT

libCM: $(OBJ_CM) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libCM.$(ES) $(OBJ_CM) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM

libCAT: $(OBJ_CAT) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libCAT.$(ES) $(OBJ_CAT) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM

libLN: $(OBJ_LN) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libLN.$(ES) $(OBJ_LN) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT 

libTWOP: $(OBJ_TWOP) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libTWOP.$(ES) $(OBJ_TWOP) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN

libTHREEP: $(OBJ_THREEP) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libTHREEP.$(ES) $(OBJ_THREEP) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP

libMODEL_GLOB: $(OBJ_MODEL_GLOB) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libMODEL_GLOB.$(ES) $(OBJ_MODEL_GLOB) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP

libMODEL_COSM: $(OBJ_MODEL_COSM) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libMODEL_COSM.$(ES) $(OBJ_MODEL_COSM) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP -lMODEL_GLOB

libMODEL_TWOP: $(OBJ_MODEL_TWOP) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libMODEL_TWOP.$(ES) $(OBJ_MODEL_TWOP) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP -lMODEL_GLOB

libMODEL_THREEP: $(OBJ_MODEL_THREEP) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libMODEL_THREEP.$(ES) $(OBJ_MODEL_THREEP) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP -lMODEL_GLOB

libGLOB: $(OBJ_GLOB) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libGLOB.$(ES) $(OBJ_GLOB) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP -lMODEL_GLOB -lMODEL_COSM -lMODEL_TWOP -lMODEL_THREEP

libREADP: $(OBJ_READP) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libREADP.$(ES) $(OBJ_READP) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC

libCBL: $(OBJ_CBL) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libCBL.$(ES) $(OBJ_CBL) $(FLAGS_CUBA) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW) -lgfortran

conv: $(dir_FUNC)conv.o
	$(F) -o $(dir_FUNC)conv $(dir_FUNC)conv.o 

CUBA: $(CUBA_LIB)

CAMB: $(PWD)/External/CAMB/camb

CLASS: $(PWD)/External/classgal_v1/class

MPTbreeze:  $(PWD)/External/MPTbreeze-v1/mptbreeze

fftlog-f90: $(PWD)/External/fftlog-f90-master/fftlog-f90

mangle: $(PWD)/External/mangle/bin/ransack

venice: $(PWD)/External/VIPERS/venice3.9/venice

allExamples:
	$(call colorecho, "\n"Compiling the example code: vector.cpp ... "\n")
	cd $(PWD)/Examples/vectors ; make
	$(call colorecho, "\n"Compiling the example code: randomNumbers.cpp ... "\n")
	cd $(PWD)/Examples/randomNumbers ; make 
	$(call colorecho, "\n"Compiling the example code: integration_gsl.cpp ... "\n")
	cd $(PWD)/Examples/wrappers ; make integration_gsl
	$(call colorecho, "\n"Compiling the example code: integration_cuba.cpp ... "\n")
	cd $(PWD)/Examples/wrappers ; make integration_cuba
	$(call colorecho, "\n"Compiling the example code: minimisation.cpp ... "\n")
	cd $(PWD)/Examples/wrappers ; make minimisation
	$(call colorecho, "\n"Compiling the example code: distances.cpp ... "\n")
	cd $(PWD)/Examples/distances ; make 
	$(call colorecho, "\n"Compiling the example code: covsample.cpp ... "\n")
	cd $(PWD)/Examples/covsample ; make 
	$(call colorecho, "\n"Compiling the example code: cosmology.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make cosmology
	$(call colorecho, "\n"Compiling the example code: fsigma8.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make fsigma8
	$(call colorecho, "\n"Compiling the example code: model_cosmology.cpp ... "\n")
	cd $(PWD)/Examples/cosmology ; make model_cosmology
	$(call colorecho, "\n"Compiling the example code: prior.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make prior 
	$(call colorecho, "\n"Compiling the example code: fit.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make fit
	$(call colorecho, "\n"Compiling the example code: sampler.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make sampler
	$(call colorecho, "\n"Compiling the example code: catalogue.cpp ... "\n")
	cd $(PWD)/Examples/catalogue ; make catalogue 
	$(call colorecho, "\n"Compiling the example code: 2pt_monopole.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_monopole
	$(call colorecho, "\n"Compiling the example code: 2pt_monopole_errors.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_monopole_errors
	$(call colorecho, "\n"Compiling the example code: 2pt_2D.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_2D
	$(call colorecho, "\n"Compiling the example code: 2pt_projected.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_projected
	$(call colorecho, "\n"Compiling the example code: 2pt_angular.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_angular
	$(call colorecho, "\n"Compiling the example code: 3pt.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 3pt
	$(call colorecho, "\n"Compiling the example code: model_2pt_monopole_BAO.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_monopole_BAO
	$(call colorecho, "\n"Compiling the example code: model_2pt_monopole_RSD.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_monopole_RSD
	$(call colorecho, "\n"Compiling the example code: model_2pt_projected.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_projected
	$(call colorecho, "\n"Compiling the example code: model_2pt_2D.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_2pt_2D
	$(call colorecho, "\n"Compiling the example code: model_3pt.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make model_3pt
	$(call colorecho, "\n"Compiling the example code: sizeFunction.cpp ... "\n")
	cd $(PWD)/Examples/cosmicVoids/codes ; make sizeFunction
	$(call colorecho, "\n"Compiling the example code: cleanVoidCatalogue.cpp ... "\n")
	cd $(PWD)/Examples/cosmicVoids/codes ; make cleanVoidCatalogue
	$(call colorecho, "\n"Compiling the example code: readParameterFile.cpp ... "\n")
	cd $(PWD)/Examples/readParameterFile/ ; make 

swig_wrapper: $(dir_Python)CBL_wrap.cxx

python: $(dir_Python)CBL_wrap.o $(OBJ_CBL) $(dir_Python)CBL.i 
	$(C) -shared $(OBJ_CBL) $(dir_Python)CBL_wrap.o -o $(dir_Python)/_CosmoBolognaLib.so $(FLAGS_CUBA) $(FLAGS_GSL) $(FLAGS_FFTW) -lgomp $(FLAGS_PY) -lgfortran
	mv $(dir_Python)/_CosmoBolognaLib.so $(dir_Python)/CosmoBolognaLib/
	mv $(dir_Python)/CosmoBolognaLib.py $(dir_Python)/CosmoBolognaLib/

doc:
	rm Doc/html/* Doc/xml/* -rf
	doxygen Doc/dconfig
	rm Doc/doxygen_sqlite3.db -f 
#	python ../bin/doxy2swig2.py Doc/xml/index.xml Doc/documentation.i
#	python ../doxy2swig/doxy2swig.py Doc/xml/index.xml Doc/documentation.i

doct:
	rm Doc/html/* Doc/xml/* -rf
	doxygen Doc/dconfigT
	rm Doc/doxygen_sqlite3.db -f 

cleanExamples:
	cd $(PWD)/Examples/vectors ; make clean && cd ../..
	cd $(PWD)/Examples/randomNumbers ; make clean && cd ../..
	cd $(PWD)/Examples/wrappers ; make clean && cd ../..
	cd $(PWD)/Examples/distances ; make clean && cd ../..
	cd $(PWD)/Examples/covsample ; make clean && cd ../..
	cd $(PWD)/Examples/cosmology ; make clean && cd ../..
	cd $(PWD)/Examples/statistics/codes ; make clean && cd ../..
	cd $(PWD)/Examples/catalogue ; make clean && cd ../..
	cd $(PWD)/Examples/clustering/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/cosmicVoids/codes ; make clean && cd ../../..
	cd $(PWD)/Examples/readParameterFile ; make clean && cd ../..
	rm -rf $(PWD)/Examples/statistics/output/* $(PWD)/Examples/clustering/output/* $(PWD)/Examples/cosmicVoids/output/*


cleanpy:
	rm -f $(dir_Python)*~ $(dir_Python)CBL_wrap.o $(dir_Python)CBL_wrap.cxx $(dir_Python)CosmoBolognaLib.py*
	rm -rf $(dir_Python)dist $(dir_Python)build $(dir_Python)CosmoBolognaLib.egg-info
	rm -f $(dir_Python)Lib/*~ $(dir_Python)Lib/*.o $(dir_Python)Lib/*.cxx $(dir_Python)Lib/*.py
	rm -f $(dir_Python)/CosmoBolognaLib/*CosmoBolognaLib* $(dir_Python)CosmoBolognaLib/*~ $(dir_Python)CosmoBolognaLib/*.pyc

purgepy:
	make cleanpy

cleanTEMP:
	rm -f $(OBJ_ALL) core* $(PWD)/*~ $(dir_FUNC)*~ $(dir_STAT)*~ $(dir_COSM)*~ $(dir_CM)*~ $(dir_CAT)*~ $(dir_LN)*~ $(dir_TWOP)*~ $(dir_MODEL_GLOB)*~ $(dir_MODEL_COSM)*~ $(dir_MODEL_TWOP)*~ $(dir_MODEL_THREEP)*~ $(dir_THREEP)*~ $(dir_GLOB)*~ $(dir_READP)*~ $(dir_H)*~ $(dir_O)*~ $(PWD)/\#* $(dir_FUNC)\#* $(dir_STAT)\#* $(dir_COSM)\#* $(dir_CM)\#* $(dir_CAT)\#* $(dir_LN)\#* $(dir_TWOP)\#* $(dir_THREEP)\#* $(dir_MODEL_GLOB)\#* $(dir_MODEL_COSM)\#* $(dir_MODEL_TWOP)\#* $(dir_MODEL_THREEP)\#* $(dir_GLOB)\#* $(dir_READP)\#* $(dir_H)\#* $(dir_O)\#* $(PWD)/Doc/WARNING_LOGFILE* $(PWD)/Doc/*~

clean:
	make cleanExamples
	make cleanTEMP

purge:
	make clean
	rm -f *.$(ES) temp*

purgeALL:
	make purge
	make cleanpy
	rm -rf Cosmology/Tables/* ;
	cd External/CAMB ; make clean ; cd .. ;
	rm -rf External/CAMB/camb
	rm -rf External/CAMB/output_linear/* ;
	rm -rf External/CAMB/output_nonlinear/* ;
	rm -rf External/CAMB/test_* ; 
	rm -rf External/VIPERS/venice3.9/venice
	rm -rf External/mangle/bin/*
	cd External/mangle/src; make cleaner ; cd .. ;
	cd External/classgal_v1/ ; make clean ; cd .. ;
	rm -rf External/classgal_v1/output_linear/* ;
	rm -rf External/classgal_v1/output_nonlinear/* ;
	cd External/fftlog-f90-master/ ; make clean ; cd .. ;
	cd External/mangle/src ; make clean ; cd .. ;
	cd External/MPTbreeze-v1/Cuba-1.4/ ; make distclean ;
	rm -rf External/MPTbreeze-v1/mptbreeze ;
	rm -rf External/MPTbreeze-v1/*~ ;
	rm -rf External/MPTbreeze-v1/output_linear/* ;
	rm -rf External/MPTbreeze-v1/output_nonlinear/* ;
	cd External/Cuba-4.2 ; rm -rf config.h config.log config.status demo-fortran.dSYM/ libcuba.a libcuba.so makefile *~ ; cd .. ;


#################################################################### 


$(CUBA_LIB):
	$(CUBA_COMPILE)

$(dir_FUNC)Func.o: $(dir_FUNC)Func.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Func.cpp -o $(dir_FUNC)Func.o 

$(dir_FUNC)FuncXi.o: $(dir_FUNC)FuncXi.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)FuncXi.cpp -o $(dir_FUNC)FuncXi.o

$(dir_FUNC)FuncMultipoles.o: $(dir_FUNC)FuncMultipoles.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)FuncMultipoles.cpp -o $(dir_FUNC)FuncMultipoles.o

$(dir_FUNC)GSLfunction.o: $(dir_FUNC)GSLfunction.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)GSLfunction.cpp -o $(dir_FUNC)GSLfunction.o

$(dir_FUNC)Data.o: $(dir_FUNC)Data.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Data.cpp -o $(dir_FUNC)Data.o

$(dir_FUNC)Data1D.o: $(dir_FUNC)Data1D.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Data1D.cpp -o $(dir_FUNC)Data1D.o

$(dir_FUNC)Data1D_collection.o: $(dir_FUNC)Data1D_collection.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Data1D_collection.cpp -o $(dir_FUNC)Data1D_collection.o

$(dir_FUNC)Data2D.o: $(dir_FUNC)Data2D.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Data2D.cpp -o $(dir_FUNC)Data2D.o

$(dir_FUNC)Data1D_extra.o: $(dir_FUNC)Data1D_extra.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Data1D_extra.cpp -o $(dir_FUNC)Data1D_extra.o

$(dir_FUNC)Data2D_extra.o: $(dir_FUNC)Data2D_extra.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Data2D_extra.cpp -o $(dir_FUNC)Data2D_extra.o

$(dir_FUNC)Field3D.o: $(dir_FUNC)Field3D.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Field3D.cpp -o $(dir_FUNC)Field3D.o

$(dir_FUNC)FuncGrid.o: $(dir_FUNC)FuncGrid.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)FuncGrid.cpp -o $(dir_FUNC)FuncGrid.o

$(dir_FUNC)GSLwrapper.o: $(dir_FUNC)GSLwrapper.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_FUNC)GSLwrapper.cpp -o $(dir_FUNC)GSLwrapper.o 

$(dir_FUNC)CUBAwrapper.o: $(dir_FUNC)CUBAwrapper.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_FUNC)CUBAwrapper.cpp -o $(dir_FUNC)CUBAwrapper.o

$(dir_FUNC)RandomNumbers.o: $(dir_FUNC)RandomNumbers.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_FUNC)RandomNumbers.cpp -o $(dir_FUNC)RandomNumbers.o 

$(dir_FUNC)Distribution.o: $(dir_FUNC)Distribution.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Distribution.cpp -o $(dir_FUNC)Distribution.o 

$(dir_FUNC)conv.o: $(dir_FUNC)conv.f90 
	$(F) -c $(dir_FUNC)conv.f90 -o $(dir_FUNC)conv.o

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

$(dir_FUNC)FFTlog.o:  $(OBJ_FFTLOG) $(dir_FUNC)FFTlog.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) $(Dvar) -c -fPIC $(FLAGS_INC) $(dir_FUNC)FFTlog.cpp -o $(dir_FUNC)FFTlog.o 


#################################################################### 


$(dir_STAT)Model.o: $(dir_STAT)Model.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Model.cpp -o $(dir_STAT)Model.o

$(dir_STAT)Model1D.o: $(dir_STAT)Model1D.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Model1D.cpp -o $(dir_STAT)Model1D.o

$(dir_STAT)Model2D.o: $(dir_STAT)Model2D.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Model2D.cpp -o $(dir_STAT)Model2D.o

$(dir_STAT)Chain.o: $(dir_STAT)Chain.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Chain.cpp -o $(dir_STAT)Chain.o

$(dir_STAT)Parameter.o: $(dir_STAT)Parameter.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Parameter.cpp -o $(dir_STAT)Parameter.o

$(dir_STAT)BaseParameter.o: $(dir_STAT)BaseParameter.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)BaseParameter.cpp -o $(dir_STAT)BaseParameter.o

$(dir_STAT)DerivedParameter.o: $(dir_STAT)DerivedParameter.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)DerivedParameter.cpp -o $(dir_STAT)DerivedParameter.o

$(dir_STAT)LikelihoodParameters.o: $(dir_STAT)LikelihoodParameters.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)LikelihoodParameters.cpp -o $(dir_STAT)LikelihoodParameters.o

$(dir_STAT)LikelihoodFunction.o: $(dir_STAT)LikelihoodFunction.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)LikelihoodFunction.cpp -o $(dir_STAT)LikelihoodFunction.o

$(dir_STAT)Sampler.o: $(dir_STAT)Sampler.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Sampler.cpp -o $(dir_STAT)Sampler.o

$(dir_STAT)Likelihood.o: $(dir_STAT)Likelihood.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Likelihood.cpp -o $(dir_STAT)Likelihood.o


#################################################################### 


$(dir_COSM)Cosmology.o: $(dir_COSM)Cosmology.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Cosmology.cpp -o $(dir_COSM)Cosmology.o

$(dir_COSM)MassFunction.o: $(dir_COSM)MassFunction.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)MassFunction.cpp -o $(dir_COSM)MassFunction.o

$(dir_COSM)SizeFunction.o: $(dir_COSM)SizeFunction.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)SizeFunction.cpp -o $(dir_COSM)SizeFunction.o

$(dir_COSM)PkXi.o: $(dir_COSM)PkXi.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)PkXi.cpp -o $(dir_COSM)PkXi.o

$(dir_COSM)PkXizSpace.o: $(dir_COSM)PkXizSpace.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)PkXizSpace.cpp -o $(dir_COSM)PkXizSpace.o

$(dir_COSM)Bias.o: $(dir_COSM)Bias.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Bias.cpp -o $(dir_COSM)Bias.o

$(dir_COSM)RSD.o: $(dir_COSM)RSD.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)RSD.cpp -o $(dir_COSM)RSD.o

$(dir_COSM)DensityProfile.o: $(dir_COSM)DensityProfile.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)DensityProfile.cpp -o $(dir_COSM)DensityProfile.o

$(dir_COSM)Sigma.o: $(dir_COSM)Sigma.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Sigma.cpp -o $(dir_COSM)Sigma.o

$(dir_COSM)Velocities.o: $(dir_COSM)Velocities.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)Velocities.cpp -o $(dir_COSM)Velocities.o

$(dir_COSM)MassGrowth.o: $(dir_COSM)MassGrowth.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)MassGrowth.cpp -o $(dir_COSM)MassGrowth.o

$(dir_COSM)NG.o: $(dir_COSM)NG.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)NG.cpp -o $(dir_COSM)NG.o

$(dir_COSM)BAO.o: $(dir_COSM)BAO.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)BAO.cpp -o $(dir_COSM)BAO.o

$(dir_COSM)3PCF.o: $(dir_COSM)3PCF.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_COSM)3PCF.cpp -o $(dir_COSM)3PCF.o

#################################################################### 


$(dir_CM)ChainMesh.o: $(dir_CM)ChainMesh.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CM)ChainMesh.cpp -o $(dir_CM)ChainMesh.o


#################################################################### 


$(dir_CAT)Object.o: $(dir_CAT)Object.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)Object.cpp -o $(dir_CAT)Object.o

$(dir_CAT)Catalogue.o: $(dir_CAT)Catalogue.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)Catalogue.cpp -o $(dir_CAT)Catalogue.o

$(dir_CAT)RandomCatalogue.o: $(dir_CAT)RandomCatalogue.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)RandomCatalogue.cpp -o $(dir_CAT)RandomCatalogue.o

$(dir_CAT)ChainMesh_Catalogue.o: $(dir_CAT)ChainMesh_Catalogue.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)ChainMesh_Catalogue.cpp -o $(dir_CAT)ChainMesh_Catalogue.o

$(dir_CAT)RandomCatalogueVIPERS.o: $(dir_CAT)RandomCatalogueVIPERS.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)RandomCatalogueVIPERS.cpp -o $(dir_CAT)RandomCatalogueVIPERS.o

$(dir_CAT)VoidCatalogue.o: $(dir_CAT)VoidCatalogue.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)VoidCatalogue.cpp -o $(dir_CAT)VoidCatalogue.o

$(dir_CAT)GadgetCatalogue.o: $(dir_CAT)GadgetCatalogue.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_CAT)GadgetCatalogue.cpp -o $(dir_CAT)GadgetCatalogue.o


#################################################################### 


$(dir_LN)LogNormal.o: $(dir_LN)LogNormal.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_LN)LogNormal.cpp -o $(dir_LN)LogNormal.o

$(dir_LN)LogNormalFull.o: $(dir_LN)LogNormalFull.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_LN)LogNormalFull.cpp -o $(dir_LN)LogNormalFull.o


#################################################################### 


$(dir_TWOP)Pair.o: $(dir_TWOP)Pair.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair.cpp -o $(dir_TWOP)Pair.o

$(dir_TWOP)Pair1D.o: $(dir_TWOP)Pair1D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair1D.cpp -o $(dir_TWOP)Pair1D.o

$(dir_TWOP)Pair2D.o: $(dir_TWOP)Pair2D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair2D.cpp -o $(dir_TWOP)Pair2D.o

$(dir_TWOP)Pair1D_extra.o: $(dir_TWOP)Pair1D_extra.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair1D_extra.cpp -o $(dir_TWOP)Pair1D_extra.o

$(dir_TWOP)Pair2D_extra.o: $(dir_TWOP)Pair2D_extra.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair2D_extra.cpp -o $(dir_TWOP)Pair2D_extra.o

$(dir_TWOP)TwoPointCorrelation.o: $(dir_TWOP)TwoPointCorrelation.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation.cpp -o $(dir_TWOP)TwoPointCorrelation.o

$(dir_TWOP)TwoPointCorrelation1D.o: $(dir_TWOP)TwoPointCorrelation1D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation1D.cpp -o $(dir_TWOP)TwoPointCorrelation1D.o

$(dir_TWOP)TwoPointCorrelation2D.o: $(dir_TWOP)TwoPointCorrelation2D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation2D.cpp -o $(dir_TWOP)TwoPointCorrelation2D.o

$(dir_TWOP)TwoPointCorrelation1D_monopole.o: $(dir_TWOP)TwoPointCorrelation1D_monopole.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation1D_monopole.cpp -o $(dir_TWOP)TwoPointCorrelation1D_monopole.o

$(dir_TWOP)TwoPointCorrelation1D_angular.o: $(dir_TWOP)TwoPointCorrelation1D_angular.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation1D_angular.cpp -o $(dir_TWOP)TwoPointCorrelation1D_angular.o

$(dir_TWOP)TwoPointCorrelation2D_cartesian.o: $(dir_TWOP)TwoPointCorrelation2D_cartesian.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation2D_cartesian.cpp -o $(dir_TWOP)TwoPointCorrelation2D_cartesian.o

$(dir_TWOP)TwoPointCorrelation2D_polar.o: $(dir_TWOP)TwoPointCorrelation2D_polar.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation2D_polar.cpp -o $(dir_TWOP)TwoPointCorrelation2D_polar.o

$(dir_TWOP)TwoPointCorrelation_projected.o: $(dir_TWOP)TwoPointCorrelation_projected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_projected.cpp -o $(dir_TWOP)TwoPointCorrelation_projected.o

$(dir_TWOP)TwoPointCorrelation_deprojected.o: $(dir_TWOP)TwoPointCorrelation_deprojected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_deprojected.cpp -o $(dir_TWOP)TwoPointCorrelation_deprojected.o

$(dir_TWOP)TwoPointCorrelation_multipoles.o: $(dir_TWOP)TwoPointCorrelation_multipoles.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_multipoles.cpp -o $(dir_TWOP)TwoPointCorrelation_multipoles.o

$(dir_TWOP)TwoPointCorrelation_wedges.o: $(dir_TWOP)TwoPointCorrelation_wedges.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation_wedges.cpp -o $(dir_TWOP)TwoPointCorrelation_wedges.o

$(dir_TWOP)TwoPointCorrelation1D_filtered.o: $(dir_TWOP)TwoPointCorrelation1D_filtered.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation1D_filtered.cpp -o $(dir_TWOP)TwoPointCorrelation1D_filtered.o

$(dir_TWOP)TwoPointCorrelationCross.o: $(dir_TWOP)TwoPointCorrelationCross.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelationCross.cpp -o $(dir_TWOP)TwoPointCorrelationCross.o

$(dir_TWOP)TwoPointCorrelationCross1D.o: $(dir_TWOP)TwoPointCorrelationCross1D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelationCross1D.cpp -o $(dir_TWOP)TwoPointCorrelationCross1D.o

$(dir_TWOP)TwoPointCorrelaation2D.o: $(dir_TWOP)TwoPointCorrelation2D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelation2D.cpp -o $(dir_TWOP)TwoPointCorrelation2D.o

$(dir_TWOP)TwoPointCorrelationCross1D_monopole.o: $(dir_TWOP)TwoPointCorrelationCross1D_monopole.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelationCross1D_monopole.cpp -o $(dir_TWOP)TwoPointCorrelationCross1D_monopole.o

$(dir_TWOP)TwoPointCorrelationCross1D_angular.o: $(dir_TWOP)TwoPointCorrelationCross1D_angular.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)TwoPointCorrelationCross1D_angular.cpp -o $(dir_TWOP)TwoPointCorrelationCross1D_angular.o


#################################################################### 


$(dir_THREEP)Triplet.o: $(dir_THREEP)Triplet.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)Triplet.cpp -o $(dir_THREEP)Triplet.o

$(dir_THREEP)ThreePointCorrelation.o: $(dir_THREEP)ThreePointCorrelation.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation.cpp -o $(dir_THREEP)ThreePointCorrelation.o

$(dir_THREEP)ThreePointCorrelation_angular_connected.o: $(dir_THREEP)ThreePointCorrelation_angular_connected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_angular_connected.cpp -o $(dir_THREEP)ThreePointCorrelation_angular_connected.o

$(dir_THREEP)ThreePointCorrelation_angular_reduced.o: $(dir_THREEP)ThreePointCorrelation_angular_reduced.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_angular_reduced.cpp -o $(dir_THREEP)ThreePointCorrelation_angular_reduced.o

$(dir_THREEP)ThreePointCorrelation_comoving_connected.o: $(dir_THREEP)ThreePointCorrelation_comoving_connected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_comoving_connected.cpp -o $(dir_THREEP)ThreePointCorrelation_comoving_connected.o

$(dir_THREEP)ThreePointCorrelation_comoving_reduced.o: $(dir_THREEP)ThreePointCorrelation_comoving_reduced.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_THREEP)ThreePointCorrelation_comoving_reduced.cpp -o $(dir_THREEP)ThreePointCorrelation_comoving_reduced.o


#################################################################### 


$(dir_MODEL_GLOB)Modelling.o: $(dir_MODEL_GLOB)Modelling.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_GLOB)Modelling.cpp -o $(dir_MODEL_GLOB)Modelling.o


#################################################################### 


$(dir_MODEL_COSM)Modelling_Cosmology.o: $(dir_MODEL_COSM)Modelling_Cosmology.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_COSM)Modelling_Cosmology.cpp -o $(dir_MODEL_COSM)Modelling_Cosmology.o

$(dir_MODEL_COSM)ModelFunction_Cosmology.o: $(dir_MODEL_COSM)ModelFunction_Cosmology.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_COSM)ModelFunction_Cosmology.cpp -o $(dir_MODEL_COSM)ModelFunction_Cosmology.o


#################################################################### 


$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_monopole.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_monopole.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_angular.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_angular.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_cartesian.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_cartesian.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation2D_polar.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation2D_polar.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_projected.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_projected.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_deprojected.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_deprojected.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_multipoles.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_multipoles.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation_wedges.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation_wedges.o

$(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.o: $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.cpp -o $(dir_MODEL_TWOP)Modelling_TwoPointCorrelation1D_filtered.o

$(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.o: $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.cpp -o $(dir_MODEL_TWOP)ModelFunction_TwoPointCorrelation1D_filtered.o


#################################################################### 


$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation.o

$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_connected.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_connected.o

$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_angular_reduced.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_angular_reduced.o

$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_connected.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_connected.o

$(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.o: $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.cpp -o $(dir_MODEL_THREEP)Modelling_ThreePointCorrelation_comoving_reduced.o

$(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.o: $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.cpp -o $(dir_MODEL_THREEP)ModelFunction_ThreePointCorrelation_comoving_reduced.o


#################################################################### 


$(dir_GLOB)FuncCosmology.o: $(dir_GLOB)FuncCosmology.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)FuncCosmology.cpp -o $(dir_GLOB)FuncCosmology.o

$(dir_GLOB)Func.o: $(dir_GLOB)Func.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)Func.cpp -o $(dir_GLOB)Func.o

$(dir_GLOB)SubSample.o: $(dir_GLOB)SubSample.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)SubSample.cpp -o $(dir_GLOB)SubSample.o

$(dir_GLOB)Reconstruction.o: $(dir_GLOB)Reconstruction.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)Reconstruction.cpp -o $(dir_GLOB)Reconstruction.o

$(dir_GLOB)Forecast.o: $(dir_GLOB)Forecast.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)Forecast.cpp -o $(dir_GLOB)Forecast.o

#################################################################### 


$(dir_READP)ReadParameters.o: $(dir_READP)ReadParameters.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_READP)ReadParameters.cpp -o $(dir_READP)ReadParameters.o


#################################################################### 


$(dir_Python)CBL_wrap.o: $(dir_Python)CBL_wrap.cxx $(dir_Python)CBL.i $(HH) 
	$(call colorecho, "\n"Compiling the python wrapper. It may take a few minutes ... "\n")
	$(C) $(FLAGST) -Wno-uninitialized $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)CBL_wrap.cxx -o $(dir_Python)CBL_wrap.o

$(dir_Python)CBL_wrap.cxx: $(dir_Python)CBL.i $(HH) 
	$(SWIG) -python -c++ -I$(dir_H) -I$(dir_O) -I$(dir_EH) $(dir_Python)CBL.i


#################################################################### 


$(PWD)/External/CAMB/camb:
	cd $(PWD)/External/CAMB ; make clean && make && make clean && cd ../..

$(PWD)/External/classgal_v1/class:
	cd $(PWD)/External/classgal_v1 ; make clean && make && make clean && cd ../..

CLASSpy:
	cd $(PWD)/External/classgal_v1/python/ ; python setup.py install --user

$(PWD)/External/MPTbreeze-v1/mptbreeze:
	cd External/MPTbreeze-v1/Cuba-1.4/ ; ./configure && make lib && cd ../../../
	cd $(PWD)/External/MPTbreeze-v1 ; ./compile.sh && cd ../..

$(PWD)/External/fftlog-f90-master/fftlog-f90:
	cd $(PWD)/External/fftlog-f90-master ; make clean && make "F90 = gfortran -g -w" && make clean && cd ../..


$(PWD)/External/mangle/bin/ransack:
	cd $(PWD)/External/mangle/src ; make cleaner ; ./configure && make && make clean && cd -


$(PWD)/External/VIPERS/venice3.9/venice:
	cd $(PWD)/External/VIPERS/venice3.9 ; make clean && make "CC = gcc" && make && make clean && cd -

