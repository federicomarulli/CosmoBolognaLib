C = g++
F = gfortran

FLAGS0 = -std=c++11 -fopenmp

FLAGS = -O3 -unroll -ftree-loop-distribution -Wall -Wextra -pedantic #-floop-nest-optimize -floop-parallelize-all 

dir_H = $(PWD)/Headers/Lib/
dir_O = $(PWD)/Headers/Objects/
dir_EH = $(PWD)/Cosmology/EH/

HH = $(dir_H)*.h

FLAGS_INC = -I$(dir_H) -isystem $(HOME)/Numerical/
FLAGS_FFTW = -lfftw3 -lfftw3_omp
FLAGS_GSL = -lgsl -lgslcblas -lm
FLAGS_LINK = -shared

ES = so

varDIR = -DSYS=\"Linux\"

ifeq ($(SYS),MAC)
	varDIR = -DSYS=\"MAC\" 
	FLAGS0 = -std=c++11 -fopenmp
	FLAGS_FFTW = -lfftw3 
	FLAGS_LINK = -dynamiclib -undefined suppress -flat_namespace
        ES = dylib
endif

FLAGST = $(FLAGS0) $(FLAGS)


####################################################################


dir_FUNC = $(PWD)/Func/
dir_STAT = $(PWD)/Statistics/
dir_COSM = $(PWD)/Cosmology/Lib/
dir_CM = $(PWD)/ChainMesh/
dir_CAT = $(PWD)/Catalogue/
dir_LN = $(PWD)/LogNormal/
dir_RANDOM = $(PWD)/CatalogueAnalysis/RandomCatalogue/
dir_TWOP = $(PWD)/CatalogueAnalysis/TwoPointCorrelation/
dir_MTWOP = $(PWD)/CatalogueAnalysis/ModelTwoPointCorrelation/
dir_THREEP = $(PWD)/CatalogueAnalysis/ThreePointCorrelation/
dir_GLOB = $(PWD)/GlobalFunc/

OBJ_FUNC = $(dir_FUNC)Func.o $(dir_FUNC)FuncXi.o $(dir_FUNC)FuncMultipoles.o
OBJ_STAT = $(dir_STAT)Prior/Prior.o $(dir_STAT)Chi2/Chi2.o $(dir_STAT)MCMC/MCMC.o 
OBJ_COSM = $(dir_EH)power_whu.o $(dir_COSM)Cosmology.o $(dir_COSM)Sigma.o $(dir_COSM)PkXi.o $(dir_COSM)PkXizSpace.o $(dir_COSM)Bias.o $(dir_COSM)RSD.o $(dir_COSM)Velocities.o $(dir_COSM)MassGrowth.o $(dir_COSM)NG.o $(dir_COSM)BAO.o $(dir_COSM)MassFunction.o $(dir_COSM)SizeFunction.o
OBJ_CM = $(dir_CM)ChainMesh.o
OBJ_CAT = $(dir_CAT)Catalogue.o $(dir_CAT)ChainMesh_Catalogue.o
OBJ_LN = $(dir_LN)LogNormal.o
OBJ_RANDOM = $(dir_RANDOM)RandomCatalogue.o
OBJ_TWOP = $(dir_TWOP)Pairs.o $(dir_TWOP)Init.o $(dir_TWOP)IO.o $(dir_TWOP)Measurements.o $(dir_TWOP)Errors.o $(dir_TWOP)RealSpaceCorrelations.o $(dir_TWOP)Bias.o $(dir_TWOP)Multipoles.o $(dir_TWOP)FuncTest.o
OBJ_MTWOP = $(dir_MTWOP)Init.o $(dir_MTWOP)IO.o $(dir_MTWOP)HaloHost.o $(dir_MTWOP)DispersionModelXiMeasured.o $(dir_MTWOP)DispersionModel.o
OBJ_THREEP = $(dir_THREEP)Triplets.o $(dir_THREEP)Init.o $(dir_THREEP)IO.o $(dir_THREEP)Measurements.o 
OBJ_GLOB = $(dir_GLOB)FuncCosmology.o $(dir_GLOB)Func.o $(dir_GLOB)SubSample.o

OBJ_CBL = $(OBJ_FUNC) $(OBJ_STAT) $(OBJ_COSM) $(OBJ_CM) $(OBJ_CAT) $(OBJ_LN) $(OBJ_RANDOM) $(OBJ_TWOP) $(OBJ_MTWOP) $(OBJ_THREEP) $(OBJ_GLOB)
OBJ_ALL = $(OBJ_CBL) $(PWD)/Cosmology/CAMB/*.o $(PWD)/Cosmology/classgal_v1/*.o


#################################################################### 

define colorecho
      @tput setaf 3
      @echo $1
      @tput sgr0
endef

ALL:
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
	$(call colorecho, "\n"Compiling the library: libRANDOM... "\n")
	make -j3 libRANDOM
	$(call colorecho, "\n"Compiling the library: libTWOP... "\n")
	make -j3 libTWOP
	$(call colorecho, "\n"Compiling the library: libMTWOP... "\n")
	make -j3 libMTWOP
	$(call colorecho, "\n"Compiling the library: libTHREEP... "\n")
	make -j3 libTHREEP
	$(call colorecho, "\n"Compiling the library: libGLOB... "\n")
	make -j3 libGLOB

libFUNC: $(OBJ_FUNC) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libFUNC.$(ES) $(OBJ_FUNC) $(FLAGS_GSL) -lgomp

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

libRANDOM: $(OBJ_RANDOM) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libRANDOM.$(ES) $(OBJ_RANDOM) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN

libTWOP: $(OBJ_TWOP) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libTWOP.$(ES) $(OBJ_TWOP) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lRANDOM

libMTWOP: $(OBJ_MTWOP) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libMTWOP.$(ES) $(OBJ_MTWOP) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lRANDOM -lTWOP

libTHREEP: $(OBJ_THREEP) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libTHREEP.$(ES) $(OBJ_THREEP) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lRANDOM -lTWOP -lMTWOP

libGLOB: $(OBJ_GLOB) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libGLOB.$(ES) $(OBJ_GLOB) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lRANDOM -lTWOP -lMTWOP -lTHREEP

conv: $(dir_FUNC)conv.o
	$(F) -o $(dir_FUNC)conv $(dir_FUNC)conv.o 

CAMB:
	cd $(PWD)/Cosmology/CAMB ; make clean && make && make clean && cd ../..

CLASS:
	cd $(PWD)/Cosmology/classgal_v1 ; make clean && make && make clean && cd ../..

CLASSpy:
	cd $(PWD)/Cosmology/classgal_v1/python/ ; python setup.py install --user

MPTbreeze:
	cd $(PWD)/Cosmology/MPTbreeze-v1 ; ./compile.sh && cd ../.. 

python:
	cd $(PWD)/Python/CosmoBolognaLib/ ; python setup.py install --user

clean:
	rm -f $(OBJ_ALL) core* $(PWD)/*~ $(dir_FUNC)*~ $(dir_STAT)Prior/*~ $(dir_STAT)Chi2/*~ $(dir_STAT)MCMC/*~ $(dir_COSM)*~ $(dir_CM)*~ $(dir_CAT)*~ $(dir_LN)*~ $(dir_RANDOM)*~ $(dir_TWOP)*~ $(dir_MTWOP)*~ $(dir_THREEP)*~ $(dir_GLOB)*~ $(dir_H)*~ $(dir_O)*~ $(PWD)/\#* $(dir_FUNC)\#* $(dir_STAT)\#* $(dir_COSM)\#* $(dir_CM)\#* $(dir_CAT)\#* $(dir_LN)\#* $(dir_RANDOM)\#* $(dir_TWOP)\#* $(dir_MTWOP)\#* $(dir_THREEP)\#* $(dir_GLOB)\#* $(dir_H)\#* $(dir_O)\#* $(PWD)/Doc/WARNING_LOGFILE*
	cd $(PWD)/Python/CosmoBolognaLib ; python setup.py clean
	cd $(PWD)/Python/CosmoBolognaLib/ ; rm -rf CosmoBolognaLib.cpp external.cpp build Cosmology/Lib
	cd $(PWD)/Cosmology/classgal_v1/python/ ; python setup.py clean --all

purge:
	make clean
	rm -f *.$(ES) test temp*

purgeALL:
	make purge
	rm -rf Cosmology/grid_SigmaM/* ;
	rm -rf Cosmology/Tables/* ;
	rm -rf Cosmology/table_dc/* ;
	rm -rf Cosmology/EisensteinHu/output_linear/* ;
	rm -rf Cosmology/CAMB/output_linear/* ;
	rm -rf Cosmology/CAMB/output_nonlinear/* ;
	rm -rf Cosmology/grid_NG/* ;
	rm -rf Cosmology/MPTbreeze-v1/output_linear/* ;
	rm -rf Cosmology/MPTbreeze-v1/output_nonlinear/* ;
	rm -rf Cosmology/classgal_v1/output_linear/* ;
	rm -rf Cosmology/classgal_v1/output_nonlinear/* ;


#################################################################### 


$(dir_FUNC)Func.o: $(dir_FUNC)Func.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) $(varDIR) -c -fPIC $(FLAGS_INC) -I$(dir_FUNC) $(dir_FUNC)Func.cpp -o $(dir_FUNC)Func.o

$(dir_FUNC)FuncXi.o: $(dir_FUNC)FuncXi.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_FUNC) $(dir_FUNC)FuncXi.cpp -o $(dir_FUNC)FuncXi.o

$(dir_FUNC)FuncMultipoles.o: $(dir_FUNC)FuncMultipoles.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_FUNC) $(dir_FUNC)FuncMultipoles.cpp -o $(dir_FUNC)FuncMultipoles.o

$(dir_FUNC)conv.o: $(dir_FUNC)conv.f90 
	$(F) -c $(dir_FUNC)conv.f90 -o $(dir_FUNC)conv.o


#################################################################### 


$(dir_STAT)Prior/Prior.o: $(dir_STAT)Prior/Prior.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_STAT) $(dir_STAT)Prior/Prior.cpp -o $(dir_STAT)Prior/Prior.o

$(dir_STAT)Chi2/Chi2.o: $(dir_STAT)Chi2/Chi2.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_STAT) $(dir_STAT)Chi2/Chi2.cpp -o $(dir_STAT)Chi2/Chi2.o

$(dir_STAT)MCMC/MCMC.o: $(dir_STAT)MCMC/MCMC.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_STAT) $(dir_STAT)MCMC/MCMC.cpp -o $(dir_STAT)MCMC/MCMC.o


#################################################################### 


$(dir_COSM)Cosmology.o: $(dir_COSM)Cosmology.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)Cosmology.cpp -o $(dir_COSM)Cosmology.o

$(dir_COSM)MassFunction.o: $(dir_COSM)MassFunction.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)MassFunction.cpp -o $(dir_COSM)MassFunction.o

$(dir_COSM)SizeFunction.o: $(dir_COSM)SizeFunction.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)SizeFunction.cpp -o $(dir_COSM)SizeFunction.o

$(dir_COSM)PkXi.o: $(dir_COSM)PkXi.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)PkXi.cpp -o $(dir_COSM)PkXi.o

$(dir_COSM)PkXizSpace.o: $(dir_COSM)PkXizSpace.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)PkXizSpace.cpp -o $(dir_COSM)PkXizSpace.o

$(dir_COSM)Bias.o: $(dir_COSM)Bias.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)Bias.cpp -o $(dir_COSM)Bias.o

$(dir_COSM)RSD.o: $(dir_COSM)RSD.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)RSD.cpp -o $(dir_COSM)RSD.o

$(dir_COSM)Sigma.o: $(dir_COSM)Sigma.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)Sigma.cpp -o $(dir_COSM)Sigma.o

$(dir_COSM)Velocities.o: $(dir_COSM)Velocities.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)Velocities.cpp -o $(dir_COSM)Velocities.o

$(dir_COSM)MassGrowth.o: $(dir_COSM)MassGrowth.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)MassGrowth.cpp -o $(dir_COSM)MassGrowth.o

$(dir_COSM)NG.o: $(dir_COSM)NG.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)NG.cpp -o $(dir_COSM)NG.o

$(dir_COSM)BAO.o: $(dir_COSM)BAO.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) -I$(dir_EH) $(dir_COSM)BAO.cpp -o $(dir_COSM)BAO.o

$(dir_EH)power_whu.o: $(dir_EH)power_whu.cpp $(dir_EH)power_whu.h
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_COSM) $(dir_EH)power_whu.cpp -o $(dir_EH)power_whu.o


#################################################################### 


$(dir_CM)ChainMesh.o: $(dir_CM)ChainMesh.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_CM) -I$(dir_EH) $(dir_CM)ChainMesh.cpp -o $(dir_CM)ChainMesh.o


#################################################################### 


$(dir_CAT)Catalogue.o: $(dir_CAT)Catalogue.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_CAT) -I$(dir_EH) -I$(dir_O) $(dir_CAT)Catalogue.cpp -o $(dir_CAT)Catalogue.o

$(dir_CAT)ChainMesh_Catalogue.o: $(dir_CAT)ChainMesh_Catalogue.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_CAT) -I$(dir_EH) -I$(dir_O) $(dir_CAT)ChainMesh_Catalogue.cpp -o $(dir_CAT)ChainMesh_Catalogue.o


#################################################################### 


$(dir_LN)LogNormal.o: $(dir_LN)LogNormal.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_LN) -I$(dir_EH) -I$(dir_O) $(dir_LN)LogNormal.cpp -o $(dir_LN)LogNormal.o


#################################################################### 


$(dir_RANDOM)RandomCatalogue.o: $(dir_RANDOM)RandomCatalogue.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_RANDOM) -I$(dir_EH) -I$(dir_O) $(dir_RANDOM)RandomCatalogue.cpp -o $(dir_RANDOM)RandomCatalogue.o


#################################################################### 


$(dir_TWOP)Pairs.o: $(dir_TWOP)Pairs.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)Pairs.cpp -o $(dir_TWOP)Pairs.o

$(dir_TWOP)Init.o: $(dir_TWOP)Init.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)Init.cpp -o $(dir_TWOP)Init.o

$(dir_TWOP)IO.o: $(dir_TWOP)IO.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)IO.cpp -o $(dir_TWOP)IO.o

$(dir_TWOP)Measurements.o: $(dir_TWOP)Measurements.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)Measurements.cpp -o $(dir_TWOP)Measurements.o

$(dir_TWOP)Errors.o: $(dir_TWOP)Errors.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)Errors.cpp -o $(dir_TWOP)Errors.o

$(dir_TWOP)RealSpaceCorrelations.o: $(dir_TWOP)RealSpaceCorrelations.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)RealSpaceCorrelations.cpp -o $(dir_TWOP)RealSpaceCorrelations.o

$(dir_TWOP)Bias.o: $(dir_TWOP)Bias.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)Bias.cpp -o $(dir_TWOP)Bias.o

$(dir_TWOP)Multipoles.o: $(dir_TWOP)Multipoles.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)Multipoles.cpp -o $(dir_TWOP)Multipoles.o

$(dir_TWOP)FuncTest.o: $(dir_TWOP)FuncTest.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_TWOP) -I$(dir_EH) -I$(dir_O) $(dir_TWOP)FuncTest.cpp -o $(dir_TWOP)FuncTest.o


#################################################################### 


$(dir_MTWOP)Init.o: $(dir_MTWOP)Init.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_MTWOP) -I$(dir_EH) -I$(dir_O) $(dir_MTWOP)Init.cpp -o $(dir_MTWOP)Init.o

$(dir_MTWOP)IO.o: $(dir_MTWOP)IO.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_MTWOP) -I$(dir_EH) -I$(dir_O) $(dir_MTWOP)IO.cpp -o $(dir_MTWOP)IO.o

$(dir_MTWOP)HaloHost.o: $(dir_MTWOP)HaloHost.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_MTWOP) -I$(dir_EH) -I$(dir_O) $(dir_MTWOP)HaloHost.cpp -o $(dir_MTWOP)HaloHost.o

$(dir_MTWOP)DispersionModelXiMeasured.o: $(dir_MTWOP)DispersionModelXiMeasured.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_MTWOP) -I$(dir_EH) -I$(dir_O) $(dir_MTWOP)DispersionModelXiMeasured.cpp -o $(dir_MTWOP)DispersionModelXiMeasured.o

$(dir_MTWOP)DispersionModel.o: $(dir_MTWOP)DispersionModel.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_MTWOP) -I$(dir_EH) -I$(dir_O) $(dir_MTWOP)DispersionModel.cpp -o $(dir_MTWOP)DispersionModel.o


#################################################################### 


$(dir_THREEP)Triplets.o: $(dir_THREEP)Triplets.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_THREEP) -I$(dir_EH) -I$(dir_O) $(dir_THREEP)Triplets.cpp -o $(dir_THREEP)Triplets.o

$(dir_THREEP)Init.o: $(dir_THREEP)Init.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_THREEP) -I$(dir_EH) -I$(dir_O) $(dir_THREEP)Init.cpp -o $(dir_THREEP)Init.o

$(dir_THREEP)IO.o: $(dir_THREEP)IO.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_THREEP) -I$(dir_EH) -I$(dir_O) $(dir_THREEP)IO.cpp -o $(dir_THREEP)IO.o

$(dir_THREEP)Measurements.o: $(dir_THREEP)Measurements.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_THREEP) -I$(dir_EH) -I$(dir_O) $(dir_THREEP)Measurements.cpp -o $(dir_THREEP)Measurements.o


#################################################################### 


$(dir_GLOB)FuncCosmology.o: $(dir_GLOB)FuncCosmology.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_GLOB) -I$(dir_EH) -I$(dir_O) $(dir_GLOB)FuncCosmology.cpp -o $(dir_GLOB)FuncCosmology.o

$(dir_GLOB)Func.o: $(dir_GLOB)Func.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_GLOB) -I$(dir_EH) -I$(dir_O) $(dir_GLOB)Func.cpp -o $(dir_GLOB)Func.o

$(dir_GLOB)SubSample.o: $(dir_GLOB)SubSample.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) -I$(dir_GLOB) -I$(dir_EH) -I$(dir_O) $(dir_GLOB)SubSample.cpp -o $(dir_GLOB)SubSample.o

