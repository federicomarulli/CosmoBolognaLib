C = g++
F = gfortran

FLAGS0 = -std=c++11 -fopenmp

FLAGS = -O3 -unroll -Wall -Wextra -pedantic -Wno-unused-parameter

Dir_H = Headers/Lib/
Dir_O = Headers/Objects/
Dir_M = Headers/Models/
Dir_EH = External/EH/

dir_H = $(addprefix $(PWD)/,$(Dir_H))
dir_O = $(addprefix $(PWD)/,$(Dir_O))
dir_M = $(addprefix $(PWD)/,$(Dir_M))
dir_EH = $(addprefix $(PWD)/,$(Dir_EH))
dir_PYLIB =  $(HOME)/.local/lib/python2.7/site-packages/

dir_Python = $(PWD)/Python/

HH = $(dir_H)*.h $(dir_O)*.h $(dir_M)*.h

FLAGS_INC = -I$(HOME)/include/ -I/usr/local/include/ -I$(dir_H) -I$(dir_O) -I$(dir_M) -I$(dir_EH) -isystem $(HOME)/Numerical/ 
FLAGS_FFTW = -lfftw3 -lfftw3_omp
FLAGS_GSL = -lgsl -lgslcblas -lm -L$(HOME)/lib 
FLAGS_LINK = -shared

PFLAGS = -I/usr/include/python2.7 

ES = so

Dvar = -DLINUX

ifeq ($(SYS),MAC)
	Dvar = -DMAC
	FLAGS0 = -std=c++11 -fopenmp
	FLAGS_FFTW = -lfftw3 
	FLAGS_LINK = -dynamiclib -undefined suppress -flat_namespace
        ES = dylib
	dir_PYLIB = $(HOME)/Library/Python/2.7/lib/python/site-packages/
	FLAGS_PY = -L/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config -lpython2.7 -ldl

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
Dir_TWOP = CatalogueAnalysis/TwoPointCorrelation/
Dir_THREEP = CatalogueAnalysis/ThreePointCorrelation/
Dir_MODEL = Modelling/
Dir_GLOB = GlobalFunc/

Dir_CBL = $(Dir_H) $(Dir_O) $(Dir_M) $(Dir_FUNC) $(Dir_STAT) $(Dir_COSM) $(Dir_CM) $(Dir_CAT) $(Dir_LN) $(Dir_TWOP) $(Dir_THREEP) $(Dir_MODEL) $(Dir_GLOB)
Dir_ALL = $(Dir_CBL) $(Dir_EH) Headers/ Cosmology/ CatalogueAnalysis/ 

dir_FUNC = $(addprefix $(PWD)/,$(Dir_FUNC))
dir_STAT = $(addprefix $(PWD)/,$(Dir_STAT))
dir_COSM = $(addprefix $(PWD)/,$(Dir_COSM))
dir_CM = $(addprefix $(PWD)/,$(Dir_CM))
dir_CAT = $(addprefix $(PWD)/,$(Dir_CAT))
dir_LN = $(addprefix $(PWD)/,$(Dir_LN))
dir_TWOP = $(addprefix $(PWD)/,$(Dir_TWOP))
dir_THREEP = $(addprefix $(PWD)/,$(Dir_THREEP))
dir_MODEL = $(addprefix $(PWD)/,$(Dir_MODEL))
dir_GLOB = $(addprefix $(PWD)/,$(Dir_GLOB))
dir_CBL = $(addprefix $(PWD)/,$(Dir_CBL))
dir_ALL = $(addprefix $(PWD)/,$(Dir_ALL))

dir_ALLP = $(addprefix $(dir_Python)CosmoBolognaLib/,$(Dir_ALL))


##### CBL object files #####

OBJ_FUNC = $(dir_FUNC)Func.o $(dir_FUNC)FuncXi.o $(dir_FUNC)FuncMultipoles.o $(dir_FUNC)GSLfunction.o  $(dir_FUNC)Data.o $(dir_FUNC)Data1D.o $(dir_FUNC)Data1D_collection.o $(dir_FUNC)Data2D.o $(dir_FUNC)Field3D.o
OBJ_STAT = $(dir_STAT)Chain.o $(dir_STAT)Prior.o $(dir_STAT)Parameter.o $(dir_STAT)Model.o $(dir_STAT)Chi2.o $(dir_STAT)Likelihood.o
OBJ_COSM = $(dir_EH)power_whu.o $(dir_COSM)Cosmology.o $(dir_COSM)Sigma.o $(dir_COSM)PkXi.o $(dir_COSM)PkXizSpace.o $(dir_COSM)Bias.o $(dir_COSM)RSD.o $(dir_COSM)Velocities.o $(dir_COSM)MassGrowth.o $(dir_COSM)NG.o $(dir_COSM)BAO.o $(dir_COSM)MassFunction.o $(dir_COSM)SizeFunction.o
OBJ_CM = $(dir_CM)ChainMesh.o
OBJ_CAT = $(dir_CAT)Object.o $(dir_CAT)Catalogue.o $(dir_CAT)RandomCatalogue.o $(dir_CAT)ChainMesh_Catalogue.o
OBJ_LN = $(dir_LN)LogNormal.o
OBJ_TWOP = $(dir_TWOP)Pair.o $(dir_TWOP)TwoPointCorrelation.o $(dir_TWOP)TwoPointCorrelation1D.o $(dir_TWOP)TwoPointCorrelation1D_monopole.o $(dir_TWOP)TwoPointCorrelation1D_angular.o $(dir_TWOP)TwoPointCorrelation2D.o $(dir_TWOP)TwoPointCorrelation2D_cartesian.o $(dir_TWOP)TwoPointCorrelation2D_polar.o $(dir_TWOP)TwoPointCorrelation_projected.o $(dir_TWOP)TwoPointCorrelation_deprojected.o $(dir_TWOP)TwoPointCorrelation_multipoles.o $(dir_TWOP)TwoPointCorrelation_wedges.o $(dir_TWOP)TwoPointCorrelation1D_filtered.o
OBJ_THREEP = $(dir_THREEP)Triplet.o $(dir_THREEP)ThreePointCorrelation.o $(dir_THREEP)ThreePointCorrelation_angular_connected.o $(dir_THREEP)ThreePointCorrelation_angular_reduced.o $(dir_THREEP)ThreePointCorrelation_comoving_connected.o $(dir_THREEP)ThreePointCorrelation_comoving_reduced.o 
OBJ_MODEL = $(dir_MODEL)ModelFunction.o $(dir_MODEL)ModelBias.o $(dir_MODEL)ModelBAO.o $(dir_MODEL)Modelling.o $(dir_MODEL)Modelling_TwoPointCorrelation.o $(dir_MODEL)Modelling_TwoPointCorrelation_monopole.o $(dir_MODEL)Modelling_TwoPointCorrelation_projected.o $(dir_MODEL)Modelling_TwoPointCorrelation_deprojected.o 
OBJ_GLOB = $(dir_GLOB)FuncCosmology.o $(dir_GLOB)Func.o $(dir_GLOB)SubSample.o

OBJ_CBL = $(OBJ_FUNC) $(OBJ_STAT) $(OBJ_COSM) $(OBJ_CM) $(OBJ_CAT) $(OBJ_LN) $(OBJ_TWOP) $(OBJ_THREEP) $(OBJ_MODEL) $(OBJ_GLOB)
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
	$(call colorecho, "\n"Compiling the library: libMODEL... "\n")
	make -j3 libMODEL
	$(call colorecho, "\n"Compiling the library: libGLOB... "\n")
	make -j3 libGLOB
	$(call colorecho, "\n"Compiling the full library: libCBL... "\n")
	make -j3 libCBL

libFUNC: $(OBJ_FUNC) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libFUNC.$(ES) $(OBJ_FUNC) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW)

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

libMODEL: $(OBJ_MODEL) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libMODEL.$(ES) $(OBJ_MODEL) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP

libGLOB: $(OBJ_GLOB) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libGLOB.$(ES) $(OBJ_GLOB) $(FLAGS_GSL) -lgomp -Wl,-rpath,$(PWD) -L$(PWD)/ -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lTWOP -lTHREEP -lMODEL 

libCBL: $(OBJ_CBL) $(PWD)/Makefile
	$(C) $(FLAGS_LINK) -o $(PWD)/libCBL.$(ES) $(OBJ_CBL) $(FLAGS_GSL) -lgomp $(FLAGS_FFTW)


conv: $(dir_FUNC)conv.o
	$(F) -o $(dir_FUNC)conv $(dir_FUNC)conv.o 

CAMB:
	cd $(PWD)/External/CAMB ; make clean && make && make clean && cd ../..

CLASS:
	cd $(PWD)/External/classgal_v1 ; make clean && make && make clean && cd ../..

CLASSpy:
	cd $(PWD)/External/classgal_v1/python/ ; python setup.py install --user

MPTbreeze:
	cd $(PWD)/External/MPTbreeze-v1 ; ./compile.sh && cd ../..

fftlog-f90:
	cd $(PWD)/External/fftlog-f90-master ; make clean && make && make clean && cd ../..

allExamples:
	$(call colorecho, "\n"Compiling the example code: vector.cpp ... "\n")
	cd $(PWD)/Examples/vectors ; make
	$(call colorecho, "\n"Compiling the example code: randomNumbers.cpp ... "\n")
	cd $(PWD)/Examples/randomNumbers ; make 
	$(call colorecho, "\n"Compiling the example code: distances.cpp ... "\n")
	cd $(PWD)/Examples/distances ; make 
	$(call colorecho, "\n"Compiling the example code: covsample.cpp ... "\n")
	cd $(PWD)/Examples/covsample ; make 
	$(call colorecho, "\n"Compiling the example code: fsigma8.cpp ... "\n")
	cd $(PWD)/Examples/fsigma8 ; make 
	$(call colorecho, "\n"Compiling the example code: prior.cpp ... "\n")
	cd $(PWD)/Examples/statistics/codes ; make prior
	$(call colorecho, "\n"Compiling the example code: catalogue.cpp ... "\n")
	cd $(PWD)/Examples/catalogue ; make catalogue 
	$(call colorecho, "\n"Compiling the example code: 2pt_monopole.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_monopole
	$(call colorecho, "\n"Compiling the example code: 2pt_2D.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_2D
	$(call colorecho, "\n"Compiling the example code: 2pt_jackknife.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_jackknife
	$(call colorecho, "\n"Compiling the example code: 2pt_projected_jackknife.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 2pt_projected_jackknife
	$(call colorecho, "\n"Compiling the example code: modelBias_2pt_projected.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make modelBias_2pt_projected
	$(call colorecho, "\n"Compiling the example code: 3pt.cpp ... "\n")
	cd $(PWD)/Examples/clustering/codes ; make 3pt

python: $(OBJ_CBL) $(dir_Python)CBL_wrap.o
	$(C) -shared $(OBJ_CBL) $(dir_Python)CBL_wrap.o -o $(dir_PYLIB)_CosmoBolognaLib.so $(FLAGS_GSL) $(FLAGS_FFTW) -lgomp $(FLAGS_PY)
	cp $(dir_Python)CosmoBolognaLib.py* $(dir_PYLIB)

python2: $(dir_Python)CBL_wrap.cxx $(dir_Python)CBL.i $(dir_Python)setup.py $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Compiling the python wrapper. It may take a few minutes ... "\n")
	cd $(dir_Python)CosmoBolognaLib ; mkdir -p $(Dir_ALL)  
	rsync -av $(dir_H)*.h $(dir_Python)CosmoBolognaLib/$(Dir_H)
	rsync -av $(dir_O)*.h $(dir_Python)CosmoBolognaLib/$(Dir_O)
	rsync -av $(dir_M)*.h $(dir_Python)CosmoBolognaLib/$(Dir_M)
	rsync -av $(dir_EH)* $(dir_Python)CosmoBolognaLib/$(Dir_EH)
	rsync -av $(dir_FUNC)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_FUNC)
	rsync -av $(dir_STAT)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_STAT)
	rsync -av $(dir_COSM)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_COSM)
	rsync -av $(dir_CM)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_CM)
	rsync -av $(dir_CAT)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_CAT)
	rsync -av $(dir_LN)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_LN)
	rsync -av $(dir_TWOP)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_TWOP)
	rsync -av $(dir_THREEP)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_THREEP)
	rsync -av $(dir_MODEL)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_MODEL)
	rsync -av $(dir_GLOB)*.cpp $(dir_Python)CosmoBolognaLib/$(Dir_GLOB)
	cd $(dir_Python) ; python setup.py install --user 

doc:
	rm Doc/html/* Doc/xml/* -rf
	doxygen .dconfig
	rm Doc/doxygen_sqlite3.db -f 
#	python ../bin/doxy2swig2.py Doc/xml/index.xml Doc/documentation.i
#	python ../doxy2swig/doxy2swig.py Doc/xml/index.xml Doc/documentation.i

doct:
	rm Doc/html/* Doc/xml/* -rf
	doxygen .dconfigT
	rm Doc/doxygen_sqlite3.db -f 

cleanExamples:
	cd $(PWD)/Examples/vectors ; make clean && cd ../..
	cd $(PWD)/Examples/randomNumbers ; make clean && cd ../..
	cd $(PWD)/Examples/distances ; make clean && cd ../..
	cd $(PWD)/Examples/covsample ; make clean && cd ../..
	cd $(PWD)/Examples/fsigma8 ; make clean && cd ../..
	cd $(PWD)/Examples/statistics/codes ; make clean && cd ../..
	cd $(PWD)/Examples/catalogue ; make clean && cd ../
	cd $(PWD)/Examples/clustering/codes ; make clean && cd ../..
	rm -rf $(PWD)/Examples/clustering/output/*

cleanpy:
	rm -f $(PWD)/Python/*~ $(PWD)/Python/CBL_wrap.o $(PWD)/Python/CBL_wrap.cxx $(PWD)/Python/CosmoBolognaLib.py*
	rm -rf $(PWD)/Python/dist $(PWD)/Python/build $(PWD)/Python/CosmoBolognaLib.egg-info
	rm -f $(PWD)/Python/Lib/*~ $(PWD)/Python/Lib/*.o $(PWD)/Python/Lib/*.cxx $(PWD)/Python/Lib/*.py
	rm -f $(PWD)/Python/CosmoBolognaLib/CosmoBolognaLib* $(PWD)/Python/CosmoBolognaLib/*pyc $(PWD)/Python/CosmoBolognaLib/*~
	rm -rf $(dir_ALLP) 

purgepy:
	make cleanpy
	rm -fr $(HOME)/.local/lib/python2.7/site-packages/*CosmoBolognaLib*

clean:
	make cleanExamples
	rm -f $(OBJ_ALL) core* $(PWD)/*~ $(dir_FUNC)*~ $(dir_STAT)/*~ $(dir_STAT)/*~  $(dir_STAT)/*~ $(dir_STAT)/*~ $(dir_STAT)/*~ $(dir_STAT)/*~ $(dir_STAT)/*~ $(dir_COSM)*~ $(dir_CM)*~ $(dir_CAT)*~ $(dir_LN)*~ $(dir_TWOP)*~ $(dir_MODEL)*~ $(dir_THREEP)*~ $(dir_GLOB)*~ $(dir_H)*~ $(dir_O)*~ $(dir_M)*~ $(PWD)/\#* $(dir_FUNC)\#* $(dir_STAT)\#* $(dir_COSM)\#* $(dir_CM)\#* $(dir_CAT)\#* $(dir_LN)\#* $(dir_TWOP)\#* $(dir_MODEL)\#* $(dir_THREEP)\#* $(dir_GLOB)\#* $(dir_H)\#* $(dir_O)\#* $(dir_M)\#* $(PWD)/Doc/WARNING_LOGFILE* $(PWD)/Doc/*~

purge:
	make clean
	rm -f *.$(ES) temp*

purgeALL:
	make purge
	make cleanpy
	rm -rf Cosmology/Tables/* ;
	cd External/CAMB ; make clean ; cd .. ;
	rm -rf External/CAMB/output_linear/* ;
	rm -rf External/CAMB/output_nonlinear/* ;
	cd External/classgal_v1/ ; make clean ; cd .. ;
	rm -rf External/classgal_v1/output_linear/* ;
	rm -rf External/classgal_v1/output_nonlinear/* ;
	rm -rf External/EH/*.o External/EH/*~ ;
	cd External/fftlog-f90-master/ ; make clean ; cd .. ;
	cd External/mangle/src ; make clean ; cd .. ;
	rm -rf External/MPTbreeze-v1/*~ ;
	rm -rf External/MPTbreeze-v1/output_linear/* ;
	rm -rf External/MPTbreeze-v1/output_nonlinear/* ;


#################################################################### 


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

$(dir_FUNC)Field3D.o: $(dir_FUNC)Field3D.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_FUNC)Field3D.cpp -o $(dir_FUNC)Field3D.o

$(dir_FUNC)conv.o: $(dir_FUNC)conv.f90 
	$(F) -c $(dir_FUNC)conv.f90 -o $(dir_FUNC)conv.o


#################################################################### 


$(dir_STAT)Chain.o: $(dir_STAT)Chain.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Chain.cpp -o $(dir_STAT)Chain.o

$(dir_STAT)Prior.o: $(dir_STAT)Prior.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Prior.cpp -o $(dir_STAT)Prior.o

$(dir_STAT)Parameter.o: $(dir_STAT)Parameter.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Parameter.cpp -o $(dir_STAT)Parameter.o

$(dir_STAT)Model.o: $(dir_STAT)Model.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Model.cpp -o $(dir_STAT)Model.o

$(dir_STAT)Chi2.o: $(dir_STAT)Chi2.cpp $(HH) $(PWD)/Makefile 
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_STAT)Chi2.cpp -o $(dir_STAT)Chi2.o

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

$(dir_EH)power_whu.o: $(dir_EH)power_whu.cpp $(dir_EH)power_whu.h
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_EH)power_whu.cpp -o $(dir_EH)power_whu.o


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


#################################################################### 


$(dir_LN)LogNormal.o: $(dir_LN)LogNormal.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_LN)LogNormal.cpp -o $(dir_LN)LogNormal.o


#################################################################### 


$(dir_TWOP)Pair.o: $(dir_TWOP)Pair.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_TWOP)Pair.cpp -o $(dir_TWOP)Pair.o

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


$(dir_MODEL)ModelFunction.o: $(dir_MODEL)ModelFunction.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL)ModelFunction.cpp -o $(dir_MODEL)ModelFunction.o

$(dir_MODEL)ModelBias.o: $(dir_MODEL)ModelBias.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL)ModelBias.cpp -o $(dir_MODEL)ModelBias.o

$(dir_MODEL)ModelBAO.o: $(dir_MODEL)ModelBAO.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL)ModelBAO.cpp -o $(dir_MODEL)ModelBAO.o

$(dir_MODEL)Modelling.o: $(dir_MODEL)Modelling.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL)Modelling.cpp -o $(dir_MODEL)Modelling.o

$(dir_MODEL)Modelling_TwoPointCorrelation.o: $(dir_MODEL)Modelling_TwoPointCorrelation.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL)Modelling_TwoPointCorrelation.cpp -o $(dir_MODEL)Modelling_TwoPointCorrelation.o

$(dir_MODEL)Modelling_TwoPointCorrelation_monopole.o: $(dir_MODEL)Modelling_TwoPointCorrelation_monopole.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL)Modelling_TwoPointCorrelation_monopole.cpp -o $(dir_MODEL)Modelling_TwoPointCorrelation_monopole.o

$(dir_MODEL)Modelling_TwoPointCorrelation_projected.o: $(dir_MODEL)Modelling_TwoPointCorrelation_projected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL)Modelling_TwoPointCorrelation_projected.cpp -o $(dir_MODEL)Modelling_TwoPointCorrelation_projected.o

$(dir_MODEL)Modelling_TwoPointCorrelation_deprojected.o: $(dir_MODEL)Modelling_TwoPointCorrelation_deprojected.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_MODEL)Modelling_TwoPointCorrelation_deprojected.cpp -o $(dir_MODEL)Modelling_TwoPointCorrelation_deprojected.o


#################################################################### 


$(dir_GLOB)FuncCosmology.o: $(dir_GLOB)FuncCosmology.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)FuncCosmology.cpp -o $(dir_GLOB)FuncCosmology.o

$(dir_GLOB)Func.o: $(dir_GLOB)Func.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)Func.cpp -o $(dir_GLOB)Func.o

$(dir_GLOB)SubSample.o: $(dir_GLOB)SubSample.cpp $(HH) $(PWD)/Makefile
	$(C) $(FLAGST) -c -fPIC $(FLAGS_INC) $(dir_GLOB)SubSample.cpp -o $(dir_GLOB)SubSample.o


#################################################################### 


$(dir_Python)CBL_wrap.o: $(dir_Python)CBL.i $(HH) $(PWD)/Makefile
	$(call colorecho, "\n"Compiling the python wrapper. It may take a few minutes ... "\n")
	swig -Wall -python -c++ -I$(dir_H) -I$(dir_O) -I$(dir_M) -I$(dir_EH) $(dir_Python)CBL.i
	$(C) $(FLAGST) -Wno-uninitialized $(PFLAGS) -c -fPIC $(FLAGS_INC) $(dir_Python)CBL_wrap.cxx -o $(dir_Python)CBL_wrap.o

$(dir_Python)CBL_wrap.cxx: $(dir_Python)CBL.i $(HH) $(PWD)/Makefile
	swig -python -c++ -I$(dir_H) -I$(dir_O) -I$(dir_M) -I$(dir_EH) $(dir_Python)CBL.i
	mv $(dir_Python)CosmoBolognaLib.py $(dir_Python)CosmoBolognaLib/
